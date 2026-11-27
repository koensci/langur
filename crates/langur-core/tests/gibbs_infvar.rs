//! Phase 4 validation suite for targeted collapsed Gibbs `Inf_Var` on complex (`K_n >= 2`) events:
//! - Ambiguous cross-inclusion/exclusion shared edges yield `Inf_Var > beta_posterior_var(psi, N)`
//! - Empirical 90% CI coverage calibration on simulated `K_2`/`K_3` multi-path read counts
//! - End-to-end `quantify_sample`: binary `K_1` events stay exact closed-form `O(1)` while complex
//!   `K_2` events use collapsed Gibbs `Inf_Var` and empirical `[CI_Lo, CI_Hi]`
//! - 100% bit-for-bit determinism across 1-thread vs 4-thread `rayon` pools.

use flate2::read::MultiGzDecoder;
use std::io::{BufReader, Cursor, Read, Write};
use tempfile::tempdir;
use langur_core::annotation::{load_gtf_reader, GtfOptions};
use langur_core::events::gibbs::{
    sample_complex_event_gibbs, DEFAULT_GIBBS_BURN_IN, DEFAULT_GIBBS_THIN,
};
use langur_core::events::psi::beta_posterior_var;
use langur_core::events::ReducedPath;
use langur_core::index::{build_index_from_memory, DecoyMode, IndexBuildConfig, MmapIndex};
use langur_core::io::{quantify_sample, QuantConfig};

struct SimpleRng(u64);

impl SimpleRng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / ((1u64 << 53) as f64))
    }
    fn poisson(&mut self, lambda: f64) -> u32 {
        // Exact Knuth for small lambda, normal approx for large lambda
        if lambda <= 0.0 {
            return 0;
        }
        if lambda < 30.0 {
            let l = (-lambda).exp();
            let mut k = 0u32;
            let mut p = 1.0;
            loop {
                k += 1;
                p *= self.uniform().max(1e-15);
                if p <= l {
                    return k - 1;
                }
            }
        } else {
            let u1 = self.uniform().max(1e-15);
            let u2 = self.uniform();
            let z = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
            (lambda + lambda.sqrt() * z).round().max(0.0) as u32
        }
    }
}

#[test]
fn test_gibbs_cross_group_ambiguity_inflates_inf_var_and_calibrates_90pct_ci() {
    // Complex K_2 event (2 inc paths, 2 exc paths) where a shared edge spans BOTH inclusion path 0
    // and exclusion path 2 (classic tandem cassette exon topology: E1-E2-E3-E4, E1-E2-E4, E1-E3-E4, E1-E4).
    // Path rates: theta = [0.25, 0.25, 0.25, 0.25] => true Psi = 0.50.
    // Unique edges:
    //   p0 (inc): edge 1-2 (eff_len = 0.3)
    //   p1 (inc): edge 1-3 (eff_len = 0.3)
    //   p2 (exc): edge 0-2 (eff_len = 0.3)
    //   p3 (exc): edge 0-3 (eff_len = 0.3)
    // Shared ambiguous edges:
    //   g0 = {0, 2} (edge 2-3 crossing inc & exc!, eff_len = 1.4 -> dominates read counts!)
    //   g1 = {0, 1} (edge 0-1 within inc, eff_len = 0.4)
    let theta_true = [0.25, 0.25, 0.25, 0.25];
    let true_psi = 0.50;
    let total_rate = 40.0; // ~80 expected reads per replicate

    let mut rng = SimpleRng(0xCAFE_BABE_1234_5678);
    let n_reps = 200usize;
    let mut covered = 0usize;
    let mut gibbs_vars = Vec::with_capacity(n_reps);
    let mut naive_vars = Vec::with_capacity(n_reps);
    let mut mean_psis = Vec::with_capacity(n_reps);

    for rep in 0..n_reps {
        let u0 = rng.poisson(total_rate * theta_true[0] * 0.3) as f64;
        let u1 = rng.poisson(total_rate * theta_true[1] * 0.3) as f64;
        let u2 = rng.poisson(total_rate * theta_true[2] * 0.3) as f64;
        let u3 = rng.poisson(total_rate * theta_true[3] * 0.3) as f64;
        // Cross-inc/exc shared edge {0, 2}
        let c_cross = rng.poisson(total_rate * (theta_true[0] + theta_true[2]) * 1.4) as f64;
        // Within-inc shared edge {0, 1}
        let c_inc = rng.poisson(total_rate * (theta_true[0] + theta_true[1]) * 0.4) as f64;

        let total_reads = u0 + u1 + u2 + u3 + c_cross + c_inc;
        if total_reads < 10.0 {
            continue;
        }

        // Effective lengths are the sum of compatible edge effective lengths per path:
        // p0: 0.3 + 1.4 + 0.4 = 2.1; p1: 0.3 + 0.4 = 0.7; p2: 0.3 + 1.4 = 1.7; p3: 0.3
        let inc = vec![
            ReducedPath {
                nodes: vec![0, 1, 2, 3],
                count: u0,
                eff_len: 2.1,
                psi: 0.25,
            },
            ReducedPath {
                nodes: vec![0, 1, 3],
                count: u1,
                eff_len: 0.7,
                psi: 0.25,
            },
        ];
        let exc = vec![
            ReducedPath {
                nodes: vec![0, 2, 3],
                count: u2,
                eff_len: 1.7,
                psi: 0.25,
            },
            ReducedPath {
                nodes: vec![0, 3],
                count: u3,
                eff_len: 0.3,
                psi: 0.25,
            },
        ];
        let ambig = vec![(vec![0, 2], c_cross), (vec![0, 1], c_inc)];
        let init_psi = vec![0.25, 0.25, 0.25, 0.25];

        let gs = sample_complex_event_gibbs(
            &inc,
            &exc,
            &ambig,
            &init_psi,
            None,
            250,
            DEFAULT_GIBBS_BURN_IN,
            DEFAULT_GIBBS_THIN,
            100,
            rep as u32,
            2,
        )
        .unwrap();

        if gs.ci_lo <= true_psi && true_psi <= gs.ci_hi {
            covered += 1;
        }
        gibbs_vars.push(gs.inf_var);
        naive_vars.push(beta_posterior_var(gs.mean_psi, total_reads));
        mean_psis.push(gs.mean_psi);
    }

    let n_done = gibbs_vars.len() as f64;
    let coverage = (covered as f64) / n_done;
    let mean_gibbs_var = gibbs_vars.iter().sum::<f64>() / n_done;
    let mean_naive_var = naive_vars.iter().sum::<f64>() / n_done;

    assert!(
        mean_gibbs_var > mean_naive_var * 1.35,
        "cross-group shared-edge ambiguity must inflate Gibbs Inf_Var above naive Beta var \
         (gibbs={mean_gibbs_var:.6}, naive={mean_naive_var:.6})"
    );
    assert!(
        (0.82..=0.97).contains(&coverage),
        "expected ~90% empirical CI coverage on K_2 ambiguous event, got {:.3}",
        coverage
    );
}

#[test]
fn test_end_to_end_k1_closed_form_k2_gibbs_and_thread_determinism() {
    // Build a deterministic 2-gene genome:
    // - G_SIMPLE: 3 exons (E1-E2-E3 and E1-E3) -> binary K_1 CE event on node 2
    // - G_COMPLEX: 4 exons (E1-E2-E3-E4, E1-E2-E4, E1-E3-E4, E1-E4) -> K_2 CE event on node 2
    //   sharing edge E3->E4 across inclusion & exclusion!
    let mut state: u64 = 0x517C_C1B7_2722_0A95;
    let bases = [b'A', b'C', b'G', b'T'];
    let mut chrom = vec![b'A'; 12_000];
    for b in &mut chrom {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        *b = bases[((state >> 60) & 3) as usize];
    }

    let gtf = "\
chr1\tSIM\texon\t501\t650\t.\t+\t.\tgene_id \"G_SIMPLE\"; transcript_id \"TS1\";
chr1\tSIM\texon\t901\t1050\t.\t+\t.\tgene_id \"G_SIMPLE\"; transcript_id \"TS1\";
chr1\tSIM\texon\t1301\t1450\t.\t+\t.\tgene_id \"G_SIMPLE\"; transcript_id \"TS1\";
chr1\tSIM\texon\t501\t650\t.\t+\t.\tgene_id \"G_SIMPLE\"; transcript_id \"TS2\";
chr1\tSIM\texon\t1301\t1450\t.\t+\t.\tgene_id \"G_SIMPLE\"; transcript_id \"TS2\";
chr1\tSIM\texon\t4001\t4150\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC1\";
chr1\tSIM\texon\t4401\t4550\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC1\";
chr1\tSIM\texon\t4801\t4950\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC1\";
chr1\tSIM\texon\t5201\t5350\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC1\";
chr1\tSIM\texon\t4001\t4150\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC2\";
chr1\tSIM\texon\t4401\t4550\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC2\";
chr1\tSIM\texon\t5201\t5350\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC2\";
chr1\tSIM\texon\t4001\t4150\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC3\";
chr1\tSIM\texon\t4801\t4950\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC3\";
chr1\tSIM\texon\t5201\t5350\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC3\";
chr1\tSIM\texon\t4001\t4150\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC4\";
chr1\tSIM\texon\t5201\t5350\t.\t+\t.\tgene_id \"G_COMPLEX\"; transcript_id \"TC4\";
";

    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("gibbs_test.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 17,
        minimizer_w: 8,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chrom)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    // Generate 100 bp spliced reads covering all junctions of G_SIMPLE and G_COMPLEX
    let exon = |s1: usize, e1: usize| -> &[u8] { &chrom[(s1 - 1)..e1] };
    let splice_read = |left: &[u8], right: &[u8], off: usize| -> Vec<u8> {
        let mut r = Vec::with_capacity(100);
        let l_take = 35 + (off % 30);
        let r_take = 100 - l_take;
        r.extend_from_slice(&left[(left.len() - l_take)..]);
        r.extend_from_slice(&right[..r_take]);
        r
    };

    let mut reads = Vec::new();
    // G_SIMPLE: 30 E1->E2, 30 E2->E3, 30 E1->E3
    let (s1, s2, s3) = (exon(501, 650), exon(901, 1050), exon(1301, 1450));
    for i in 0..30 {
        reads.push(splice_read(s1, s2, i));
        reads.push(splice_read(s2, s3, i));
        reads.push(splice_read(s1, s3, i));
    }
    // G_COMPLEX: isoform-coherent reads from all 4 isoforms (TC1=10, TC2=10, TC3=10, TC4=10):
    //   TC1 (E1-E2-E3-E4): E1->E2, E2->E3, E3->E4
    //   TC2 (E1-E2-E4):    E1->E2, E2->E4
    //   TC3 (E1-E3-E4):    E1->E3, E3->E4
    //   TC4 (E1-E4):       E1->E4
    // Note E3->E4 (20 reads) is shared across inclusion (TC1) & exclusion (TC3), and E1->E2 (20 reads)
    // is shared across inclusion paths TC1 & TC2!
    let (c1, c2, c3, c4) = (
        exon(4001, 4150),
        exon(4401, 4550),
        exon(4801, 4950),
        exon(5201, 5350),
    );
    for i in 0..10 {
        // TC1
        reads.push(splice_read(c1, c2, i));
        reads.push(splice_read(c2, c3, i));
        reads.push(splice_read(c3, c4, i));
        // TC2
        reads.push(splice_read(c1, c2, i + 10));
        reads.push(splice_read(c2, c4, i));
        // TC3
        reads.push(splice_read(c1, c3, i));
        reads.push(splice_read(c3, c4, i + 10));
        // TC4
        reads.push(splice_read(c1, c4, i));
    }

    let fq_path = dir.path().join("reads.fq");
    {
        let mut f = std::fs::File::create(&fq_path).unwrap();
        let qual = "I".repeat(100);
        for (i, r) in reads.iter().enumerate() {
            writeln!(f, "@r{i}\n{}\n+\n{qual}", std::str::from_utf8(r).unwrap()).unwrap();
        }
    }

    let read_psi_gz = |p: &std::path::Path| -> String {
        let f = std::fs::File::open(p).unwrap();
        let mut dec = MultiGzDecoder::new(f);
        let mut s = String::new();
        dec.read_to_string(&mut s).unwrap();
        s
    };

    // Run 1: Default Gibbs (100 samples) on 1-thread rayon pool
    let out_1t = dir.path().join("out_1t");
    let mut cfg = QuantConfig {
        r1: fq_path.clone(),
        out_prefix: out_1t.clone(),
        ..Default::default()
    };
    cfg.events.gibbs_samples = 100;
    cfg.events.gibbs_seed = 42;
    cfg.chunk_size = 32;

    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| quantify_sample(&idx, &cfg).unwrap());
    let psi_1t = read_psi_gz(&dir.path().join("out_1t.psi.gz"));

    // Run 2: Default Gibbs (100 samples) on 4-thread rayon pool -> MUST be byte-for-byte identical!
    let out_4t = dir.path().join("out_4t");
    let mut cfg_4t = cfg.clone();
    cfg_4t.out_prefix = out_4t.clone();
    rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap()
        .install(|| quantify_sample(&idx, &cfg_4t).unwrap());
    let psi_4t = read_psi_gz(&dir.path().join("out_4t.psi.gz"));
    assert_eq!(
        psi_1t, psi_4t,
        "Gibbs Inf_Var output must be 100% deterministic across 1 vs 4 threads"
    );

    // Run 3: --gibbs-samples 0 (closed-form fallback)
    let out_nogibbs = dir.path().join("out_nogibbs");
    let mut cfg_ng = cfg.clone();
    cfg_ng.out_prefix = out_nogibbs.clone();
    cfg_ng.events.gibbs_samples = 0;
    quantify_sample(&idx, &cfg_ng).unwrap();
    let psi_ng = read_psi_gz(&dir.path().join("out_nogibbs.psi.gz"));

    let find_ce = |txt: &str, gene: &str, node: &str| -> (String, f64, f64, f64, String) {
        for line in txt.lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            if f[0] == gene && f[1] == node && f[4] == "CE" {
                let psi = f[5].parse::<f64>().unwrap();
                let reads = f[8].parse::<f64>().unwrap();
                let iv = f[14].parse::<f64>().unwrap();
                return (f[9].to_string(), psi, reads, iv, f[7].to_string());
            }
        }
        panic!("missing CE row for {gene} node {node}");
    };

    // Binary K_1 event (G_SIMPLE node 2): MUST stay exact closed-form O(1) (identical with or without Gibbs!)
    let (k1_comp, k1_psi, k1_reads, k1_iv_g, k1_ci_g) = find_ce(&psi_1t, "G_SIMPLE", "2");
    let (_, _, _, k1_iv_ng, k1_ci_ng) = find_ce(&psi_ng, "G_SIMPLE", "2");
    assert_eq!(k1_comp, "K1");
    assert_eq!(k1_ci_g, k1_ci_ng);
    assert!((k1_iv_g - k1_iv_ng).abs() < 1e-12);
    let k1_beta = beta_posterior_var(k1_psi, k1_reads);
    assert!(
        (k1_iv_g / k1_beta - 1.0).abs() < 1e-4,
        "K_1 event should use exact beta_posterior_var"
    );

    // Complex K_2 event (G_COMPLEX node 2): uses collapsed Gibbs Inf_Var > naive beta_posterior_var!
    let (k2_comp, k2_psi, k2_reads, k2_iv_g, k2_ci_g) = find_ce(&psi_1t, "G_COMPLEX", "2");
    let (_, _, _, k2_iv_ng, k2_ci_ng) = find_ce(&psi_ng, "G_COMPLEX", "2");
    assert_eq!(k2_comp, "K2");
    let k2_beta = beta_posterior_var(k2_psi, k2_reads);
    assert!(
        k2_iv_g > k2_beta,
        "K_2 complex event Gibbs Inf_Var ({k2_iv_g:.6e}) should exceed naive beta_posterior_var ({k2_beta:.6e})"
    );
    assert!(
        k2_ci_g != k2_ci_ng || (k2_iv_g - k2_iv_ng).abs() > 1e-7,
        "K_2 event with --gibbs-samples 100 should use empirical Gibbs CI/Inf_Var rather than heuristic"
    );
}
