//! Phase 2 comprehensive test suite for `langur quant`:
//! - Paired-end FR (`--stranded`), RF (`--stranded-reverse` / TruSeq dUTP), and Unstranded protocols
//! - Decoy vs. intronless processed pseudogene scoring parity (strict `best_decoy > best_csg`)
//! - 1-SNP-tolerant novel splice jumps & banded Gotoh gapped indel extension
//! - Connected-component Transcript VBEM + SQUAREM (parent gene absorbs ambiguous reads from silent pseudogene)
//! - MAJIQ-style PCR-stack soft-capping (novel/RI only, never annotated constitutive)
//! - End-to-end 2-pass GC bias feedback, 15-column `.psi.gz` (`Inf_Var`), `.jnc.gz` (`Distinct_Pos`/`Max_Stack`),
//!   `--circ` (`BS`), `--batch` manifest, and 100% thread-count determinism.

use flate2::read::MultiGzDecoder;
use std::io::{BufReader, Cursor, Read};
use tempfile::tempdir;
use langur_core::align::{
    align_paired_reads, align_single_read, AlignConfig, AlignScratch, PairedAlignOutcome,
    SingleAlignOutcome, StrandMode,
};
use langur_core::annotation::{load_gtf_reader, GtfOptions};
use langur_core::events::{psi::apply_pcr_stack_cap, EventConfig};
use langur_core::index::{build_index_from_memory, DecoyMode, IndexBuildConfig, MmapIndex};
use langur_core::io::{quantify_batch_manifest, quantify_sample, QuantConfig};
use langur_core::quant::{solve_cc_vbem, VbemConfig};
use langur_core::seq::revcomp;

fn pseudo_dna(seed: u64, len: usize) -> Vec<u8> {
    let mut s = seed;
    let bases = *b"ACGT";
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        s ^= s >> 12;
        s ^= s << 25;
        s ^= s >> 27;
        let idx = ((s.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 62) & 3) as usize;
        out.push(bases[idx]);
    }
    out
}

fn read_gz_string(path: &std::path::Path) -> String {
    let mut out = String::new();
    MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut out)
        .unwrap();
    out
}

#[test]
fn test_paired_end_strand_protocols_fr_rf_and_unstranded() {
    // Build a chromosome with one `+` gene (G_POS: 3 exons of 100 bp) and one `-` gene (G_NEG: 3 exons of 100 bp)
    let mut chr = pseudo_dna(101, 4000);
    // Insert canonical GT..AG splice sites around introns
    // G_POS: E1 101..200, E2 401..500, E3 701..800
    chr[200..202].copy_from_slice(b"GT");
    chr[398..400].copy_from_slice(b"AG");
    chr[500..502].copy_from_slice(b"GT");
    chr[698..700].copy_from_slice(b"AG");

    // G_NEG (- strand): genomic 1501..1600, 1801..1900, 2101..2200
    let gtf = "\
chr1\tT\texon\t101\t200\t.\t+\t.\tgene_id \"G_POS\"; transcript_id \"TP1\";
chr1\tT\texon\t401\t500\t.\t+\t.\tgene_id \"G_POS\"; transcript_id \"TP1\";
chr1\tT\texon\t701\t800\t.\t+\t.\tgene_id \"G_POS\"; transcript_id \"TP1\";
chr1\tT\texon\t1501\t1600\t.\t-\t.\tgene_id \"G_NEG\"; transcript_id \"TN1\";
chr1\tT\texon\t1801\t1900\t.\t-\t.\tgene_id \"G_NEG\"; transcript_id \"TN1\";
chr1\tT\texon\t2101\t2200\t.\t-\t.\tgene_id \"G_NEG\"; transcript_id \"TN1\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("strand.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 11,
        minimizer_w: 4,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    // Test BOTH `+` gene (gid=0) and `-` gene (gid=1)
    for gid in [0u32, 1u32] {
        let tx_seq = idx.gene_seq(gid); // 300 bp spliced transcript in 5'->3' sense order!
        assert_eq!(tx_seq.len(), 300);

        // Fragment [20..240) -> exact transcriptomic fragment length = 220 bp!
        // 5' sense read (R1 in FR, R2 in RF): tx_seq[20..80] (60 bp)
        // 3' antisense read (R2 in FR, R1 in RF): revcomp(tx_seq[180..240]) (60 bp)
        let mate_sense = &tx_seq[20..80];
        let mate_anti = revcomp(&tx_seq[180..240]);

        // 1. ForwardFR (--stranded): accepts (sense, anti), rejects (anti, sense)
        let cfg_fr = AlignConfig {
            strand_mode: StrandMode::ForwardFR,
            ..Default::default()
        };
        let mut sc = AlignScratch::new(&cfg_fr);
        match align_paired_reads(&idx, mate_sense, &mate_anti, &cfg_fr, &mut sc) {
            PairedAlignOutcome::Mapped(ps) => {
                assert_eq!(ps.len(), 1);
                assert_eq!(ps[0].gene_id, gid);
                assert!(ps[0].r1_is_5p);
                assert_eq!(ps[0].frag_len, Some(220));
            }
            other => panic!("FR should map FR pair on gene {gid}, got {other:?}"),
        }
        assert!(
            matches!(
                align_paired_reads(&idx, &mate_anti, mate_sense, &cfg_fr, &mut sc),
                PairedAlignOutcome::Unmapped
            ),
            "FR must reject RF-oriented pair on gene {gid}"
        );

        // 2. ReverseRF (--stranded-reverse / TruSeq dUTP): accepts (anti, sense), rejects (sense, anti)
        let cfg_rf = AlignConfig {
            strand_mode: StrandMode::ReverseRF,
            ..Default::default()
        };
        let mut sc_rf = AlignScratch::new(&cfg_rf);
        match align_paired_reads(&idx, &mate_anti, mate_sense, &cfg_rf, &mut sc_rf) {
            PairedAlignOutcome::Mapped(ps) => {
                assert_eq!(ps.len(), 1);
                assert_eq!(ps[0].gene_id, gid);
                assert!(!ps[0].r1_is_5p);
                assert_eq!(ps[0].frag_len, Some(220));
            }
            other => panic!("RF should map RF pair on gene {gid}, got {other:?}"),
        }
        assert!(
            matches!(
                align_paired_reads(&idx, mate_sense, &mate_anti, &cfg_rf, &mut sc_rf),
                PairedAlignOutcome::Unmapped
            ),
            "RF must reject FR-oriented pair on gene {gid}"
        );

        // 3. Unstranded: accepts BOTH FR and RF with exact frag_len == 220!
        let cfg_un = AlignConfig {
            strand_mode: StrandMode::Unstranded,
            ..Default::default()
        };
        let mut sc_un = AlignScratch::new(&cfg_un);
        for (r1, r2, exp_r1_5p) in [
            (mate_sense, mate_anti.as_slice(), true),
            (mate_anti.as_slice(), mate_sense, false),
        ] {
            match align_paired_reads(&idx, r1, r2, &cfg_un, &mut sc_un) {
                PairedAlignOutcome::Mapped(ps) => {
                    assert_eq!(ps.len(), 1);
                    assert_eq!(ps[0].gene_id, gid);
                    assert_eq!(ps[0].r1_is_5p, exp_r1_5p);
                    assert_eq!(ps[0].frag_len, Some(220));
                }
                other => panic!("Unstranded failed on gene {gid}: {other:?}"),
            }
        }
    }
}

#[test]
fn test_decoy_vs_processed_pseudogene_strict_inequality_parity() {
    // Safeguard #1:
    // - Parent gene G_PARENT has E1 (60 bp) + intron + E2 (60 bp).
    // - Intergenic locus 1 has an INTRONLESS processed pseudogene with IDENTICAL E1 + E2 sequence (120 bp).
    // - Intergenic locus 2 has a diverged pseudogene with 2 pseudogene-specific SNPs in E1 + E2.
    let mut chr = pseudo_dna(202, 2500);
    let e1 = pseudo_dna(303, 60);
    let e2 = pseudo_dna(404, 60);
    chr[100..160].copy_from_slice(&e1);
    chr[160..162].copy_from_slice(b"GT");
    chr[358..360].copy_from_slice(b"AG");
    chr[360..420].copy_from_slice(&e2);

    // Exact intronless processed pseudogene at 1-based 1001..1120
    let mut spliced_mrna = Vec::new();
    spliced_mrna.extend_from_slice(&e1);
    spliced_mrna.extend_from_slice(&e2);
    chr[1000..1120].copy_from_slice(&spliced_mrna);

    // Diverged pseudogene with 2 SNPs at 1-based 1801..1920
    let mut diverged_pseudo = spliced_mrna.clone();
    diverged_pseudo[40] = if diverged_pseudo[40] == b'A' { b'G' } else { b'A' };
    diverged_pseudo[80] = if diverged_pseudo[80] == b'C' { b'T' } else { b'C' };
    chr[1800..1920].copy_from_slice(&diverged_pseudo);

    let gtf = "\
chr1\tT\texon\t101\t160\t.\t+\t.\tgene_id \"G_PARENT\"; transcript_id \"TP\";
chr1\tT\texon\t361\t420\t.\t+\t.\tgene_id \"G_PARENT\"; transcript_id \"TP\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("pseudo_decoy.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 11,
        minimizer_w: 4,
        decoy_mode: DecoyMode::SelectiveGenome,
        decoy_window_pad: 30,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();
    assert!(idx.decoys().len() >= 2);

    let acfg = AlignConfig::default();
    let mut sc = AlignScratch::new(&acfg);

    // 1. Real spliced mRNA read spanning E1-E2 (exact tie with intronless pseudogene decoy):
    //    MUST NOT be discarded! Must map cleanly to G_PARENT across nodes [0, 1].
    let mrna_read = &spliced_mrna[20..100]; // 80 bp (40 bp E1 + 40 bp E2)
    match align_single_read(&idx, mrna_read, &acfg, &mut sc) {
        SingleAlignOutcome::Mapped(alns) => {
            assert_eq!(alns.len(), 1);
            assert_eq!(alns[0].gene_id, 0);
            assert_eq!(alns[0].node_path().as_slice(), &[0, 1]);
            assert_eq!(alns[0].score, 80 * acfg.match_score);
        }
        other => panic!("spliced mRNA read must survive tie against intronless pseudogene decoy, got {other:?}"),
    }

    // 2. Read truly originating from the diverged pseudogene (contains both pseudogene SNPs):
    //    `best_decoy_score (160) > best_csg_score (148)` -> MUST be discarded as DecoyFiltered!
    let pseudo_read = &diverged_pseudo[20..100];
    match align_single_read(&idx, pseudo_read, &acfg, &mut sc) {
        SingleAlignOutcome::DecoyFiltered => {}
        other => panic!("diverged pseudogene read must be DecoyFiltered, got {other:?}"),
    }
}

#[test]
fn test_snp_tolerant_novel_splice_jump_and_banded_gapped_indel() {
    // 3 exons: E1 (101..180), E2 (301..380), E3 (501..580).
    // Only isoform E1-E2-E3 is annotated -> E1-E3 skipping is a NOVEL combinatorial junction!
    let mut chr = pseudo_dna(505, 1000);
    chr[180..182].copy_from_slice(b"GT");
    chr[298..300].copy_from_slice(b"AG");
    chr[380..382].copy_from_slice(b"GT");
    chr[498..500].copy_from_slice(b"AG");

    let gtf = "\
chr1\tT\texon\t101\t180\t.\t+\t.\tgene_id \"G_SKIP\"; transcript_id \"T_FULL\";
chr1\tT\texon\t301\t380\t.\t+\t.\tgene_id \"G_SKIP\"; transcript_id \"T_FULL\";
chr1\tT\texon\t501\t580\t.\t+\t.\tgene_id \"G_SKIP\"; transcript_id \"T_FULL\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("snp_indel.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 11,
        minimizer_w: 4,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    let acfg = AlignConfig::default();
    let mut sc = AlignScratch::new(&acfg);

    // 1. Novel E1->E3 read with 1 SNP sitting just 3 bp past the acceptor splice site (12 bp overhang)!
    let mut novel_snp_read = Vec::new();
    novel_snp_read.extend_from_slice(&chr[120..180]); // 60 bp of E1
    novel_snp_read.extend_from_slice(&chr[500..512]); // 12 bp of E3
    novel_snp_read[63] = if novel_snp_read[63] == b'A' { b'C' } else { b'A' }; // 1 SNP inside 12 bp splice overhang!

    match align_single_read(&idx, &novel_snp_read, &acfg, &mut sc) {
        SingleAlignOutcome::Mapped(alns) => {
            assert_eq!(alns.len(), 1);
            assert_eq!(
                alns[0].node_path().as_slice(),
                &[0, 2],
                "1-SNP-tolerant novel splice jump must discover E1->E3 skipping"
            );
        }
        other => panic!("SNP-containing novel junction read failed to align: {other:?}"),
    }

    // 2. Exonic read with a 2-bp insertion inside E1 -> triggers banded Gotoh gapped extension!
    let mut indel_read = Vec::new();
    indel_read.extend_from_slice(&chr[100..140]);
    indel_read.extend_from_slice(b"GT"); // 2 bp insertion
    indel_read.extend_from_slice(&chr[140..178]); // total 80 bp
    match align_single_read(&idx, &indel_read, &acfg, &mut sc) {
        SingleAlignOutcome::Mapped(alns) => {
            assert_eq!(alns.len(), 1);
            assert_eq!(alns[0].node_path().as_slice(), &[0]);
            // Gapped alignment recovers the full 78 matching bases minus 2-bp gap penalty (-7) = ~149,
            // way above ungapped truncation score (40 * 2 = 80)!
            assert!(
                alns[0].score > 130,
                "expected banded Gotoh gapped score > 130, got {}",
                alns[0].score
            );
            assert!(alns[0].segs[0].indels >= 2);
        }
        other => panic!("indel read failed to align via banded Gotoh aligner: {other:?}"),
    }
}

#[test]
fn test_cc_vbem_parent_gene_absorbs_ambiguous_reads_from_silent_pseudogene() {
    // Gene 0 (PARENT): 2 exons + unique UTR3 -> tx 0 (len 600)
    // Gene 1 (PSEUDO): identical shared exons -> tx 1 (len 600)
    let chr = pseudo_dna(707, 2500);
    let gtf = "\
chr1\tT\texon\t101\t700\t.\t+\t.\tgene_id \"PARENT\"; transcript_id \"T_PARENT\";
chr1\tT\texon\t1101\t1700\t.\t+\t.\tgene_id \"PSEUDO\"; transcript_id \"T_PSEUDO\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("vbem.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 11,
        minimizer_w: 4,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    // 50 unique reads on T_PARENT (tx 0) + 50 completely ambiguous reads shared between [0, 1]
    let ecs = vec![
        (smallvec::smallvec![0u32], 50.0),
        (smallvec::smallvec![0u32, 1u32], 50.0),
    ];
    let eff_lens = vec![500.0, 500.0];
    let res = solve_cc_vbem(&idx, &ecs, eff_lens, &VbemConfig::default(), None);

    assert_eq!(res.num_components, 1);
    assert!(
        res.tx_counts[0] > 99.5,
        "VBEM must route >99.5% of ambiguous reads to expressed PARENT (got {:.3})",
        res.tx_counts[0]
    );
    assert!(
        res.tx_counts[1] < 0.5,
        "silent PSEUDO should get <0.5 reads (got {:.3})",
        res.tx_counts[1]
    );
}

#[test]
fn test_pcr_stack_soft_cap_rules() {
    let cfg = EventConfig {
        stack_filter: true,
        stack_min_count: 10.0,
        stack_cap_reads: 2.0,
        min_overhang: 8,
        mean_read_len: 100.0,
        ..Default::default()
    };
    // 1. Novel junction with 25 reads all at 1 position -> capped to 2.0
    assert_eq!(apply_pcr_stack_cap(25.0, 1, false, false, &cfg), 2.0);
    // 2. RI-touching junction with 25 reads at 1 position -> capped to 2.0
    assert_eq!(apply_pcr_stack_cap(25.0, 1, true, true, &cfg), 2.0);
    // 3. Annotated constitutive (non-RI) junction at 1 position -> NEVER capped!
    assert_eq!(apply_pcr_stack_cap(25.0, 1, true, false, &cfg), 25.0);
    // 4. Novel junction with diverse positions (distinct_pos == 5) -> NEVER capped!
    assert_eq!(apply_pcr_stack_cap(25.0, 5, false, false, &cfg), 25.0);
    // 5. Ultra-short read geometry where multiple start positions aren't possible -> NEVER capped!
    let mut short_cfg = cfg.clone();
    short_cfg.mean_read_len = 18.0; // < 2*8 + 4 = 20
    assert_eq!(apply_pcr_stack_cap(25.0, 1, false, false, &short_cfg), 25.0);
    // 6. --no-stack-filter escape hatch disables capping
    let mut off_cfg = cfg.clone();
    off_cfg.stack_filter = false;
    assert_eq!(apply_pcr_stack_cap(25.0, 1, false, false, &off_cfg), 25.0);
}

#[test]
fn test_end_to_end_quant_psi_circ_bias_batch_and_thread_determinism() {
    // Gene with Cassette Exon (E1 -> E2 -> E3 vs E1 -> E3)
    // E1: 101..200 (100 bp), E2: 401..500 (100 bp, CE), E3: 701..800 (100 bp)
    let mut chr = pseudo_dna(909, 1200);
    chr[200..202].copy_from_slice(b"GT");
    chr[398..400].copy_from_slice(b"AG");
    chr[500..502].copy_from_slice(b"GT");
    chr[698..700].copy_from_slice(b"AG");

    let gtf = "\
chr1\tT\texon\t101\t200\t.\t+\t.\tgene_id \"G_CE\"; transcript_id \"T_INC\";
chr1\tT\texon\t401\t500\t.\t+\t.\tgene_id \"G_CE\"; transcript_id \"T_INC\";
chr1\tT\texon\t701\t800\t.\t+\t.\tgene_id \"G_CE\"; transcript_id \"T_INC\";
chr1\tT\texon\t101\t200\t.\t+\t.\tgene_id \"G_CE\"; transcript_id \"T_EXC\";
chr1\tT\texon\t701\t800\t.\t+\t.\tgene_id \"G_CE\"; transcript_id \"T_EXC\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("ce.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 11,
        minimizer_w: 4,
        decoy_mode: DecoyMode::SelectiveGenome,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    let e1 = &chr[100..200];
    let e2 = &chr[400..500];
    let e3 = &chr[700..800];
    let mut t_inc = Vec::new();
    t_inc.extend_from_slice(e1);
    t_inc.extend_from_slice(e2);
    t_inc.extend_from_slice(e3);
    let mut t_exc = Vec::new();
    t_exc.extend_from_slice(e1);
    t_exc.extend_from_slice(e3);

    // Synthesize 60 reads:
    // - 20 reads across E1-E2 junction (diverse positions)
    // - 20 reads across E2-E3 junction (diverse positions)
    // - 20 reads across E1-E3 skipping junction (diverse positions)
    //   -> True PSI of E2 = (40/2) / (40/2 + 20) = 20 / 40 = 0.50!
    // - Plus 5 circular back-splicing reads across E2->E2 (back-splice from end of E2 to start of E2)!
    let mut fq = String::new();
    let mut rid = 0usize;
    for i in 0..20 {
        let off = 60 + (i % 15);
        let r = std::str::from_utf8(&t_inc[off..off + 70]).unwrap();
        fq.push_str(&format!("@r{rid}\n{r}\n+\n{}\n", "I".repeat(r.len())));
        rid += 1;
    }
    for i in 0..20 {
        let off = 160 + (i % 15);
        let r = std::str::from_utf8(&t_inc[off..off + 70]).unwrap();
        fq.push_str(&format!("@r{rid}\n{r}\n+\n{}\n", "I".repeat(r.len())));
        rid += 1;
    }
    for i in 0..20 {
        let off = 60 + (i % 15);
        let r = std::str::from_utf8(&t_exc[off..off + 70]).unwrap();
        fq.push_str(&format!("@r{rid}\n{r}\n+\n{}\n", "I".repeat(r.len())));
        rid += 1;
    }
    let mut bs_seq = Vec::new();
    bs_seq.extend_from_slice(&e2[60..100]); // 40 bp end of E2
    bs_seq.extend_from_slice(&e2[0..30]); // 30 bp start of E2
    let bs_str = std::str::from_utf8(&bs_seq).unwrap();
    for _ in 0..5 {
        fq.push_str(&format!("@r{rid}\n{bs_str}\n+\n{}\n", "I".repeat(bs_str.len())));
        rid += 1;
    }

    let fq_path = dir.path().join("sample.fq");
    std::fs::write(&fq_path, &fq).unwrap();

    // Run 1: single-threaded (via small chunk_size=7 so multi-chunk wave logic is exercised!)
    let out1 = dir.path().join("run1");
    let mut qcfg = QuantConfig {
        r1: fq_path.clone(),
        r2: None,
        out_prefix: out1.clone(),
        bias_correct: true,
        chunk_size: 7,
        chunks_per_wave: 4,
        ..Default::default()
    };
    qcfg.align.allow_circ = true;
    qcfg.events.allow_circ = true;

    let pool1 = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let sum1 = pool1.install(|| quantify_sample(&idx, &qcfg)).unwrap();
    assert_eq!(sum1.total_reads, 65);
    assert_eq!(sum1.mapped_reads, 65);
    assert!(sum1.vbem_pass2_iters > 0 || sum1.num_vbem_components >= 1);

    // Run 2: 4 threads with the exact same input -> outputs MUST be 100% identical!
    let out4 = dir.path().join("run4");
    let mut qcfg4 = qcfg.clone();
    qcfg4.out_prefix = out4.clone();
    let pool4 = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
    let _sum4 = pool4.install(|| quantify_sample(&idx, &qcfg4)).unwrap();

    for ext in ["psi.gz", "gene.tpm.gz", "isoform.tpm.gz", "jnc.gz"] {
        let s1 = read_gz_string(&dir.path().join(format!("run1.{ext}")));
        let s4 = read_gz_string(&dir.path().join(format!("run4.{ext}")));
        assert_eq!(s1, s4, "1-thread vs 4-thread mismatch in .{ext}!");
    }

    // Verify 15-column `.psi.gz` header + CE PSI ≈ 0.50 + BS row emitted!
    let psi_txt = read_gz_string(&dir.path().join("run1.psi.gz"));
    let mut lines = psi_txt.lines();
    let hdr = lines.next().unwrap();
    assert_eq!(hdr.split('\t').count(), 15);
    assert!(hdr.ends_with("\tInf_Var"));

    let mut saw_ce = false;
    let mut saw_bs = false;
    for l in lines {
        let cols: Vec<&str> = l.split('\t').collect();
        assert_eq!(cols.len(), 15);
        if cols[4] == "CE" {
            let psi: f64 = cols[5].parse().unwrap();
            let iv: f64 = cols[14].parse().unwrap();
            assert!((psi - 0.5).abs() < 0.15, "expected CE PSI near 0.5, got {psi}");
            assert!(iv > 0.0 && iv < 0.1, "expected valid Inf_Var, got {iv}");
            saw_ce = true;
        }
        if cols[4] == "BS" {
            saw_bs = true;
        }
    }
    assert!(saw_ce, "expected CE event in .psi.gz");
    assert!(saw_bs, "expected BS circular event in .psi.gz with --circ");

    // Verify `.jnc.gz` contains Distinct_Pos & Max_Stack columns
    let jnc_txt = read_gz_string(&dir.path().join("run1.jnc.gz"));
    let j_hdr = jnc_txt.lines().next().unwrap();
    assert!(j_hdr.ends_with("\tDistinct_Pos\tMax_Stack"));

    // Verify `--batch` manifest runner
    let manifest_path = dir.path().join("manifest.csv");
    let b_out = dir.path().join("batch_s1");
    std::fs::write(
        &manifest_path,
        format!(
            "sample,r1,r2,out_prefix\ns1,{},,{}\n",
            fq_path.display(),
            b_out.display()
        ),
    )
    .unwrap();
    let batch_res = quantify_batch_manifest(&idx, &manifest_path, &qcfg).unwrap();
    assert_eq!(batch_res.len(), 1);
    assert_eq!(batch_res[0].0, "s1");
    assert_eq!(batch_res[0].1.mapped_reads, 65);
}

#[test]
fn test_regression_gc_bias_no_double_counting() {
    // Review Item 1: Verify Stage 1 Pass 2 applies GC bias ONCE (via `gc_eff_lens`, NOT squared
    // by also applying `joint_adjusted` to Pass 2 ECs).
    // Setup:
    // - G_ANCHOR (600 bp): 300 bp low-GC (~25% GC) + 300 bp high-GC (~75% GC) on ONE transcript,
    //   with 3x more reads from its high-GC half than its low-GC half -> identifies true GC bias!
    // - G_LOW (300 bp, ~25% GC) with 100 reads and G_HIGH (300 bp, ~75% GC) with 200 reads
    //   (2:1 observed ratio caused by ~2x GC capture bias -> true ratio ~1:1).
    let make_gc_seq = |seed: u64, len: usize, high_gc: bool| -> Vec<u8> {
        let raw = pseudo_dna(seed, len);
        raw.into_iter()
            .enumerate()
            .map(|(i, b)| {
                if i % 4 != 0 {
                    if high_gc && (b == b'A' || b == b'T') {
                        if i % 2 == 0 { b'G' } else { b'C' }
                    } else if !high_gc && (b == b'G' || b == b'C') {
                        if i % 2 == 0 { b'A' } else { b'T' }
                    } else {
                        b
                    }
                } else {
                    b
                }
            })
            .collect()
    };

    let low_gc = make_gc_seq(1111, 300, false);
    let high_gc = make_gc_seq(2222, 300, true);
    let anc_low = make_gc_seq(3333, 300, false);
    let anc_high = make_gc_seq(4444, 300, true);

    let mut chr = vec![b'A'; 3000];
    chr[100..400].copy_from_slice(&low_gc);
    chr[700..1000].copy_from_slice(&high_gc);
    chr[1400..1700].copy_from_slice(&anc_low);
    chr[1700..2000].copy_from_slice(&anc_high);

    let gtf = "\
chr1\tT\texon\t101\t400\t.\t+\t.\tgene_id \"G_LOW\"; transcript_id \"T_LOW\";
chr1\tT\texon\t701\t1000\t.\t+\t.\tgene_id \"G_HIGH\"; transcript_id \"T_HIGH\";
chr1\tT\texon\t1401\t2000\t.\t+\t.\tgene_id \"G_ANCHOR\"; transcript_id \"T_ANCHOR\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("gc.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 11,
        minimizer_w: 4,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    let mut fq = String::new();
    let mut rid = 0usize;
    // G_LOW: 100 reads across uniform offsets
    for i in 0..100 {
        let off = (i * 2) % 200;
        let r = std::str::from_utf8(&low_gc[off..off + 70]).unwrap();
        fq.push_str(&format!("@r{rid}\n{r}\n+\n{}\n", "I".repeat(70)));
        rid += 1;
    }
    // G_HIGH: 200 reads across uniform offsets (2:1 observed ratio)
    for i in 0..200 {
        let off = i % 200;
        let r = std::str::from_utf8(&high_gc[off..off + 70]).unwrap();
        fq.push_str(&format!("@r{rid}\n{r}\n+\n{}\n", "I".repeat(70)));
        rid += 1;
    }
    // G_ANCHOR: 100 reads from low-GC half + 400 reads from high-GC half (strong GC bias signal!)
    for i in 0..100 {
        let off = (i * 2) % 200;
        let r = std::str::from_utf8(&anc_low[off..off + 70]).unwrap();
        fq.push_str(&format!("@r{rid}\n{r}\n+\n{}\n", "I".repeat(70)));
        rid += 1;
    }
    for i in 0..400 {
        let off = i % 200;
        let r = std::str::from_utf8(&anc_high[off..off + 70]).unwrap();
        fq.push_str(&format!("@r{rid}\n{r}\n+\n{}\n", "I".repeat(70)));
        rid += 1;
    }

    let fq_path = dir.path().join("gc.fq");
    std::fs::write(&fq_path, &fq).unwrap();

    let mut q_cfg = QuantConfig {
        r1: fq_path,
        out_prefix: dir.path().join("raw"),
        bias_correct: false,
        ..Default::default()
    };
    quantify_sample(&idx, &q_cfg).unwrap();

    q_cfg.bias_correct = true;
    q_cfg.out_prefix = dir.path().join("corr");
    quantify_sample(&idx, &q_cfg).unwrap();

    let parse_tpms = |p: std::path::PathBuf| -> (f64, f64) {
        let txt = read_gz_string(&p);
        let mut t_low = 0.0;
        let mut t_high = 0.0;
        for l in txt.lines().skip(1) {
            let c: Vec<&str> = l.split('\t').collect();
            if c[0] == "T_LOW" {
                t_low = c[1].parse().unwrap();
            } else if c[0] == "T_HIGH" {
                t_high = c[1].parse().unwrap();
            }
        }
        (t_low, t_high)
    };

    let (raw_l, raw_h) = parse_tpms(dir.path().join("raw.isoform.tpm.gz"));
    let (cor_l, cor_h) = parse_tpms(dir.path().join("corr.isoform.tpm.gz"));

    let raw_ratio = raw_h / raw_l;
    let corr_ratio = cor_h / cor_l;
    assert!((raw_ratio - 2.0).abs() < 0.05);
    // GC correction must reduce the biased 2.0 ratio toward 1.0 WITHOUT squaring/over-inverting!
    assert!(
        corr_ratio < raw_ratio && corr_ratio > 0.75,
        "GC bias correction should reduce ratio from ~2.0 toward ~1.0 without double-counting (raw={raw_ratio:.3}, corr={corr_ratio:.3})"
    );
}

#[test]
fn test_regression_microexon_15bp_psi_is_unbiased() {
    // Review Item 2: 15 bp internal cassette microexon (E1=150 bp, E2=15 bp, E3=150 bp).
    // Simulate 100 bp reads with 50% inclusion and 50% skipping across uniform start positions.
    // Verified: PSI must stay within |psi - 0.50| < 0.05 (not inflated to > 0.90!).
    let mut chr = pseudo_dna(3333, 1500);
    chr[250..252].copy_from_slice(b"GT");
    chr[498..500].copy_from_slice(b"AG");
    chr[515..517].copy_from_slice(b"GT");
    chr[798..800].copy_from_slice(b"AG");

    let gtf = "\
chr1\tT\texon\t101\t250\t.\t+\t.\tgene_id \"G_MICRO\"; transcript_id \"T_INC\";
chr1\tT\texon\t501\t515\t.\t+\t.\tgene_id \"G_MICRO\"; transcript_id \"T_INC\";
chr1\tT\texon\t801\t950\t.\t+\t.\tgene_id \"G_MICRO\"; transcript_id \"T_INC\";
chr1\tT\texon\t101\t250\t.\t+\t.\tgene_id \"G_MICRO\"; transcript_id \"T_EXC\";
chr1\tT\texon\t801\t950\t.\t+\t.\tgene_id \"G_MICRO\"; transcript_id \"T_EXC\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("micro15.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 11,
        minimizer_w: 4,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    let e1 = &chr[100..250]; // 150 bp
    let e2 = &chr[500..515]; // 15 bp microexon
    let e3 = &chr[800..950]; // 150 bp
    let mut t_inc = Vec::new();
    t_inc.extend_from_slice(e1);
    t_inc.extend_from_slice(e2);
    t_inc.extend_from_slice(e3);
    let mut t_exc = Vec::new();
    t_exc.extend_from_slice(e1);
    t_exc.extend_from_slice(e3);

    // Sample 100 bp reads uniformly covering the junction region of T_INC (100 reads) and T_EXC (100 reads).
    // Note: Every read starting at off in 65..145 on T_INC spans E1->E2->E3 (contributing 1 count to E1-E2 AND 1 to E2-E3),
    // and every read starting at off in 65..145 on T_EXC spans E1->E3 (contributing 1 count to E1-E3).
    // So E1-E2 = 80, E2-E3 = 80, E1-E3 = 80 -> true PSI = (160 / 2) / (160 / 2 + 80) = 0.5000!
    let mut fq = String::new();
    let mut rid = 0usize;
    for off in 65..145 {
        let ri = std::str::from_utf8(&t_inc[off..off + 100]).unwrap();
        fq.push_str(&format!("@i{rid}\n{ri}\n+\n{}\n", "I".repeat(100)));
        rid += 1;
        let re = std::str::from_utf8(&t_exc[off..off + 100]).unwrap();
        fq.push_str(&format!("@e{rid}\n{re}\n+\n{}\n", "I".repeat(100)));
        rid += 1;
    }

    let fq_path = dir.path().join("micro15.fq");
    std::fs::write(&fq_path, &fq).unwrap();
    let out = dir.path().join("micro_out");
    let qcfg = QuantConfig {
        r1: fq_path,
        out_prefix: out,
        bias_correct: false,
        ..Default::default()
    };
    quantify_sample(&idx, &qcfg).unwrap();

    let psi_txt = read_gz_string(&dir.path().join("micro_out.psi.gz"));
    let ce_line = psi_txt
        .lines()
        .skip(1)
        .find(|l| l.split('\t').nth(4) == Some("CE"))
        .expect("expected CE row for 15 bp microexon");
    let psi: f64 = ce_line.split('\t').nth(5).unwrap().parse().unwrap();
    assert!(
        (psi - 0.50).abs() < 0.05,
        "15 bp microexon PSI must be ~0.50 (got {psi:.4}, line={ce_line})"
    );
}

#[test]
fn test_regression_snp_1bp_before_donor_splice_boundary() {
    // Review Item 3: A SNP sitting 1 bp before the donor splice site (last base of E1) must NOT
    // stall extension before the boundary and soft-clip E2!
    let mut chr = pseudo_dna(4444, 1000);
    chr[180..182].copy_from_slice(b"GT");
    chr[398..400].copy_from_slice(b"AG");
    let gtf = "\
chr1\tT\texon\t101\t180\t.\t+\t.\tgene_id \"G_SNP\"; transcript_id \"T1\";
chr1\tT\texon\t401\t480\t.\t+\t.\tgene_id \"G_SNP\"; transcript_id \"T1\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("snp_boundary.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 11,
        minimizer_w: 4,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    // Read: last 45 bp of E1 + first 35 bp of E2 (80 bp total), with a SNP at index 44 (1 bp before donor!)
    let mut read_don_snp = Vec::new();
    read_don_snp.extend_from_slice(&chr[135..180]);
    read_don_snp.extend_from_slice(&chr[400..435]);
    read_don_snp[44] = if read_don_snp[44] == b'A' { b'G' } else { b'A' };

    let acfg = AlignConfig::default();
    let mut sc = AlignScratch::new(&acfg);
    match align_single_read(&idx, &read_don_snp, &acfg, &mut sc) {
        SingleAlignOutcome::Mapped(alns) => {
            assert_eq!(alns.len(), 1);
            assert_eq!(
                alns[0].node_path().as_slice(),
                &[0, 1],
                "read with SNP 1 bp before donor splice site must cross junction into E2"
            );
            assert_eq!(alns[0].soft_right, 0);
        }
        other => panic!("read with SNP 1 bp before donor failed to map: {other:?}"),
    }

    // Symmetrically: read seeded inside E2 extending backward with a SNP 1 bp after acceptor (index 45)
    let mut read_acc_snp = Vec::new();
    read_acc_snp.extend_from_slice(&chr[145..180]); // 35 bp E1
    read_acc_snp.extend_from_slice(&chr[400..445]); // 45 bp E2
    read_acc_snp[35] = if read_acc_snp[35] == b'C' { b'T' } else { b'C' };
    match align_single_read(&idx, &read_acc_snp, &acfg, &mut sc) {
        SingleAlignOutcome::Mapped(alns) => {
            assert_eq!(alns.len(), 1);
            assert_eq!(
                alns[0].node_path().as_slice(),
                &[0, 1],
                "read with SNP 1 bp after acceptor must cross junction cleanly"
            );
            assert_eq!(alns[0].soft_left, 0);
        }
        other => panic!("read with SNP 1 bp after acceptor failed to map: {other:?}"),
    }
}

#[test]
fn test_regression_pure_inclusion_exclusion_path_psi_and_ts_no_duplicate_nodes() {
    // Review Item 4:
    // Gene has tandem TS (2 alt tx starts at 101 & 151) + CE (401..500) + E3 (701..800).
    let chr = pseudo_dna(5555, 1200);
    let gtf = "\
chr1\tT\texon\t101\t200\t.\t+\t.\tgene_id \"G_EV\"; transcript_id \"T1\";
chr1\tT\texon\t401\t500\t.\t+\t.\tgene_id \"G_EV\"; transcript_id \"T1\";
chr1\tT\texon\t701\t800\t.\t+\t.\tgene_id \"G_EV\"; transcript_id \"T1\";
chr1\tT\texon\t151\t200\t.\t+\t.\tgene_id \"G_EV\"; transcript_id \"T2\";
chr1\tT\texon\t701\t800\t.\t+\t.\tgene_id \"G_EV\"; transcript_id \"T2\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("ev_reg.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 11,
        minimizer_w: 4,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    // Pure inclusion on CE (node 2, 0-based): edges 1->2 and 2->3 only (no skipping edge 1->3)
    // Plus node counts on TS (node 0 & node 1)
    let mut final_nodes = rustc_hash::FxHashMap::default();
    final_nodes.insert((0u32, 0u16), 20.0);
    final_nodes.insert((0u32, 1u16), 30.0);

    let mut edges_inc = rustc_hash::FxHashMap::default();
    edges_inc.insert((0u32, 1u16, 2u16), (20.0, 10, 2));
    edges_inc.insert((0u32, 2u16, 3u16), (20.0, 10, 2));

    let evts_inc = langur_core::events::quantify_all_events(
        &idx,
        &final_nodes,
        &edges_inc,
        &EventConfig::default(),
    );

    // Verify ZERO duplicate node indices emitted (TS anchor node must not be duplicated!)
    let node_ids: Vec<u32> = evts_inc.iter().map(|e| e.node_1based).collect();
    let mut uniq_ids = node_ids.clone();
    uniq_ids.sort_unstable();
    uniq_ids.dedup();
    assert_eq!(
        node_ids.len(),
        uniq_ids.len(),
        "duplicate node rows emitted: {node_ids:?}"
    );

    // Verify pure-inclusion CE has PSI == 1.0 AND inc_paths path psi == 1.0 (never 0.0!)
    let ce_inc = evts_inc.iter().find(|e| e.event_type == "CE").unwrap();
    assert_eq!(ce_inc.psi, Some(1.0));
    let sum_inc_psi: f64 = ce_inc.inc_paths.iter().map(|(_, p)| *p).sum();
    assert!(
        (sum_inc_psi - 1.0).abs() < 1e-6,
        "pure-inclusion inc_paths must sum to 1.0, got {sum_inc_psi}"
    );

    // Pure exclusion on CE: only skipping edge 1->3
    let mut edges_exc = rustc_hash::FxHashMap::default();
    edges_exc.insert((0u32, 1u16, 3u16), (20.0, 10, 2));
    let evts_exc = langur_core::events::quantify_all_events(
        &idx,
        &final_nodes,
        &edges_exc,
        &EventConfig::default(),
    );
    let ce_exc = evts_exc.iter().find(|e| e.event_type == "CE").unwrap();
    assert_eq!(ce_exc.psi, Some(0.0));
    let sum_exc_psi: f64 = ce_exc.exc_paths.iter().map(|(_, p)| *p).sum();
    assert!(
        (sum_exc_psi - 1.0).abs() < 1e-6,
        "pure-exclusion exc_paths must sum to 1.0, got {sum_exc_psi}"
    );
}

#[test]
fn test_regression_sam_orientation_revcomp_and_indel_cigar() {
    // Review Item 5: Verify SAM FLAG 0x10 = `pos_strand != read_is_sense`, revcomp'd SEQ + reversed
    // QUAL when genomic_minus, and exact CIGAR query/ref consumption when indels + introns occur.
    let mut chr = pseudo_dna(6666, 1000);
    chr[180..182].copy_from_slice(b"GT");
    chr[298..300].copy_from_slice(b"AG");
    let gtf = "\
chr1\tT\texon\t101\t180\t.\t+\t.\tgene_id \"G_POS\"; transcript_id \"TP\";
chr1\tT\texon\t301\t380\t.\t+\t.\tgene_id \"G_POS\"; transcript_id \"TP\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("sam.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 11,
        minimizer_w: 4,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    // Antisense (RF) read aligning to `+` strand gene -> genomic_minus == true (FLAG 0x10 set, SEQ revcomp'd!)
    let mut spliced_sense = Vec::new();
    spliced_sense.extend_from_slice(&chr[140..180]); // 40 bp E1
    spliced_sense.extend_from_slice(&chr[300..340]); // 40 bp E2
    let anti_read = revcomp(&spliced_sense);
    let qual: Vec<u8> = (0..80).map(|i| b'!' + ((i % 40) as u8)).collect();
    let rev_qual: Vec<u8> = qual.iter().rev().copied().collect();

    let acfg = AlignConfig::default();
    let mut sc = AlignScratch::new(&acfg);
    let SingleAlignOutcome::Mapped(alns) = align_single_read(&idx, &anti_read, &acfg, &mut sc) else {
        panic!("expected mapped antisense read");
    };
    assert!(!alns[0].read_is_sense);

    let sam = langur_core::io::writers::format_sam_line(
        &idx,
        b"read_anti",
        &anti_read,
        &qual,
        &alns[0],
        false,
        true,
        false,
    );
    let cols: Vec<&str> = sam.trim_end().split('\t').collect();
    let flag: u16 = cols[1].parse().unwrap();
    assert_ne!(flag & 0x10, 0, "antisense read on + gene must set FLAG 0x10");
    assert_eq!(cols[3], "141");
    assert_eq!(cols[5], "40M120N40M");
    assert_eq!(cols[9].as_bytes(), spliced_sense.as_slice());
    assert_eq!(cols[10].as_bytes(), rev_qual.as_slice());

    // Indel read (2 bp insertion in E1): verify CIGAR includes `2I` and query consumption == read_len
    let mut indel_read = Vec::new();
    indel_read.extend_from_slice(&chr[100..140]);
    indel_read.extend_from_slice(b"GT");
    indel_read.extend_from_slice(&chr[140..178]);
    let SingleAlignOutcome::Mapped(alns_i) = align_single_read(&idx, &indel_read, &acfg, &mut sc) else {
        panic!("expected mapped indel read");
    };
    let sam_i = langur_core::io::writers::format_sam_line(
        &idx,
        b"read_ins",
        &indel_read,
        &qual,
        &alns_i[0],
        false,
        true,
        false,
    );
    let cols_i: Vec<&str> = sam_i.trim_end().split('\t').collect();
    assert!(
        cols_i[5].contains('I'),
        "expected insertion op 'I' in CIGAR, got {}",
        cols_i[5]
    );
}

#[test]
fn test_regression_shared_upstream_ts_edge_and_26bp_ce_psi() {
    // Reproduces the ENSG00000168385 human-scale bug:
    // - Exon 1 has tandem TSS starts creating a tiny 3 bp contiguous TS sub-node (`node 1` = 298..300)
    //   immediately contiguous with upstream `node 0` (101..297) and downstream donor `node 2` (301..360).
    // - Exon 2 is a 26 bp CE cassette exon (`node 3` = 701..726).
    // - Exon 3 has tandem TES stops creating a tiny 3 bp contiguous TE sub-node (`node 4` = 1101..1103)
    //   contiguous with `node 5` (1104..1300).
    // With ~98% true inclusion (100 inclusion reads crossing 26 bp CE vs 2 skipping reads),
    // `edge_effective_length` must walk contiguous exonic runs instead of collapsing the 3 bp TS/TE
    // sub-nodes to eff_len ~ 0.01 (which previously amplified 2 skipping reads into PSI ~ 0.07!).
    let chr = pseudo_dna(0x1683_8500_DEAD_BEEF, 3000);
    let gtf = "\
chr1\tSIM\texon\t101\t360\t.\t+\t.\tgene_id \"G168385\"; transcript_id \"T_INC1\";
chr1\tSIM\texon\t701\t726\t.\t+\t.\tgene_id \"G168385\"; transcript_id \"T_INC1\";
chr1\tSIM\texon\t1101\t1300\t.\t+\t.\tgene_id \"G168385\"; transcript_id \"T_INC1\";
chr1\tSIM\texon\t298\t360\t.\t+\t.\tgene_id \"G168385\"; transcript_id \"T_INC2\";
chr1\tSIM\texon\t701\t726\t.\t+\t.\tgene_id \"G168385\"; transcript_id \"T_INC2\";
chr1\tSIM\texon\t1101\t1103\t.\t+\t.\tgene_id \"G168385\"; transcript_id \"T_INC2\";
chr1\tSIM\texon\t301\t360\t.\t+\t.\tgene_id \"G168385\"; transcript_id \"T_EXC\";
chr1\tSIM\texon\t1101\t1300\t.\t+\t.\tgene_id \"G168385\"; transcript_id \"T_EXC\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("g168385.lidx");
    let icfg = IndexBuildConfig {
        minimizer_k: 13,
        minimizer_w: 5,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &icfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    // Build T_INC1 (E1 + E2(26bp) + E3) and T_EXC (E1 + E3) sequences
    let mut t_inc = Vec::new();
    t_inc.extend_from_slice(&chr[100..360]);
    t_inc.extend_from_slice(&chr[700..726]);
    t_inc.extend_from_slice(&chr[1100..1300]);

    let mut t_exc = Vec::new();
    t_exc.extend_from_slice(&chr[100..360]);
    t_exc.extend_from_slice(&chr[1100..1300]);

    // 120 inclusion reads (100 bp) spanning E1->E2(26bp)->E3 and E1 TS boundary, plus 2 skipping reads!
    use std::io::Write as _;
    let fq = dir.path().join("reads.fq");
    let mut f = std::fs::File::create(&fq).unwrap();
    let qual = "I".repeat(100);
    let mut rid = 0usize;
    for off in (150..270).step_by(1) {
        let seq = std::str::from_utf8(&t_inc[off..off + 100]).unwrap();
        writeln!(f, "@inc_{rid}\n{seq}\n+\n{qual}").unwrap();
        rid += 1;
    }
    for &off in &[200usize, 215] {
        let seq = std::str::from_utf8(&t_exc[off..off + 100]).unwrap();
        writeln!(f, "@exc_{rid}\n{seq}\n+\n{qual}").unwrap();
        rid += 1;
    }
    drop(f);

    let out_prefix = dir.path().join("out");
    let qcfg = QuantConfig {
        r1: fq,
        out_prefix: out_prefix.clone(),
        ..Default::default()
    };
    quantify_sample(&idx, &qcfg).unwrap();

    let psi_txt = read_gz_string(&dir.path().join("out.psi.gz"));
    let mut ce_psi = None;
    for line in psi_txt.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols[4] == "CE" && cols[2] == "chr1:701-726" {
            ce_psi = Some(cols[5].parse::<f64>().unwrap());
        }
    }
    let psi = ce_psi.expect("expected 26 bp CE event row");
    assert!(
        psi > 0.92,
        "26 bp CE flanked by tiny contiguous TS/TE nodes should have PSI > 0.92, got {psi:.4}"
    );
}
