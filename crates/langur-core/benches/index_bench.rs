use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::io::{BufReader, Cursor};
use tempfile::tempdir;
use langur_core::align::{align_single_read, AlignConfig, AlignScratch, SingleAlignOutcome};
use langur_core::annotation::{load_gtf_reader, GtfOptions};
use langur_core::delta::{build_two_group_design, fit_glm_cohort, GlmEventData, GlmOptions};
use langur_core::events::gibbs::{
    sample_complex_event_gibbs, DEFAULT_GIBBS_BURN_IN, DEFAULT_GIBBS_THIN,
};
use langur_core::events::ReducedPath;
use langur_core::index::{build_index_from_memory, DecoyMode, IndexBuildConfig, MmapIndex};
use langur_core::seq::extract_minimizers;

fn synthetic_dataset(num_genes: usize) -> (Vec<u8>, String) {
    // Deterministic pseudo-random DNA chromosome with `num_genes` multi-exon genes
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let bases = *b"ACGT";
    let chrom_len = num_genes * 1200 + 2000;
    let mut chrom = Vec::with_capacity(chrom_len);
    for _ in 0..chrom_len {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        let b = bases[((state.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 62) & 3) as usize];
        chrom.push(b);
    }

    let mut gtf = String::with_capacity(num_genes * 300);
    for g in 0..num_genes {
        let base = g * 1200 + 100;
        let strand = if g % 2 == 0 { '+' } else { '-' };
        // Isoform 1: 3 exons (150 bp each)
        gtf.push_str(&format!(
            "chr1\tSIM\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"G{g}\"; transcript_id \"T{g}_1\";\n",
            base,
            base + 149
        ));
        gtf.push_str(&format!(
            "chr1\tSIM\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"G{g}\"; transcript_id \"T{g}_1\";\n",
            base + 300,
            base + 449
        ));
        gtf.push_str(&format!(
            "chr1\tSIM\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"G{g}\"; transcript_id \"T{g}_1\";\n",
            base + 600,
            base + 749
        ));
        // Isoform 2: skips middle exon
        gtf.push_str(&format!(
            "chr1\tSIM\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"G{g}\"; transcript_id \"T{g}_2\";\n",
            base,
            base + 149
        ));
        gtf.push_str(&format!(
            "chr1\tSIM\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"G{g}\"; transcript_id \"T{g}_2\";\n",
            base + 600,
            base + 749
        ));
    }
    (chrom, gtf)
}

fn bench_phase1_index(c: &mut Criterion) {
    let (chrom, gtf) = synthetic_dataset(500);
    let refset = load_gtf_reader(BufReader::new(Cursor::new(&gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("bench.lidx");
    let cfg = IndexBuildConfig {
        minimizer_k: 17,
        minimizer_w: 8,
        decoy_mode: DecoyMode::SelectiveGenome,
        ..Default::default()
    };

    c.bench_function("index_build_500_genes_with_decoys", |b| {
        b.iter(|| {
            build_index_from_memory(&[("chr1", &chrom)], &refset, &idx_path, &cfg).unwrap();
        });
    });

    build_index_from_memory(&[("chr1", &chrom)], &refset, &idx_path, &cfg).unwrap();

    c.bench_function("index_mmap_open_crc32", |b| {
        b.iter(|| {
            let idx = MmapIndex::open(black_box(&idx_path)).unwrap();
            black_box(idx.header().num_genes);
        });
    });

    let idx = MmapIndex::open(&idx_path).unwrap();
    let mut query_hashes = Vec::new();
    let mut mins = Vec::new();
    for g in 0..50 {
        extract_minimizers(idx.gene_seq(g), 17, 8, idx.hash_bits(), &mut mins);
        query_hashes.extend(mins.iter().map(|m| m.hash));
    }

    c.bench_function("minimizer_lookup_batch", |b| {
        b.iter(|| {
            let mut hits = 0usize;
            for &h in &query_hashes {
                if idx.lookup_minimizer(black_box(h)).is_some() {
                    hits += 1;
                }
            }
            black_box(hits)
        });
    });

    // Phase 2 short-read spliced alignment throughput benchmark (1,000 100-bp spliced & exonic reads)
    let mut reads: Vec<Vec<u8>> = Vec::with_capacity(1000);
    for g in 0..500u32 {
        let gseq = idx.gene_seq(g);
        if gseq.len() >= 300 {
            // Spliced read crossing exon1 -> exon2 (offset 100..200)
            reads.push(gseq[100..200].to_vec());
            // Spliced read skipping exon2: last 50 bp of exon1 + first 50 bp of exon3
            let mut skip = Vec::with_capacity(100);
            skip.extend_from_slice(&gseq[100..150]);
            skip.extend_from_slice(&gseq[300..350]);
            reads.push(skip);
        }
    }
    let acfg = AlignConfig::default();
    let mut scratch = AlignScratch::new(&acfg);

    c.bench_function("quant_align_1000_spliced_100bp_reads_1thread", |b| {
        b.iter(|| {
            let mut mapped = 0usize;
            for r in &reads {
                if matches!(
                    align_single_read(&idx, black_box(r), &acfg, &mut scratch),
                    SingleAlignOutcome::Mapped(_)
                ) {
                    mapped += 1;
                }
            }
            black_box(mapped)
        });
    });

    // Phase 3 Beta-Binomial 3-pass EB GLM benchmark (500 events, 12 samples, parallel rayon)
    let files_a: Vec<_> = (0..6).map(|i| format!("a{i}.psi.gz").into()).collect();
    let files_b: Vec<_> = (0..6).map(|i| format!("b{i}.psi.gz").into()).collect();
    let design = build_two_group_design(files_a, files_b, "A", "B").unwrap();
    let glm_events: Vec<GlmEventData> = (0..500)
        .map(|i| {
            let psi = if i % 5 == 0 {
                vec![0.75, 0.72, 0.78, 0.74, 0.76, 0.73, 0.25, 0.28, 0.22, 0.26, 0.24, 0.27]
            } else {
                vec![0.48, 0.52, 0.49, 0.51, 0.47, 0.53, 0.50, 0.49, 0.51, 0.48, 0.52, 0.50]
            };
            GlmEventData {
                gene: format!("G{i}"),
                node: "2".into(),
                coord: "chr1:100-200".into(),
                strand: "+".into(),
                etype: "CE".into(),
                complexity: 1,
                entropy: 1.0,
                valid_idx: (0..12).collect(),
                psi,
                eff_reads: vec![80.0; 12],
                samp_a: 6,
                samp_b: 6,
                mean_reads_a: 80.0,
                mean_reads_b: 80.0,
            }
        })
        .collect();
    let gopts = GlmOptions::default();

    c.bench_function("delta_eb_glm_500_events_12samples", |b| {
        b.iter(|| {
            let (res, fdr, _) = fit_glm_cohort(black_box(&glm_events), &design, &gopts);
            black_box((res.len(), fdr.len()))
        });
    });

    // Phase 4 targeted collapsed Gibbs sampler benchmark (complex K_3 event: 4 inc + 4 exc paths, 3 shared edge groups, 100 draws)
    let k3_inc: Vec<ReducedPath> = (0..4)
        .map(|i| ReducedPath {
            nodes: vec![0, 1 + i, 9],
            count: 15.0 + (i as f64) * 5.0,
            eff_len: 120.0,
            psi: 0.125,
        })
        .collect();
    let k3_exc: Vec<ReducedPath> = (0..4)
        .map(|i| ReducedPath {
            nodes: vec![0, 5 + i, 9],
            count: 10.0 + (i as f64) * 4.0,
            eff_len: 100.0,
            psi: 0.125,
        })
        .collect();
    let k3_ambig = vec![
        (vec![0, 1, 4, 5], 180.5),
        (vec![2, 3, 6], 45.0),
        (vec![0, 2, 4, 7], 92.25),
    ];
    let k3_init = vec![0.125; 8];

    c.bench_function("gibbs_complex_k3_event_100draws", |b| {
        b.iter(|| {
            let s = sample_complex_event_gibbs(
                black_box(&k3_inc),
                black_box(&k3_exc),
                black_box(&k3_ambig),
                black_box(&k3_init),
                None,
                100,
                DEFAULT_GIBBS_BURN_IN,
                DEFAULT_GIBBS_THIN,
                42,
                17,
                4,
            );
            black_box(s)
        });
    });
}

criterion_group!(benches, bench_phase1_index);
criterion_main!(benches);
