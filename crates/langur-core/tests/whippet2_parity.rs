//! Phase 1 integration, property, and Whippet2.jl topology parity test suite.
//! Reproduces the canonical splice-graph fixtures from `Whippet2.jl/test/runtests.jl`
//! plus targeted tests for short-node contiguous splice flanks, bijective 54-bit minimizer
//! hashing, `±(k-1)` selective decoy masking, and zero-copy `.lidx` mmap integrity.

use flate2::read::MultiGzDecoder;
use proptest::prelude::*;
use std::io::{BufReader, Cursor, Read};
use tempfile::tempdir;
use langur_core::annotation::{load_gtf_reader, GtfOptions};
use langur_core::graph::{build_splice_graph, EdgeMotif, EdgeType};
use langur_core::index::{
    build_index_from_files, build_index_from_memory, write_exons_tsv_gz, write_paths_tsv_gz,
    DecoyMode, IndexBuildConfig, MinimizerLookup, MmapIndex, SeedHitKind, FLAG_HASH_MAY_COLLIDE,
};
use langur_core::seq::{
    encode_flank_lsb_from_end, encode_flank_msb, extract_minimizers, hash_kmer, is_bijective_hash,
    mix_bits, revcomp,
};

/// Exact GTF fixture from `Whippet2.jl/test/runtests.jl` lines 142-161.
const WHIPPET2_TEST_GTF: &str = "\
# gtf file test
chr0\tTEST\texon\t6\t20\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t31\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t54\t62\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t76\t85\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t6\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t51\t62\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t76\t85\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t6\t20\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t31\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t54\t65\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t76\t90\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t11\t20\t.\t+\t.\tgene_id \"single\"; transcript_id \"ex1_single\";
chr0\tTEST\texon\t11\t20\t.\t-\t.\tgene_id \"single_rev\"; transcript_id \"single_rev\";
chr0\tTEST\texon\t11\t20\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"def_kiss\";
chr0\tTEST\texon\t21\t30\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"def_kiss\";
chr0\tTEST\texon\t11\t30\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"ret_kiss\";
chr0\tTEST\texon\t11\t30\t.\t-\t.\tgene_id \"spl_neg\"; transcript_id \"spl_negative\";
chr0\tTEST\texon\t92\t97\t.\t-\t.\tgene_id \"spl_neg\"; transcript_id \"spl_negative\";
";

/// Exact synthetic genome (`fwd * rev`, 200 bp) from `Whippet2.jl/test/runtests.jl` lines 205-225.
fn whippet2_test_genome() -> Vec<u8> {
    let mut fwd = Vec::new();
    fwd.extend_from_slice(b"AAAAA"); // 1-5   buffer1
    fwd.extend_from_slice(b"TTATT"); // 6-10  utr5
    fwd.extend_from_slice(b"GCGGATTACA"); // 11-20 exon1
    fwd.extend_from_slice(b"TTTTTTTTTT"); // 21-30 int1
    fwd.extend_from_slice(b"GCATTAGAAG"); // 31-40 exon2
    fwd.extend_from_slice(b"GGGGGGGGGG"); // 41-50 int2
    fwd.extend_from_slice(b"CCT"); // 51-53 exon3alt3
    fwd.extend_from_slice(b"CTATGCTAG"); // 54-62 exon3def
    fwd.extend_from_slice(b"TTC"); // 63-65 exon3alt5
    fwd.extend_from_slice(b"CCCCCCCCCC"); // 66-75 int3
    fwd.extend_from_slice(b"TTAGACAAGA"); // 76-85 exon4
    fwd.extend_from_slice(b"AATAA"); // 86-90 apa
    fwd.extend_from_slice(b"AAAAAAAAAA"); // 91-100 buffer2
    assert_eq!(fwd.len(), 100);
    let rev = revcomp(&fwd);
    let mut genome = fwd;
    genome.extend_from_slice(&rev);
    genome
}

#[test]
fn test_whippet2_csg_topology_and_sequence_parity() {
    let refset = load_gtf_reader(BufReader::new(Cursor::new(WHIPPET2_TEST_GTF)), &GtfOptions::default())
        .expect("GTF parse failed");

    // 1. Gene annotations parity (Whippet2.jl/test/runtests.jl L184-187)
    let g_one = refset.get("one").expect("missing gene 'one'");
    assert_eq!(g_one.don, vec![20, 40, 62, 65]);
    assert_eq!(g_one.acc, vec![31, 51, 54, 76]);
    assert!(g_one.novel_don.is_empty());
    assert!(g_one.novel_acc.is_empty());

    // 2. Graph building parity (Whippet2.jl/test/runtests.jl L227-255)
    let genome = whippet2_test_genome();
    let sg_one = build_splice_graph(g_one, &genome).unwrap();
    let sg_sin = build_splice_graph(refset.get("single").unwrap(), &genome).unwrap();
    let sg_rev = build_splice_graph(refset.get("single_rev").unwrap(), &genome).unwrap();
    let sg_kis = build_splice_graph(refset.get("kissing").unwrap(), &genome).unwrap();
    let sg_neg = build_splice_graph(refset.get("spl_neg").unwrap(), &genome).unwrap();

    let expected_one = b"TTATTGCGGATTACATTTTTTTTTTGCATTAGAAGCCTCTATGCTAGTTCTTAGACAAGAAATAA";
    let expected_sin = b"GCGGATTACA";
    let expected_kis = b"AAAAAAAAAATGTAATCCGC";

    assert_eq!(sg_one.seq, expected_one);
    assert_eq!(sg_sin.seq, expected_sin);
    assert_eq!(sg_rev.seq, b"TGTAATCCGC");
    assert_eq!(sg_kis.seq, expected_kis);
    assert_eq!(sg_neg.nodes.len(), 2);

    // Annotated paths (convert 0-based local node IDs to 1-based to match Whippet2 BitSets)
    assert_eq!(sg_one.transcripts.len(), 3);
    let paths_1based: Vec<Vec<u16>> = sg_one
        .transcripts
        .iter()
        .map(|t| t.nodes.iter().map(|&n| n + 1).collect())
        .collect();
    assert_eq!(paths_1based[0], vec![1, 3, 5, 7]); // def
    assert_eq!(paths_1based[1], vec![1, 2, 3, 4, 5, 7]); // int1_alt3
    assert_eq!(paths_1based[2], vec![1, 3, 5, 6, 7, 8]); // apa_alt5

    // Verify exact CSG node coordinates, lengths, EdgeTypes, and Event Motifs on gene 'one'
    assert_eq!(sg_one.nodes.len(), 8);
    // Node 1 (0-based 0): 6..20 (len 15), SL..LL
    assert_eq!((sg_one.nodes[0].genomic_coord, sg_one.nodes[0].len), (6, 15));
    // Node 2 (0-based 1): int1 21..30 (len 10), LL..RR -> RI!
    assert_eq!((sg_one.nodes[1].genomic_coord, sg_one.nodes[1].len), (21, 10));
    assert_eq!(
        EdgeMotif::from_edges(sg_one.nodes[1].left_edge, sg_one.nodes[1].right_edge),
        EdgeMotif::RI
    );
    // Node 4 (0-based 3): exon3alt3 51..53 (len 3), LR..RR -> AA!
    assert_eq!((sg_one.nodes[3].genomic_coord, sg_one.nodes[3].len), (51, 3));
    assert_eq!(
        EdgeMotif::from_edges(sg_one.nodes[3].left_edge, sg_one.nodes[3].right_edge),
        EdgeMotif::AA
    );
    // Node 5 (0-based 4): exon3def 54..62 (len 9), RR..LL -> CE!
    assert_eq!((sg_one.nodes[4].genomic_coord, sg_one.nodes[4].len), (54, 9));
    assert_eq!(
        EdgeMotif::from_edges(sg_one.nodes[4].left_edge, sg_one.nodes[4].right_edge),
        EdgeMotif::CE
    );
    // Node 6 (0-based 5): exon3alt5 63..65 (len 3), LL..LR -> AD!
    assert_eq!((sg_one.nodes[5].genomic_coord, sg_one.nodes[5].len), (63, 3));
    assert_eq!(
        EdgeMotif::from_edges(sg_one.nodes[5].left_edge, sg_one.nodes[5].right_edge),
        EdgeMotif::AD
    );
    // Node 8 (0-based 7): apa 86..90 (len 5), RS..RS -> TE!
    assert_eq!((sg_one.nodes[7].genomic_coord, sg_one.nodes[7].len), (86, 5));
    assert_eq!(
        EdgeMotif::from_edges(sg_one.nodes[7].left_edge, sg_one.nodes[7].right_edge),
        EdgeMotif::TE
    );

    // Refinement 1 check: short 3 bp AA/AD nodes extend flanks across contiguous exonic sequence!
    // Node 3 (exon3alt3, 3 bp "CCT") + Node 4 (exon3def, 9 bp "CTATGCTAG") + Node 5 (exon3alt5, 3 bp "TTC")
    // are genomically contiguous (51..53, 54..62, 63..65) -> full 15 bp flank "CCTCTATGCTAGTTC"!
    assert_eq!(sg_one.nodes[3].acceptor_flank_len, 15);
    let (exp_acc, exp_acc_len) = encode_flank_msb(b"CCTCTATGCTAGTTC");
    assert_eq!(sg_one.nodes[3].acceptor_flank_2bit, exp_acc);
    assert_eq!(sg_one.nodes[3].acceptor_flank_len, exp_acc_len);

    assert_eq!(sg_one.nodes[5].donor_flank_len, 15);
    let (exp_don, exp_don_len) = encode_flank_lsb_from_end(b"CCTCTATGCTAGTTC");
    assert_eq!(sg_one.nodes[5].donor_flank_2bit, exp_don);
    assert_eq!(sg_one.nodes[5].donor_flank_len, exp_don_len);
}

#[test]
fn test_standalone_microexon_flank_never_bleeds_across_intron_or_zero_pads() {
    // 6 bp microexon flanked by hard introns on both sides (non-contiguous with neighbors)
    let gtf = "\
chr1\tT\texon\t10\t29\t.\t+\t.\tgene_id \"micro\"; transcript_id \"tx1\";
chr1\tT\texon\t50\t55\t.\t+\t.\tgene_id \"micro\"; transcript_id \"tx1\";
chr1\tT\texon\t80\t99\t.\t+\t.\tgene_id \"micro\"; transcript_id \"tx1\";
";
    let mut chr = vec![b'A'; 120];
    chr[9..29].copy_from_slice(b"CGTACGTACGTACGTACGTA");
    chr[49..55].copy_from_slice(b"GGCCCT"); // 6 bp microexon
    chr[79..99].copy_from_slice(b"TTGCATTGCATTGCATTGCA");

    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let sg = build_splice_graph(refset.get("micro").unwrap(), &chr).unwrap();
    assert_eq!(sg.nodes.len(), 3);
    let micro = &sg.nodes[1];
    assert_eq!(micro.len, 6);
    assert_eq!(micro.left_edge, EdgeType::LR);
    assert_eq!(micro.right_edge, EdgeType::LR);
    // Must NOT bleed into exon 1 or exon 3 across introns!
    assert_eq!(micro.acceptor_flank_len, 6);
    assert_eq!(micro.donor_flank_len, 6);

    // Serialize into .lidx and test NodeRecord::acceptor_flank_mismatches / donor_flank_mismatches
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("micro.lidx");
    let cfg = IndexBuildConfig {
        minimizer_k: 5,
        minimizer_w: 3,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &cfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();
    let nd = &idx.gene_nodes(0)[1];
    assert_eq!(nd.acceptor_flank_len, 6);
    assert_eq!(nd.donor_flank_len, 6);

    // Query with 12 bp where bases 7..12 are TTTTTT (not AAAAAA): must clamp to 6 bp -> 0 mismatches!
    let (q_exact, _) = encode_flank_msb(b"GGCCCTTTTTTT");
    assert_eq!(nd.acceptor_flank_mismatches(q_exact, 12), Some((6, 0)));
    // 1 SNP inside the 6 bp microexon -> exactly 1 mismatch
    let (q_snp, _) = encode_flank_msb(b"GGCACTTTTTTT");
    assert_eq!(nd.acceptor_flank_mismatches(q_snp, 12), Some((6, 1)));
    // Donor side
    let (q_don, _) = encode_flank_lsb_from_end(b"TTTTTTGGCCCT");
    assert_eq!(nd.donor_flank_mismatches(q_don, 12), Some((6, 0)));
    let (q_don_snp, _) = encode_flank_lsb_from_end(b"TTTTTTGGCACT");
    assert_eq!(nd.donor_flank_mismatches(q_don_snp, 12), Some((6, 1)));
}

#[test]
fn test_selective_decoy_boundary_masking_and_ri_protection() {
    // Construct a 600 bp chromosome:
    // - Gene G1 (+):
    //     Exon 1:  50..89   (40 bp, distinct sequence E1)
    //     Intron1: 90..129  (40 bp, annotated RI node via tx_ri -> MUST NOT become a decoy!)
    //     Exon 2: 130..169  (40 bp, distinct sequence E2)
    //     Intron2: 170..369 (200 bp, unannotated intron):
    //         - At 170..178 (within k-1=8 bp of Exon 2 boundary): overlaps masked pad -> MUST NOT be decoy
    //         - At 240..279: exact paralogous copy of E1 -> MUST be extracted as INTRONIC decoy!
    //     Exon 3: 370..409  (40 bp, distinct sequence E3)
    // - Intergenic region 450..580:
    //     - At 500..539: exact pseudogene copy of E2 -> MUST be extracted as INTERGENIC decoy!
    let k = 9usize;
    let w = 4usize;

    let e1 = b"GCTAGCTTAGCCATGCTAGCATCGATCGATGCTAGCTAGG"; // 40 bp
    let ri = b"CCCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCC"; // 40 bp (disjoint from e1/e2/e3)
    let e2 = b"TGCATGCATTGCATGCATGGCCGGAATTCCAAGGTTAACC"; // 40 bp
    let e3 = b"AACCGGAATTCCGGTTAACCGGAATTCCGGTTCCAATTGG"; // 40 bp

    // Background sequence made of non-colliding poly-A/T-like low-complexity distinct chunks
    let mut chr = vec![b'A'; 600];
    chr[49..89].copy_from_slice(e1);
    chr[89..129].copy_from_slice(ri);
    chr[129..169].copy_from_slice(e2);
    // Put a copy of E1 deep inside unannotated Intron 2 (1-based 240..279)
    chr[239..279].copy_from_slice(e1);
    chr[369..409].copy_from_slice(e3);
    // Put a copy of E2 in intergenic region (1-based 500..539)
    chr[499..539].copy_from_slice(e2);

    let gtf = "\
chr1\tT\texon\t50\t89\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx_main\";
chr1\tT\texon\t130\t169\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx_main\";
chr1\tT\texon\t370\t409\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx_main\";
chr1\tT\texon\t50\t169\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx_ri\";
chr1\tT\texon\t370\t409\t.\t+\t.\tgene_id \"G1\"; transcript_id \"tx_ri\";
";

    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("decoy_test.lidx");
    let cfg = IndexBuildConfig {
        minimizer_k: k,
        minimizer_w: w,
        decoy_mode: DecoyMode::SelectiveGenome,
        decoy_window_pad: 30,
        ..Default::default()
    };
    let stats = build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &cfg).unwrap();
    assert_eq!(stats.num_hash_collisions, 0);
    assert!(stats.num_decoys >= 2, "expected both intronic and intergenic decoys");

    let idx = MmapIndex::open(&idx_path).unwrap();
    assert_eq!(idx.header().flags & FLAG_HASH_MAY_COLLIDE, 0);

    // 1. Verify NO decoy overlaps ANY CSG node or its ±(k-1) boundary pad!
    let pad = (k - 1) as u32;
    let csg_nodes = idx.gene_nodes(0);
    for d in idx.decoys() {
        let ds = d.genomic_start;
        let de = d.genomic_start + d.seq_len - 1;
        for nd in csg_nodes {
            let ms = nd.genomic_coord.saturating_sub(pad);
            let me = nd.genomic_coord + nd.len - 1 + pad;
            assert!(
                de < ms || ds > me,
                "decoy [{ds}..{de}] violated ±(k-1) mask [{ms}..{me}] of CSG node!"
            );
        }
        // Verify 2-bit unpacked decoy sequence matches the chromosome exactly
        let mut unpacked = Vec::new();
        idx.unpack_decoy_slice(
            // find index
            idx.decoys()
                .iter()
                .position(|x| x.genomic_start == d.genomic_start)
                .unwrap() as u32,
            0,
            d.seq_len as usize,
            &mut unpacked,
        );
        assert_eq!(unpacked, &chr[(ds - 1) as usize..de as usize]);
    }

    // 2. Verify E1 minimizers return BOTH CSG hit and decoy hit!
    let mut mins = Vec::new();
    extract_minimizers(e1, k, w, idx.hash_bits(), &mut mins);
    assert!(!mins.is_empty());
    let mut found_decoy_on_e1 = false;
    for m in &mins {
        if let Some(MinimizerLookup::Postings { csg, decoys, .. }) = idx.lookup_minimizer(m.hash) {
            if !csg.is_empty() && !decoys.is_empty() {
                assert_eq!(decoys[0].kind(), SeedHitKind::Decoy);
                found_decoy_on_e1 = true;
            }
        }
    }
    assert!(found_decoy_on_e1, "expected E1 minimizer to return both CSG and decoy hits");

    // 3. Verify RI node minimizers have ZERO decoy hits (annotated RI is protected!)
    extract_minimizers(ri, k, w, idx.hash_bits(), &mut mins);
    for m in &mins {
        match idx.lookup_minimizer(m.hash) {
            Some(MinimizerLookup::Inline(hit)) => assert_eq!(hit.kind(), SeedHitKind::CsgNode),
            Some(MinimizerLookup::Postings { decoys, .. }) => {
                assert!(decoys.is_empty(), "RI minimizer must not have decoy hits");
            }
            None => {}
        }
    }
}

#[test]
fn test_lidx_roundtrip_crc_corruption_and_aux_tsv_outputs() {
    let dir = tempdir().unwrap();
    let fa_path = dir.path().join("genome.fa");
    let gtf_path = dir.path().join("anno.gtf");
    let idx_path = dir.path().join("graph.lidx");
    let exons_path = dir.path().join("graph.exons.tsv.gz");
    let paths_path = dir.path().join("graph.paths.tsv.gz");

    let genome = whippet2_test_genome();
    std::fs::write(&fa_path, format!(">chr0\n{}\n", std::str::from_utf8(&genome).unwrap())).unwrap();
    std::fs::write(&gtf_path, WHIPPET2_TEST_GTF).unwrap();

    let cfg = IndexBuildConfig {
        minimizer_k: 5,
        minimizer_w: 3,
        decoy_mode: DecoyMode::SelectiveGenome,
        ..Default::default()
    };
    let stats = build_index_from_files(&fa_path, &gtf_path, &idx_path, &cfg).unwrap();
    assert_eq!(stats.num_chroms, 1);
    assert_eq!(stats.num_genes, 5);

    let idx = MmapIndex::open(&idx_path).unwrap();
    assert_eq!(idx.chroms().len(), 1);
    assert_eq!(idx.genes().len(), 5);
    assert_eq!(idx.header().flags & FLAG_HASH_MAY_COLLIDE, 0);

    // Test O(log G + k) genomic interval query on chr0
    let mut hits = Vec::new();
    idx.query_overlapping_genes(0, 70, 95, &mut hits);
    let hit_names: Vec<&str> = hits
        .iter()
        .map(|&g| idx.get_str(idx.genes()[g as usize].gene_id))
        .collect();
    assert!(hit_names.contains(&"one"));
    assert!(hit_names.contains(&"spl_neg"));
    assert!(!hit_names.contains(&"single"));

    // Test auxiliary TSV outputs (.exons.tsv.gz and .paths.tsv.gz)
    write_exons_tsv_gz(&exons_path, &idx).unwrap();
    write_paths_tsv_gz(&paths_path, &idx).unwrap();

    let mut paths_txt = String::new();
    MultiGzDecoder::new(std::fs::File::open(&paths_path).unwrap())
        .read_to_string(&mut paths_txt)
        .unwrap();
    assert!(paths_txt.contains("one\tdef\t1-3-5-7\t1-3,3-5,5-7\n"));
    assert!(paths_txt.contains("one\tint1_alt3\t1-2-3-4-5-7\t"));
    assert!(paths_txt.contains("one\tapa_alt5\t1-3-5-6-7-8\t"));

    // Corrupt a byte in the header and verify CRC32 detects it cleanly
    let mut raw = std::fs::read(&idx_path).unwrap();
    raw[20] ^= 0xFF;
    let bad_path = dir.path().join("corrupt.lidx");
    std::fs::write(&bad_path, &raw).unwrap();
    assert!(MmapIndex::open(&bad_path).is_err());
}

#[test]
fn test_k_gt_27_sets_collision_guard_flag_and_preserves_all_seeds() {
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("k29.lidx");
    let seq = b"GCTAGCTTAGCCATGCTAGCATCGATCGATGCTAGCTAGGGCTAGCTTAGCCATGCTAGCATCGATCG";
    let gtf = "chr1\tT\texon\t1\t68\t.\t+\t.\tgene_id \"G29\"; transcript_id \"T29\";\n";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();

    let cfg = IndexBuildConfig {
        minimizer_k: 29,
        minimizer_w: 5,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", seq)], &refset, &idx_path, &cfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();
    assert_ne!(
        idx.header().flags & FLAG_HASH_MAY_COLLIDE,
        0,
        "k=29 (> 27) must set FLAG_HASH_MAY_COLLIDE"
    );

    let mut mins = Vec::new();
    extract_minimizers(seq, 29, 5, idx.hash_bits(), &mut mins);
    assert!(!mins.is_empty());
    for m in &mins {
        assert!(idx.lookup_minimizer(m.hash).is_some());
    }
}

#[test]
fn test_unsorted_interleaved_gtf_and_basic_tsl_filtering() {
    // Interleaved exon lines (Whippet2 crashes on this; Langur handles cleanly!)
    let gtf = "\
chr1\tG\texon\t10\t20\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; tag \"basic\"; transcript_support_level \"1\";
chr1\tG\texon\t10\t25\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\"; transcript_support_level \"4\";
chr1\tG\texon\t40\t50\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; tag \"basic\"; transcript_support_level \"1\";
chr1\tG\texon\t40\t60\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\"; transcript_support_level \"4\";
";
    // 1. suppress_low_tsl drops T2 (TSL 4)
    let r_tsl = load_gtf_reader(
        BufReader::new(Cursor::new(gtf)),
        &GtfOptions {
            basic_or_ss_only: false,
            suppress_low_tsl: true,
        },
    )
    .unwrap();
    let g1 = r_tsl.get("G1").unwrap();
    assert_eq!(g1.reftx.len(), 1);
    assert_eq!(g1.don, vec![20]);

    // 2. basic_or_ss_only keeps only T1 in reftx/txen, but includes T2's donor 25!
    let r_basic = load_gtf_reader(
        BufReader::new(Cursor::new(gtf)),
        &GtfOptions {
            basic_or_ss_only: true,
            suppress_low_tsl: false,
        },
    )
    .unwrap();
    let g1b = r_basic.get("G1").unwrap();
    assert_eq!(g1b.reftx.len(), 1);
    assert_eq!(g1b.don, vec![20, 25]);
    assert_eq!(g1b.txen, vec![50]);
}

#[test]
fn test_short_src_node_anno_junction_left_overhang_exact() {
    // k = 9 (flank = 8), but src_node (Exon 1) is only 5 bp (< k - 1 = 8)!
    // Exon 1 = "GCTAG" (5 bp, 10..14), Intron = 15..59, Exon 2 = "CCAATTGGCCGGAATT" (16 bp, 60..75).
    // Bridge sequence = "GCTAG" (left_len = 5) + "CCAATTGG" (8 bp) = "GCTAGCCAATTGG" (13 bp).
    // Any 9-mer straddling the junction starts at pos in 0..=4 and has left_overhang = 5 - pos in 1..=5.
    let k = 9usize;
    let w = 2usize;
    let mut chr = vec![b'A'; 120];
    chr[9..14].copy_from_slice(b"GCTAG");
    chr[59..75].copy_from_slice(b"CCAATTGGCCGGAATT");

    let gtf = "\
chr1\tT\texon\t10\t14\t.\t+\t.\tgene_id \"short_src\"; transcript_id \"t1\";
chr1\tT\texon\t60\t75\t.\t+\t.\tgene_id \"short_src\"; transcript_id \"t1\";
";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("short_jnc.lidx");
    let cfg = IndexBuildConfig {
        minimizer_k: k,
        minimizer_w: w,
        include_jnc_seeds: true,
        decoy_mode: DecoyMode::None,
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &cfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    // Synthesize a read spanning the exact junction: Exon1 (5 bp) + Exon2 (16 bp) = 21 bp
    // Junction lies right at index 5 in the read (bases 0..5 = donor, 5..21 = acceptor).
    let read = b"GCTAGCCAATTGGCCGGAATT";
    let mut mins = Vec::new();
    extract_minimizers(read, k, w, idx.hash_bits(), &mut mins);

    let mut verified_jnc_hit = false;
    for m in &mins {
        let Some(lookup) = idx.lookup_minimizer(m.hash) else {
            continue;
        };
        let hits: Vec<_> = match lookup {
            MinimizerLookup::Inline(h) => vec![h],
            MinimizerLookup::Postings { csg, .. } => csg.to_vec(),
        };
        for h in hits {
            if h.kind() == SeedHitKind::AnnoJunction {
                let (gid, jidx, left_overhang, hit_rc) = h.unpack_anno_junction();
                assert_eq!(gid, 0);
                assert_eq!(jidx, 0);
                assert!(left_overhang >= 1 && (left_overhang as usize) < k);
                // Same-strand hits satisfy `seed_pos + left_overhang == 5`!
                if m.is_rc == hit_rc {
                    assert_eq!(
                        m.pos + (left_overhang as u32),
                        5,
                        "seed_pos ({}) + left_overhang ({}) must equal true junction pos 5",
                        m.pos,
                        left_overhang
                    );
                    verified_jnc_hit = true;
                }
            }
        }
    }
    assert!(
        verified_jnc_hit,
        "expected at least one same-strand AnnoJunction minimizer hit on the junction-spanning read"
    );

    // Also verify reverse-complemented read: junction lies at index 16 (21 - 5),
    // and opposite-strand hits satisfy `seed_pos + (k - left_overhang) == 16`!
    let rc_read = revcomp(read);
    extract_minimizers(&rc_read, k, w, idx.hash_bits(), &mut mins);
    let mut verified_rc_jnc_hit = false;
    for m in &mins {
        let Some(lookup) = idx.lookup_minimizer(m.hash) else {
            continue;
        };
        let hits: Vec<_> = match lookup {
            MinimizerLookup::Inline(h) => vec![h],
            MinimizerLookup::Postings { csg, .. } => csg.to_vec(),
        };
        for h in hits {
            if h.kind() == SeedHitKind::AnnoJunction {
                let (_gid, _jidx, left_overhang, hit_rc) = h.unpack_anno_junction();
                if m.is_rc != hit_rc {
                    assert_eq!(m.pos + ((k as u32) - (left_overhang as u32)), 16);
                    verified_rc_jnc_hit = true;
                }
            }
        }
    }
    assert!(
        verified_rc_jnc_hit,
        "expected at least one opposite-strand AnnoJunction hit on RC read"
    );
}

#[test]
fn test_decoy_splits_around_poly_n_scaffold_gaps() {
    // Gene exon has distinct 30 bp sequence E1; intergenic region has E1 flanked immediately by poly-N!
    // Verify decoy contig is split cleanly at the 'N' run and never encodes 'N' as 'A'/'T'.
    let k = 9usize;
    let w = 4usize;
    let e1 = b"GCTAGCTTAGCCATGCTAGCATCGATCGAT"; // 30 bp
    let mut chr = vec![b'C'; 400];
    chr[20..50].copy_from_slice(e1);
    // Intergenic region 180..320: poly-N before and after an exact copy of E1
    for b in &mut chr[180..220] {
        *b = b'N';
    }
    chr[220..250].copy_from_slice(e1);
    for b in &mut chr[250..290] {
        *b = b'N';
    }

    let gtf = "chr1\tT\texon\t21\t50\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n";
    let refset = load_gtf_reader(BufReader::new(Cursor::new(gtf)), &GtfOptions::default()).unwrap();
    let dir = tempdir().unwrap();
    let idx_path = dir.path().join("poly_n_decoy.lidx");
    let cfg = IndexBuildConfig {
        minimizer_k: k,
        minimizer_w: w,
        decoy_mode: DecoyMode::SelectiveGenome,
        decoy_window_pad: 60, // would cross into both poly-N runs if not split!
        ..Default::default()
    };
    build_index_from_memory(&[("chr1", &chr)], &refset, &idx_path, &cfg).unwrap();
    let idx = MmapIndex::open(&idx_path).unwrap();

    assert_eq!(idx.decoys().len(), 1);
    let d = &idx.decoys()[0];
    // Must be bounded strictly to 1-based [221..=250] (length 30), excluding all N's!
    assert_eq!(d.genomic_start, 221);
    assert_eq!(d.seq_len, 30);
    let mut unpacked = Vec::new();
    idx.unpack_decoy_slice(0, 0, d.seq_len as usize, &mut unpacked);
    assert_eq!(unpacked, e1);
}

proptest! {
    #[test]
    fn prop_mix54_is_bijective_and_collision_free_for_k_le_27(
        a in 0u64..(1u64 << 54),
        b in 0u64..(1u64 << 54),
    ) {
        prop_assume!(a != b);
        prop_assert!(is_bijective_hash(27, 54));
        prop_assert!(!is_bijective_hash(29, 54));
        prop_assert_ne!(mix_bits(a, 54), mix_bits(b, 54));
        prop_assert_ne!(hash_kmer(a, 27, 54), hash_kmer(b, 27, 54));
    }
}
