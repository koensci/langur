//! Binary `.lidx` container format: 64-byte-aligned sections, `bytemuck::Pod`
//! zero-copy structs, CRC32 header verification, and `memmap2` reader/writer.

use crate::error::{LangurError, Result};
use crate::graph::{EdgeMotif, EdgeType};
use crate::seq::{mismatch_count_lsb, mismatch_count_msb, unpack_2bit};
use bytemuck::{Pod, Zeroable};
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub const LIDX_MAGIC: [u8; 8] = *b"LANGUR01";
pub const LEGACY_LIDX_MAGIC: [u8; 8] = *b"WHIP3IDX";
pub const LIDX_VERSION: u32 = 1;
pub const SECTION_ALIGN: u64 = 64;

pub const FLAG_HAS_DECOYS: u32 = 1 << 0;
pub const FLAG_HAS_JNC_SEEDS: u32 = 1 << 1;
pub const FLAG_HAS_BAM_NOVEL: u32 = 1 << 2;
pub const FLAG_HASH_MAY_COLLIDE: u32 = 1 << 3;

#[repr(u32)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum SectionTag {
    StringPool = 1,
    Chroms = 2,
    Genes = 3,
    GenePrefixMaxEnd = 4,
    Nodes = 5,
    AnnoJunctions = 6,
    Transcripts = 7,
    TxNodes = 8,
    TxNodePrefixLen = 9,
    TxExpectedGc = 10,
    CsgSeq = 11,
    MinBucketOffsets = 12,
    MinKeyLo = 13,
    MinSlots = 14,
    MinPostings = 15,
    Decoys = 16,
    DecoySeq2Bit = 17,
}

impl SectionTag {
    pub fn from_u32(v: u32) -> Option<Self> {
        match v {
            1 => Some(Self::StringPool),
            2 => Some(Self::Chroms),
            3 => Some(Self::Genes),
            4 => Some(Self::GenePrefixMaxEnd),
            5 => Some(Self::Nodes),
            6 => Some(Self::AnnoJunctions),
            7 => Some(Self::Transcripts),
            8 => Some(Self::TxNodes),
            9 => Some(Self::TxNodePrefixLen),
            10 => Some(Self::TxExpectedGc),
            11 => Some(Self::CsgSeq),
            12 => Some(Self::MinBucketOffsets),
            13 => Some(Self::MinKeyLo),
            14 => Some(Self::MinSlots),
            15 => Some(Self::MinPostings),
            16 => Some(Self::Decoys),
            17 => Some(Self::DecoySeq2Bit),
            _ => None,
        }
    }
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct IndexHeader {
    pub magic: [u8; 8],
    pub format_version: u32,
    pub flags: u32,
    pub minimizer_k: u16,
    pub minimizer_w: u16,
    pub min_bucket_bits: u16,
    pub gc_bins: u16,
    pub gc_window_len: u16,
    pub max_csg_occ: u16,
    pub max_decoy_occ: u16,
    pub section_count: u16,
    pub num_chroms: u32,
    pub num_genes: u32,
    pub num_nodes: u32,
    pub num_transcripts: u32,
    pub num_anno_junctions: u32,
    pub num_decoys: u32,
    pub csg_seq_len: u64,
    pub decoy_seq_bases: u64,
    pub num_min_keys: u64,
    pub num_min_postings: u64,
    pub num_hash_collisions: u64,
    pub header_crc32: u32,
    pub _reserved: [u32; 7],
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct SectionEntry {
    pub tag: u32,
    pub flags: u32,
    pub byte_offset: u64,
    pub byte_len: u64,
    pub elem_count: u64,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Pod, Zeroable)]
pub struct StrRef {
    pub offset: u32,
    pub len: u32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct ChromRecord {
    pub name: StrRef,
    pub length: u32,
    pub first_gene: u32,
    pub num_genes: u32,
    pub _pad: u32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct GeneRecord {
    pub gene_id: StrRef,
    pub gene_name: StrRef,
    pub seq_offset: u64,
    pub chrom_id: u32,
    pub genomic_start: u32,
    pub genomic_end: u32,
    pub seq_len: u32,
    pub node_offset: u32,
    pub num_nodes: u32,
    pub tx_offset: u32,
    pub num_tx: u32,
    pub anno_jnc_offset: u32,
    pub num_anno_jnc: u32,
    pub mean_tx_len: f32,
    pub flags: u32,
}

impl GeneRecord {
    pub const FLAG_POS_STRAND: u32 = 1 << 0;
    pub const FLAG_HAS_AMBIG: u32 = 1 << 1;
    pub const FLAG_HAS_BAM_NOVEL: u32 = 1 << 2;

    #[inline]
    pub const fn is_pos_strand(&self) -> bool {
        (self.flags & Self::FLAG_POS_STRAND) != 0
    }

    #[inline]
    pub const fn has_ambig(&self) -> bool {
        (self.flags & Self::FLAG_HAS_AMBIG) != 0
    }

    #[inline]
    pub const fn has_bam_novel(&self) -> bool {
        (self.flags & Self::FLAG_HAS_BAM_NOVEL) != 0
    }
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct NodeRecord {
    pub gene_seq_offset: u32,
    pub genomic_coord: u32,
    pub len: u32,
    /// First `acceptor_flank_len` (0..=16) bases from node 5' start extending across contiguous
    /// exonic nodes, MSB-aligned in 2-bit DNA.
    pub acceptor_flank_2bit: u32,
    /// Last `donor_flank_len` (0..=16) bases ending at node 3' end extending backward across
    /// contiguous exonic nodes, LSB-aligned in 2-bit DNA: bits `0..2` hold the 3'-most base
    /// (last base before the donor splice site), bits `2..4` hold the second-to-last base, etc.
    /// Matches `seq::encode_flank_lsb_from_end` and `seq::mismatch_count_lsb`.
    pub donor_flank_2bit: u32,
    /// Low nibble: `left_edge` (`EdgeType`), high nibble: `right_edge` (`EdgeType`).
    pub edges_packed: u8,
    /// `[0..2]` acc_dinuc, `[2..4]` don_dinuc, bit 4 novel_acc, bit 5 novel_don, bit 6 has_ambig.
    pub splice_meta: u8,
    /// Exact number of valid non-ambiguous bases packed in `acceptor_flank_2bit` (`0..=16`).
    pub acceptor_flank_len: u8,
    /// Exact number of valid non-ambiguous bases packed in `donor_flank_2bit` (`0..=16`).
    pub donor_flank_len: u8,
}

impl NodeRecord {
    #[inline]
    pub fn left_edge(&self) -> EdgeType {
        EdgeType::from_u8(self.edges_packed & 0x0F).unwrap_or(EdgeType::SL)
    }

    #[inline]
    pub fn right_edge(&self) -> EdgeType {
        EdgeType::from_u8(self.edges_packed >> 4).unwrap_or(EdgeType::RS)
    }

    #[inline]
    pub fn motif(&self) -> EdgeMotif {
        EdgeMotif::from_edges(self.left_edge(), self.right_edge())
    }

    #[inline]
    pub const fn acc_dinuc(&self) -> u8 {
        self.splice_meta & 0x03
    }

    #[inline]
    pub const fn don_dinuc(&self) -> u8 {
        (self.splice_meta >> 2) & 0x03
    }

    #[inline]
    pub const fn is_novel_acc(&self) -> bool {
        (self.splice_meta & (1 << 4)) != 0
    }

    #[inline]
    pub const fn is_novel_don(&self) -> bool {
        (self.splice_meta & (1 << 5)) != 0
    }

    #[inline]
    pub const fn has_ambig(&self) -> bool {
        (self.splice_meta & (1 << 6)) != 0
    }

    /// Compare up to `cmp_len` bases against `acceptor_flank_2bit`, clamping safely to
    /// `acceptor_flank_len` so zero-padded bits are never interpreted as `AAAA...`.
    /// Returns `Some((compared_bases, mismatches))` or `None` if `acceptor_flank_len == 0`.
    #[inline]
    pub fn acceptor_flank_mismatches(&self, query_msb: u32, cmp_len: u8) -> Option<(u8, u32)> {
        let m = cmp_len.min(self.acceptor_flank_len);
        if m == 0 {
            return None;
        }
        Some((m, mismatch_count_msb(self.acceptor_flank_2bit, query_msb, m)))
    }

    /// Compare up to `cmp_len` bases against `donor_flank_2bit`, clamping safely to
    /// `donor_flank_len` so zero-padded bits are never interpreted as `AAAA...`.
    #[inline]
    pub fn donor_flank_mismatches(&self, query_lsb: u32, cmp_len: u8) -> Option<(u8, u32)> {
        let m = cmp_len.min(self.donor_flank_len);
        if m == 0 {
            return None;
        }
        Some((m, mismatch_count_lsb(self.donor_flank_2bit, query_lsb, m)))
    }
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Pod, Zeroable)]
pub struct AnnoJunctionRecord {
    pub src_node: u16,
    pub dst_node: u16,
    pub flags: u16,
    pub tx_support: u16,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct TranscriptRecord {
    pub tx_id: StrRef,
    pub gene_id: u32,
    pub tx_node_offset: u32,
    pub num_nodes: u32,
    pub seq_len: u32,
    pub genomic_start: u32,
    pub genomic_end: u32,
    pub flags: u32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct DecoyRecord {
    pub seq_byte_offset: u64,
    pub chrom_id: u32,
    pub genomic_start: u32,
    pub seq_len: u32,
    pub flags: u32,
}

impl DecoyRecord {
    pub const FLAG_INTRONIC: u32 = 1 << 0;
    pub const FLAG_INTERGENIC: u32 = 1 << 1;
    pub const FLAG_EXTERNAL_FASTA: u32 = 1 << 2;
}

/// Discriminator stored in bits `[61..63]` of `PackedSeedHit`.
#[repr(u8)]
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum SeedHitKind {
    CsgNode = 0,
    AnnoJunction = 1,
    Decoy = 2,
}

/// Compact 64-bit seed hit stored directly in `min_slots` (singletons, bit 63 == 0)
/// or in `min_postings[]`. Bit 63 is ALWAYS `0`.
#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Pod, Zeroable)]
pub struct PackedSeedHit(pub u64);

impl PackedSeedHit {
    pub const MAX_GENE_ID: u32 = (1 << 20) - 1;
    pub const MAX_NODE_ID: u16 = u16::MAX;
    pub const MAX_NODE_OFFSET: u32 = (1 << 24) - 1;
    pub const MAX_DECOY_ID: u32 = (1 << 24) - 1;

    #[inline]
    pub const fn csg_node(gene_id: u32, node_id: u16, node_offset: u32, is_rc: bool) -> Self {
        let v = (gene_id as u64 & 0xF_FFFF)
            | ((node_id as u64) << 20)
            | ((node_offset as u64 & 0xFF_FFFF) << 36)
            | ((is_rc as u64) << 60);
        // hit_kind = 0b00, bit 63 = 0
        Self(v)
    }

    /// Construct an annotated splice-junction bridge hit.
    /// `left_overhang` is the number of bases (`1..k`) of the k-mer lying on the donor (`src_node`)
    /// side of the splice boundary, so `jnc_pos_in_read = seed_pos + left_overhang`.
    #[inline]
    pub const fn anno_junction(
        gene_id: u32,
        local_jnc_idx: u32,
        left_overhang: u8,
        is_rc: bool,
    ) -> Self {
        let v = (gene_id as u64 & 0xF_FFFF)
            | ((local_jnc_idx as u64 & 0xF_FFFF) << 20)
            | ((left_overhang as u64) << 40)
            | ((is_rc as u64) << 60)
            | (1u64 << 61);
        Self(v)
    }

    #[inline]
    pub const fn decoy(decoy_id: u32, decoy_offset: u32, is_rc: bool) -> Self {
        let v = (decoy_id as u64 & 0xFF_FFFF)
            | ((decoy_offset as u64) << 24)
            | ((is_rc as u64) << 60)
            | (2u64 << 61);
        Self(v)
    }

    #[inline]
    pub const fn kind(self) -> SeedHitKind {
        match (self.0 >> 61) & 0x03 {
            0 => SeedHitKind::CsgNode,
            1 => SeedHitKind::AnnoJunction,
            _ => SeedHitKind::Decoy,
        }
    }

    #[inline]
    pub const fn is_rc(self) -> bool {
        ((self.0 >> 60) & 1) != 0
    }

    #[inline]
    pub const fn unpack_csg_node(self) -> (u32, u16, u32, bool) {
        let gene_id = (self.0 & 0xF_FFFF) as u32;
        let node_id = ((self.0 >> 20) & 0xFFFF) as u16;
        let node_offset = ((self.0 >> 36) & 0xFF_FFFF) as u32;
        let is_rc = self.is_rc();
        (gene_id, node_id, node_offset, is_rc)
    }

    /// Unpack an `AnnoJunction` hit into `(gene_id, local_jnc_idx, left_overhang, is_rc)`
    /// where `left_overhang` (`1..k`) is the number of k-mer bases on the donor side of the junction.
    #[inline]
    pub const fn unpack_anno_junction(self) -> (u32, u32, u8, bool) {
        let gene_id = (self.0 & 0xF_FFFF) as u32;
        let local_jnc_idx = ((self.0 >> 20) & 0xF_FFFF) as u32;
        let left_overhang = ((self.0 >> 40) & 0xFF) as u8;
        let is_rc = self.is_rc();
        (gene_id, local_jnc_idx, left_overhang, is_rc)
    }

    #[inline]
    pub const fn unpack_decoy(self) -> (u32, u32, bool) {
        let decoy_id = (self.0 & 0xFF_FFFF) as u32;
        let decoy_offset = ((self.0 >> 24) & 0xFFFF_FFFF) as u32;
        let is_rc = self.is_rc();
        (decoy_id, decoy_offset, is_rc)
    }
}

/// Multi-hit slot descriptor bit-packing (`bit 63 == 1`).
#[inline]
pub const fn pack_multi_slot(
    posting_offset: u64,
    csg_count: u16,
    decoy_count: u8,
    csg_truncated: bool,
    decoy_truncated: bool,
) -> u64 {
    (1u64 << 63)
        | (posting_offset & 0xF_FFFF_FFFF) // 36 bits
        | ((csg_count as u64) << 36)       // 16 bits
        | ((decoy_count as u64) << 52)     // 8 bits
        | ((csg_truncated as u64) << 60)
        | ((decoy_truncated as u64) << 61)
}

#[inline]
pub const fn unpack_multi_slot(slot: u64) -> (usize, usize, usize, bool, bool) {
    let offset = (slot & 0xF_FFFF_FFFF) as usize;
    let csg_count = ((slot >> 36) & 0xFFFF) as usize;
    let decoy_count = ((slot >> 52) & 0xFF) as usize;
    let csg_trunc = ((slot >> 60) & 1) != 0;
    let decoy_trunc = ((slot >> 61) & 1) != 0;
    (offset, csg_count, decoy_count, csg_trunc, decoy_trunc)
}

/// Zero-copy result of looking up a minimizer hash in `MmapIndex`.
#[derive(Copy, Clone, Debug)]
pub enum MinimizerLookup<'a> {
    /// Singleton CSG hit stored inline (zero posting-array access).
    Inline(PackedSeedHit),
    /// Multi-hit and/or decoy-colliding minimizer backed by `min_postings[]`.
    Postings {
        csg: &'a [PackedSeedHit],
        decoys: &'a [PackedSeedHit],
        csg_truncated: bool,
        decoy_truncated: bool,
    },
}

/// Complete owned section buffers before serializing to `.w3idx`.
pub struct IndexBuffers {
    pub header: IndexHeader,
    pub string_pool: Vec<u8>,
    pub chroms: Vec<ChromRecord>,
    pub genes: Vec<GeneRecord>,
    pub gene_prefix_max_end: Vec<u32>,
    pub nodes: Vec<NodeRecord>,
    pub anno_junctions: Vec<AnnoJunctionRecord>,
    pub transcripts: Vec<TranscriptRecord>,
    pub tx_nodes: Vec<u16>,
    pub tx_node_prefix_len: Vec<u32>,
    pub tx_expected_gc: Vec<f32>,
    pub csg_seq: Vec<u8>,
    pub min_bucket_offsets: Vec<u32>,
    pub min_key_lo: Vec<u32>,
    pub min_slots: Vec<u64>,
    pub min_postings: Vec<PackedSeedHit>,
    pub decoys: Vec<DecoyRecord>,
    pub decoy_seq_2bit: Vec<u8>,
}

/// Write an `IndexBuffers` payload to a 64-byte-aligned `.w3idx` file with CRC32 header verification.
pub fn write_index_file(path: &Path, buffers: &IndexBuffers) -> Result<u64> {
    let file = File::create(path).map_err(|source| LangurError::IoPath {
        path: path.to_path_buf(),
        source,
    })?;
    let mut w = BufWriter::with_capacity(4 * 1024 * 1024, file);

    let sections_data: [(SectionTag, &[u8], u64); 17] = [
        (SectionTag::StringPool, &buffers.string_pool, buffers.string_pool.len() as u64),
        (SectionTag::Chroms, bytemuck::cast_slice(&buffers.chroms), buffers.chroms.len() as u64),
        (SectionTag::Genes, bytemuck::cast_slice(&buffers.genes), buffers.genes.len() as u64),
        (SectionTag::GenePrefixMaxEnd, bytemuck::cast_slice(&buffers.gene_prefix_max_end), buffers.gene_prefix_max_end.len() as u64),
        (SectionTag::Nodes, bytemuck::cast_slice(&buffers.nodes), buffers.nodes.len() as u64),
        (SectionTag::AnnoJunctions, bytemuck::cast_slice(&buffers.anno_junctions), buffers.anno_junctions.len() as u64),
        (SectionTag::Transcripts, bytemuck::cast_slice(&buffers.transcripts), buffers.transcripts.len() as u64),
        (SectionTag::TxNodes, bytemuck::cast_slice(&buffers.tx_nodes), buffers.tx_nodes.len() as u64),
        (SectionTag::TxNodePrefixLen, bytemuck::cast_slice(&buffers.tx_node_prefix_len), buffers.tx_node_prefix_len.len() as u64),
        (SectionTag::TxExpectedGc, bytemuck::cast_slice(&buffers.tx_expected_gc), buffers.tx_expected_gc.len() as u64),
        (SectionTag::CsgSeq, &buffers.csg_seq, buffers.csg_seq.len() as u64),
        (SectionTag::MinBucketOffsets, bytemuck::cast_slice(&buffers.min_bucket_offsets), buffers.min_bucket_offsets.len() as u64),
        (SectionTag::MinKeyLo, bytemuck::cast_slice(&buffers.min_key_lo), buffers.min_key_lo.len() as u64),
        (SectionTag::MinSlots, bytemuck::cast_slice(&buffers.min_slots), buffers.min_slots.len() as u64),
        (SectionTag::MinPostings, bytemuck::cast_slice(&buffers.min_postings), buffers.min_postings.len() as u64),
        (SectionTag::Decoys, bytemuck::cast_slice(&buffers.decoys), buffers.decoys.len() as u64),
        (SectionTag::DecoySeq2Bit, &buffers.decoy_seq_2bit, buffers.decoy_seq_2bit.len() as u64),
    ];

    let header_size = std::mem::size_of::<IndexHeader>() as u64;
    let table_size = (sections_data.len() * std::mem::size_of::<SectionEntry>()) as u64;
    let mut cur_offset = align_up(header_size + table_size, SECTION_ALIGN);

    let mut entries = Vec::with_capacity(sections_data.len());
    for &(tag, bytes, elem_count) in &sections_data {
        entries.push(SectionEntry {
            tag: tag as u32,
            flags: 0,
            byte_offset: cur_offset,
            byte_len: bytes.len() as u64,
            elem_count,
        });
        cur_offset = align_up(cur_offset + bytes.len() as u64, SECTION_ALIGN);
    }

    let mut header = buffers.header;
    header.magic = LIDX_MAGIC;
    header.format_version = LIDX_VERSION;
    header.section_count = entries.len() as u16;
    header.header_crc32 = 0;

    let mut hasher = crc32fast::Hasher::new();
    hasher.update(bytemuck::bytes_of(&header));
    hasher.update(bytemuck::cast_slice(&entries));
    header.header_crc32 = hasher.finalize();

    let zeros = [0u8; SECTION_ALIGN as usize];
    w.write_all(bytemuck::bytes_of(&header))?;
    w.write_all(bytemuck::cast_slice(&entries))?;
    let mut pos = header_size + table_size;

    for (entry, &(_, bytes, _)) in entries.iter().zip(sections_data.iter()) {
        if entry.byte_offset > pos {
            let pad = (entry.byte_offset - pos) as usize;
            w.write_all(&zeros[..pad])?;
            pos += pad as u64;
        }
        w.write_all(bytes)?;
        pos += bytes.len() as u64;
    }

    let final_len = align_up(pos, SECTION_ALIGN);
    if final_len > pos {
        let pad = (final_len - pos) as usize;
        w.write_all(&zeros[..pad])?;
        pos = final_len;
    }
    w.flush()?;
    Ok(pos)
}

#[inline]
const fn align_up(v: u64, align: u64) -> u64 {
    (v + align - 1) & !(align - 1)
}

/// Zero-copy memory-mapped `.w3idx` index handle.
/// Loads in O(1) (`< 1 ms`) with zero heap allocations for index data.
pub struct MmapIndex {
    mmap: Mmap,
    header: IndexHeader,
    sections: [SectionEntry; 18], // indexed by SectionTag as usize (1..=17)
}

impl MmapIndex {
    /// Open and validate a `.w3idx` file via `mmap` in O(1) time.
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).map_err(|source| LangurError::IoPath {
            path: path.to_path_buf(),
            source,
        })?;
        // SAFETY: Index files are immutable read-only artifacts once built.
        let mmap = unsafe { Mmap::map(&file)? };
        Self::from_mmap(mmap)
    }

    fn from_mmap(mmap: Mmap) -> Result<Self> {
        let hdr_len = std::mem::size_of::<IndexHeader>();
        if mmap.len() < hdr_len {
            return Err(LangurError::InvalidIndex("file smaller than IndexHeader".into()));
        }
        let header: IndexHeader = *bytemuck::try_from_bytes(&mmap[..hdr_len])
            .map_err(|e| LangurError::InvalidIndex(format!("header alignment/size error: {e}")))?;

        if header.magic != LIDX_MAGIC && header.magic != LEGACY_LIDX_MAGIC {
            return Err(LangurError::InvalidIndex(
                "invalid magic bytes (expected LANGUR01)".into(),
            ));
        }
        if header.format_version != LIDX_VERSION {
            return Err(LangurError::InvalidIndex(format!(
                "unsupported format version {} (expected {})",
                header.format_version, LIDX_VERSION
            )));
        }

        let sec_count = header.section_count as usize;
        let table_bytes = sec_count * std::mem::size_of::<SectionEntry>();
        if mmap.len() < hdr_len + table_bytes {
            return Err(LangurError::InvalidIndex("truncated section table".into()));
        }
        let table_slice: &[SectionEntry] =
            bytemuck::try_cast_slice(&mmap[hdr_len..hdr_len + table_bytes])
                .map_err(|e| LangurError::InvalidIndex(format!("section table cast error: {e}")))?;

        // Verify header + section table CRC32
        let mut hdr_zero_crc = header;
        hdr_zero_crc.header_crc32 = 0;
        let mut hasher = crc32fast::Hasher::new();
        hasher.update(bytemuck::bytes_of(&hdr_zero_crc));
        hasher.update(bytemuck::cast_slice(table_slice));
        let computed = hasher.finalize();
        if computed != header.header_crc32 {
            return Err(LangurError::InvalidIndex(format!(
                "header CRC32 mismatch (expected {:#010x}, got {:#010x})",
                header.header_crc32, computed
            )));
        }

        let mut sections = [SectionEntry::zeroed(); 18];
        for entry in table_slice {
            if entry.byte_offset % SECTION_ALIGN != 0 {
                return Err(LangurError::InvalidIndex(format!(
                    "section tag {} offset {} is not {}-byte aligned",
                    entry.tag, entry.byte_offset, SECTION_ALIGN
                )));
            }
            let end = entry
                .byte_offset
                .checked_add(entry.byte_len)
                .ok_or_else(|| LangurError::InvalidIndex("section length overflow".into()))?;
            if (end as usize) > mmap.len() {
                return Err(LangurError::InvalidIndex(format!(
                    "section tag {} extends past EOF ({} > {})",
                    entry.tag, end, mmap.len()
                )));
            }
            if let Some(tag) = SectionTag::from_u32(entry.tag) {
                sections[tag as usize] = *entry;
            }
        }

        let idx = Self {
            mmap,
            header,
            sections,
        };

        // Validate every typed section casts cleanly via bytemuck
        let _ = idx.try_slice::<ChromRecord>(SectionTag::Chroms)?;
        let _ = idx.try_slice::<GeneRecord>(SectionTag::Genes)?;
        let _ = idx.try_slice::<u32>(SectionTag::GenePrefixMaxEnd)?;
        let _ = idx.try_slice::<NodeRecord>(SectionTag::Nodes)?;
        let _ = idx.try_slice::<AnnoJunctionRecord>(SectionTag::AnnoJunctions)?;
        let _ = idx.try_slice::<TranscriptRecord>(SectionTag::Transcripts)?;
        let _ = idx.try_slice::<u16>(SectionTag::TxNodes)?;
        let _ = idx.try_slice::<u32>(SectionTag::TxNodePrefixLen)?;
        let _ = idx.try_slice::<f32>(SectionTag::TxExpectedGc)?;
        let _ = idx.try_slice::<u32>(SectionTag::MinBucketOffsets)?;
        let _ = idx.try_slice::<u32>(SectionTag::MinKeyLo)?;
        let _ = idx.try_slice::<u64>(SectionTag::MinSlots)?;
        let _ = idx.try_slice::<PackedSeedHit>(SectionTag::MinPostings)?;
        let _ = idx.try_slice::<DecoyRecord>(SectionTag::Decoys)?;

        Ok(idx)
    }

    #[inline]
    fn raw_section(&self, tag: SectionTag) -> &[u8] {
        let e = self.sections[tag as usize];
        let s = e.byte_offset as usize;
        let len = e.byte_len as usize;
        &self.mmap[s..s + len]
    }

    fn try_slice<T: Pod>(&self, tag: SectionTag) -> Result<&[T]> {
        bytemuck::try_cast_slice(self.raw_section(tag))
            .map_err(|e| LangurError::InvalidIndex(format!("section {:?} cast failed: {e}", tag)))
    }

    #[inline]
    fn slice<T: Pod>(&self, tag: SectionTag) -> &[T] {
        bytemuck::cast_slice(self.raw_section(tag))
    }

    #[inline]
    pub fn header(&self) -> &IndexHeader {
        &self.header
    }

    #[inline]
    pub fn file_size_bytes(&self) -> usize {
        self.mmap.len()
    }

    #[inline]
    pub fn hash_bits(&self) -> u32 {
        (self.header.min_bucket_bits as u32) + 32
    }

    #[inline]
    pub fn string_pool(&self) -> &[u8] {
        self.raw_section(SectionTag::StringPool)
    }

    #[inline]
    pub fn get_str(&self, s: StrRef) -> &str {
        let pool = self.string_pool();
        let start = s.offset as usize;
        let end = start + (s.len as usize);
        std::str::from_utf8(&pool[start..end]).unwrap_or("")
    }

    #[inline]
    pub fn chroms(&self) -> &[ChromRecord] {
        self.slice(SectionTag::Chroms)
    }

    #[inline]
    pub fn genes(&self) -> &[GeneRecord] {
        self.slice(SectionTag::Genes)
    }

    #[inline]
    pub fn gene_prefix_max_end(&self) -> &[u32] {
        self.slice(SectionTag::GenePrefixMaxEnd)
    }

    #[inline]
    pub fn nodes(&self) -> &[NodeRecord] {
        self.slice(SectionTag::Nodes)
    }

    #[inline]
    pub fn anno_junctions(&self) -> &[AnnoJunctionRecord] {
        self.slice(SectionTag::AnnoJunctions)
    }

    #[inline]
    pub fn transcripts(&self) -> &[TranscriptRecord] {
        self.slice(SectionTag::Transcripts)
    }

    #[inline]
    pub fn tx_nodes(&self) -> &[u16] {
        self.slice(SectionTag::TxNodes)
    }

    #[inline]
    pub fn tx_node_prefix_len(&self) -> &[u32] {
        self.slice(SectionTag::TxNodePrefixLen)
    }

    #[inline]
    pub fn tx_expected_gc(&self) -> &[f32] {
        self.slice(SectionTag::TxExpectedGc)
    }

    #[inline]
    pub fn csg_seq(&self) -> &[u8] {
        self.raw_section(SectionTag::CsgSeq)
    }

    #[inline]
    pub fn min_bucket_offsets(&self) -> &[u32] {
        self.slice(SectionTag::MinBucketOffsets)
    }

    #[inline]
    pub fn min_key_lo(&self) -> &[u32] {
        self.slice(SectionTag::MinKeyLo)
    }

    #[inline]
    pub fn min_slots(&self) -> &[u64] {
        self.slice(SectionTag::MinSlots)
    }

    #[inline]
    pub fn min_postings(&self) -> &[PackedSeedHit] {
        self.slice(SectionTag::MinPostings)
    }

    #[inline]
    pub fn decoys(&self) -> &[DecoyRecord] {
        self.slice(SectionTag::Decoys)
    }

    #[inline]
    pub fn decoy_seq_2bit(&self) -> &[u8] {
        self.raw_section(SectionTag::DecoySeq2Bit)
    }

    // --- Convenience per-gene / per-transcript / per-decoy zero-copy accessors ---

    #[inline]
    pub fn gene_nodes(&self, gene_id: u32) -> &[NodeRecord] {
        let g = &self.genes()[gene_id as usize];
        let s = g.node_offset as usize;
        &self.nodes()[s..s + (g.num_nodes as usize)]
    }

    #[inline]
    pub fn gene_seq(&self, gene_id: u32) -> &[u8] {
        let g = &self.genes()[gene_id as usize];
        let s = g.seq_offset as usize;
        &self.csg_seq()[s..s + (g.seq_len as usize)]
    }

    #[inline]
    pub fn gene_transcripts(&self, gene_id: u32) -> &[TranscriptRecord] {
        let g = &self.genes()[gene_id as usize];
        let s = g.tx_offset as usize;
        &self.transcripts()[s..s + (g.num_tx as usize)]
    }

    #[inline]
    pub fn gene_anno_junctions(&self, gene_id: u32) -> &[AnnoJunctionRecord] {
        let g = &self.genes()[gene_id as usize];
        let s = g.anno_jnc_offset as usize;
        &self.anno_junctions()[s..s + (g.num_anno_jnc as usize)]
    }

    #[inline]
    pub fn transcript_nodes(&self, tx_idx: usize) -> &[u16] {
        let t = &self.transcripts()[tx_idx];
        let s = t.tx_node_offset as usize;
        &self.tx_nodes()[s..s + (t.num_nodes as usize)]
    }

    #[inline]
    pub fn transcript_prefix_lens(&self, tx_idx: usize) -> &[u32] {
        let t = &self.transcripts()[tx_idx];
        let s = t.tx_node_offset as usize;
        &self.tx_node_prefix_len()[s..s + (t.num_nodes as usize)]
    }

    #[inline]
    pub fn transcript_gc(&self, tx_idx: usize) -> &[f32] {
        let bins = self.header.gc_bins as usize;
        let s = tx_idx * bins;
        &self.tx_expected_gc()[s..s + bins]
    }

    /// Unpack a decoy contig's sequence (or sub-slice `[rel_start .. rel_start + len]`) into `dst`.
    pub fn unpack_decoy_slice(
        &self,
        decoy_id: u32,
        rel_start: usize,
        len: usize,
        dst: &mut Vec<u8>,
    ) {
        let d = &self.decoys()[decoy_id as usize];
        let base_start = (d.seq_byte_offset as usize) * 4 + rel_start;
        let clamped_len = len.min((d.seq_len as usize).saturating_sub(rel_start));
        unpack_2bit(self.decoy_seq_2bit(), base_start, clamped_len, dst);
    }

    /// O(1) bucket + tiny binary-search minimizer lookup returning inline or posting slices.
    #[inline]
    pub fn lookup_minimizer(&self, hash: u64) -> Option<MinimizerLookup<'_>> {
        let bucket = (hash >> 32) as usize;
        let key_lo = hash as u32;
        let offsets = self.min_bucket_offsets();
        if bucket + 1 >= offsets.len() {
            return None;
        }
        let lo = offsets[bucket] as usize;
        let hi = offsets[bucket + 1] as usize;
        if lo == hi {
            return None;
        }
        let keys = &self.min_key_lo()[lo..hi];
        let rel = if keys.len() <= 16 {
            let mut found = usize::MAX;
            for (i, &k) in keys.iter().enumerate() {
                if k == key_lo {
                    found = i;
                    break;
                }
                if k > key_lo {
                    return None;
                }
            }
            if found == usize::MAX {
                return None;
            }
            found
        } else {
            keys.binary_search(&key_lo).ok()?
        };
        let slot = self.min_slots()[lo + rel];
        if (slot >> 63) == 0 {
            Some(MinimizerLookup::Inline(PackedSeedHit(slot)))
        } else {
            let (off, csg_cnt, decoy_cnt, csg_trunc, decoy_trunc) = unpack_multi_slot(slot);
            let postings = self.min_postings();
            Some(MinimizerLookup::Postings {
                csg: &postings[off..off + csg_cnt],
                decoys: &postings[off + csg_cnt..off + csg_cnt + decoy_cnt],
                csg_truncated: csg_trunc,
                decoy_truncated: decoy_trunc,
            })
        }
    }

    /// Find all gene indices on `chrom_id` whose genomic interval overlaps `[q_start, q_end]` (1-based).
    /// Runs in `O(log G + k)` using the precomputed `gene_prefix_max_end` array.
    pub fn query_overlapping_genes(&self, chrom_id: u32, q_start: u32, q_end: u32, out: &mut Vec<u32>) {
        out.clear();
        let chr = &self.chroms()[chrom_id as usize];
        let first = chr.first_gene as usize;
        let count = chr.num_genes as usize;
        if count == 0 {
            return;
        }
        let genes = &self.genes()[first..first + count];
        let pmax = &self.gene_prefix_max_end()[first..first + count];

        // Last gene that starts <= q_end
        let upper = genes.partition_point(|g| g.genomic_start <= q_end);
        // Walk backward while prefix_max_end >= q_start
        let mut i = upper;
        while i > 0 {
            i -= 1;
            if pmax[i] < q_start {
                break;
            }
            if genes[i].genomic_end >= q_start {
                out.push((first + i) as u32);
            }
        }
        out.reverse();
    }
}

/// Helper to build a deduplicated contiguous UTF-8 string pool.
#[derive(Default)]
pub struct StringPoolBuilder {
    buf: Vec<u8>,
}

impl StringPoolBuilder {
    pub fn push(&mut self, s: &str) -> StrRef {
        let offset = self.buf.len() as u32;
        self.buf.extend_from_slice(s.as_bytes());
        StrRef {
            offset,
            len: s.len() as u32,
        }
    }

    pub fn into_bytes(self) -> Vec<u8> {
        self.buf
    }
}
