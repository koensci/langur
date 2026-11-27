//! Zero-copy `mmap` `.lidx` splice-graph + minimizer + selective-decoy index.

pub mod build;
pub mod decoy;
pub mod format;
pub mod minimizer;

pub use build::{
    build_index_from_files, build_index_from_memory, write_exons_tsv_gz, write_paths_tsv_gz,
    DecoyMode, IndexBuildConfig, IndexBuildStats,
};
pub use format::{
    AnnoJunctionRecord, ChromRecord, DecoyRecord, GeneRecord, IndexHeader, MinimizerLookup,
    MmapIndex, NodeRecord, PackedSeedHit, SectionEntry, SectionTag, SeedHitKind, StrRef,
    TranscriptRecord, FLAG_HAS_BAM_NOVEL, FLAG_HAS_DECOYS, FLAG_HAS_JNC_SEEDS,
    FLAG_HASH_MAY_COLLIDE, LEGACY_LIDX_MAGIC, LIDX_MAGIC, LIDX_VERSION,
};
