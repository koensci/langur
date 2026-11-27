//! GTF2.2 gene annotation loader and optional BAM splice-site supplementation.

pub mod bam;
pub mod gtf;

pub use bam::BamSupplementOptions;
pub use gtf::{load_gtf_file, load_gtf_reader, GtfOptions, RefGene, RefSet, RefTx};
