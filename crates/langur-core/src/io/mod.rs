//! High-throughput FASTQ (`zlib-ng`) streaming, parallel RFC 1952 multi-member
//! gzip output writers (`.psi.gz`, `.tpm.gz`, `.jnc.gz`, `.map.gz`), SAM output,
//! and end-to-end `quantify_sample` / `--batch` pipeline orchestration.

pub mod fastq;
pub mod pipeline;
pub mod writers;

pub use fastq::{FastqBatchReader, FastqChunk, FastqRec, PairedFastqBatchReader};
pub use pipeline::{quantify_batch_manifest, quantify_sample, QuantConfig, QuantSummary};
pub use writers::{
    write_gene_tpm_gz, write_isoform_tpm_gz, write_jnc_gz, write_map_gz, write_psi_gz,
};
