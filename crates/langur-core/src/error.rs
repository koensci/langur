use std::path::PathBuf;
use thiserror::Error;

pub type Result<T> = std::result::Result<T, LangurError>;

#[derive(Debug, Error)]
pub enum LangurError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("I/O error at {path}: {source}")]
    IoPath {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("GTF parse error on line {line}: {message}")]
    GtfParse { line: usize, message: String },

    #[error("BAM error: {0}")]
    Bam(String),

    #[error("FASTA error: {0}")]
    Fasta(String),

    #[error("Invalid index parameter: {0}")]
    InvalidParam(String),

    #[error("Invalid or corrupted .lidx file: {0}")]
    InvalidIndex(String),

    #[error("Gene {gene_id} exceeds maximum supported node count ({count} > {max})")]
    TooManyNodes {
        gene_id: String,
        count: usize,
        max: usize,
    },

    #[error("No genes in GTF matched chromosome names in FASTA")]
    NoMatchingGenes,
}
