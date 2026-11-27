//! Langur core library: clean-room splice-graph construction, zero-copy `mmap`
//! minimizer + selective-decoy index, and downstream quantification/differential
//! splicing modules.

pub mod align;
pub mod annotation;
pub mod delta;
pub mod error;
pub mod events;
pub mod graph;
pub mod index;
pub mod io;
pub mod quant;
pub mod seq;

pub use error::{Result, LangurError};
