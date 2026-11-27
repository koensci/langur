//! Contiguous Splice Graph (CSG) construction, Whippet2-compatible edge/motif
//! taxonomy, per-gene contiguous splice flanks, and ExpectedGC modeling.

pub mod builder;
pub mod edge_type;
pub mod gc;

pub use builder::{build_splice_graph, BuiltAnnoJunction, BuiltNode, BuiltSpliceGraph, BuiltTx};
pub use edge_type::{EdgeMotif, EdgeType};
pub use gc::{expected_gc, GC_BINS, GC_WINDOW_LEN};
