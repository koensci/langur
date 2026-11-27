//! Quantification engine: deterministic sparse sharding, empirical fragment-length
//! modeling, unified 5'-hexamer + 20-bin GC bias (2-pass feedback), Union-Find
//! connected-component transcript VBEM + SQUAREM, and Stage 1 -> Stage 2 routing.

pub mod bias;
pub mod fld;
pub mod shard;
pub mod vbem;

pub use bias::{BiasAccum, BiasModel, FeatureBias};
pub use fld::{FragLenAccum, FragLenDist};
pub use shard::{
    JunctionStats, MultiAlignCandidate, MultiGroupKey, QuantShard, SampleCounts, SampleMapStats,
};
pub use vbem::{route_multimappers_to_stage2, solve_cc_vbem, VbemConfig, VbemResult};
