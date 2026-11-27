//! Stage 2: Local splice-graph path reduction (`reduce_graph`), complexity (`K_n`),
//! Shannon entropy, PCR-stack soft-capping, local event PSI EM/VBEM, and
//! Beta/Dirichlet inferential uncertainty (`Inf_Var`).

pub mod gibbs;
pub mod psi;
pub mod reduce;

pub use gibbs::{sample_complex_event_gibbs, GibbsPosteriorSummary};
pub use psi::{quantify_all_events, EventConfig, EventRecord};
pub use reduce::{complexity, reduce_graph, shannon_entropy, PathEdge, ReducedPath};
