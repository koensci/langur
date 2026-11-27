//! Exact Whippet2 local splice-graph path enumeration (`reduce_graph`),
//! complexity (`K_n = ceil(log2(paths))`), and Shannon entropy, made 100%
//! deterministic across all event complexities.

use rustc_hash::FxHashSet;

/// An input edge `(src, dst)` with observed count and effective length.
#[derive(Copy, Clone, Debug)]
pub struct PathEdge {
    pub src: u16,
    pub dst: u16,
    pub count: f64,
    pub eff_len: f64,
}

/// A minimal enumerated path through local CSG nodes.
#[derive(Clone, Debug)]
pub struct ReducedPath {
    /// Strictly ascending 0-based local node IDs along this path.
    pub nodes: Vec<u16>,
    pub count: f64,
    pub eff_len: f64,
    pub psi: f64,
}

impl ReducedPath {
    #[inline]
    pub fn first(&self) -> u16 {
        self.nodes[0]
    }

    #[inline]
    pub fn last(&self) -> u16 {
        self.nodes[self.nodes.len() - 1]
    }

    /// Returns `true` iff `(src, dst)` are adjacent consecutive nodes in this path
    /// (matching Whippet2 `in(first, last, is::BitSet)`).
    #[inline]
    pub fn contains_edge(&self, src: u16, dst: u16) -> bool {
        if let Ok(pos) = self.nodes.binary_search(&src) {
            pos + 1 < self.nodes.len() && self.nodes[pos + 1] == dst
        } else {
            false
        }
    }
}

/// Returns `true` iff two monotonic paths share a terminal endpoint (`first(a) == last(b)` or `first(b) == last(a)`).
#[inline]
fn has_terminal_intersect(a: &[u16], b: &[u16]) -> bool {
    a[0] == b[b.len() - 1] || b[0] == a[a.len() - 1]
}

/// Merge two sorted unique node vectors into a sorted unique union.
fn union_sorted(a: &[u16], b: &[u16]) -> Vec<u16> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        if a[i] < b[j] {
            out.push(a[i]);
            i += 1;
        } else if b[j] < a[i] {
            out.push(b[j]);
            j += 1;
        } else {
            out.push(a[i]);
            i += 1;
            j += 1;
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

/// Build the minimal set of full paths explaining a local splice-graph edge set,
/// matching `Whippet2.jl/src/events.jl::reduce_graph` deterministically.
pub fn reduce_graph(edges: &[PathEdge], max_it: usize, max_iso: usize) -> Vec<ReducedPath> {
    if edges.is_empty() {
        return Vec::new();
    }

    let base_edges: Vec<ReducedPath> = edges
        .iter()
        .map(|e| ReducedPath {
            nodes: if e.src <= e.dst {
                vec![e.src, e.dst]
            } else {
                vec![e.dst, e.src]
            },
            count: e.count,
            eff_len: e.eff_len,
            psi: 0.0,
        })
        .collect();

    let mut cur = base_edges.clone();
    let mut next: Vec<ReducedPath> = Vec::new();
    let mut seen: FxHashSet<Vec<u16>> = FxHashSet::default();

    for it in 1..=max_it {
        next.clear();
        seen.clear();

        if cur.len() > max_iso {
            // Deterministic truncation by descending count, then lexicographic node path
            cur.sort_by(|a, b| {
                b.count
                    .partial_cmp(&a.count)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.nodes.cmp(&b.nodes))
            });
            cur.truncate(max_iso);
        }

        for gi in &cur {
            let mut extended = false;
            for ej in &base_edges {
                if has_terminal_intersect(&gi.nodes, &ej.nodes) {
                    extended = true;
                    let joint = union_sorted(&gi.nodes, &ej.nodes);
                    if seen.insert(joint.clone()) {
                        next.push(ReducedPath {
                            nodes: joint,
                            count: gi.count + ej.count,
                            eff_len: gi.eff_len + ej.eff_len,
                            psi: 0.0,
                        });
                    }
                }
            }
            if !extended && seen.insert(gi.nodes.clone()) {
                next.push(gi.clone());
            }
        }

        let unchanged = it > 1
            && next.len() == cur.len()
            && next.iter().zip(cur.iter()).all(|(a, b)| a.nodes == b.nodes);
        std::mem::swap(&mut cur, &mut next);
        if unchanged {
            break;
        }
    }

    cur.sort_by(|a, b| a.nodes.cmp(&b.nodes));
    cur
}

/// Whippet event complexity: `K_n = ceil(log2(num_paths))` (`0` when `num_paths <= 1`).
#[inline]
pub fn complexity(num_paths: usize) -> usize {
    if num_paths <= 1 {
        0
    } else {
        (num_paths as f64).log2().ceil() as usize
    }
}

/// Shannon entropy `-sum(p_i * log2(p_i))` in bits.
#[inline]
pub fn shannon_entropy(probs: &[f64]) -> f64 {
    let mut h = 0.0f64;
    for &p in probs {
        if p > 0.0 {
            h -= p * p.log2();
        }
    }
    if h.abs() < 1e-12 { 0.0 } else { h }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_edges(s: &str) -> Vec<PathEdge> {
        s.split(',')
            .map(|tok| {
                let mut parts = tok.split(['-', ':']);
                let u: u16 = parts.next().unwrap().parse().unwrap();
                let v: u16 = parts.next().unwrap().parse().unwrap();
                let c: f64 = parts.next().unwrap().parse().unwrap();
                PathEdge {
                    src: u,
                    dst: v,
                    count: c,
                    eff_len: 1.0,
                }
            })
            .collect()
    }

    #[test]
    fn test_reduce_graph_whippet2_fixtures() {
        // Exact test 1 from Whippet2.jl/test/runtests.jl L664-695 (ENSMUSG00000038685): 18 paths!
        let edgestr = "1-2:616.0,1-3:4.0,2-3:569.0,2-4:20.0,3-4:629.0,3-6:1.0,3-10:3.0,4-5:1.0,4-6:664.0,5-6:5.0,6-7:789.0,7-8:790.0,8-9:2.0,8-10:606.0,9-10:15.0";
        let res = reduce_graph(&parse_edges(edgestr), 64, 1024);
        assert_eq!(res.len(), 18);
        assert_eq!(complexity(res.len()), 5); // ceil(log2(18)) = 5 -> K5

        // Exact test 2 from Whippet2.jl/test/runtests.jl L697-710 (nrf2test): 6 paths!
        let nrf2 = "6-14:129.0,6-17:3.0,7-8:2.0,8-9:5.0,9-10:8.0,10-11:4.0,11-12:26.0,12-13:22.0,13-14:36.0,13-17:11.0,14-15:3.0,14-17:382.0,15-16:7.0";
        let res2 = reduce_graph(&parse_edges(nrf2), 64, 1024);
        assert_eq!(res2.len(), 6);
        let paths: Vec<Vec<u16>> = res2.into_iter().map(|p| p.nodes).collect();
        assert!(paths.contains(&vec![6, 14, 15, 16]));
        assert!(paths.contains(&vec![6, 14, 17]));
        assert!(paths.contains(&vec![6, 17]));
        assert!(paths.contains(&vec![7, 8, 9, 10, 11, 12, 13, 14, 15, 16]));
        assert!(paths.contains(&vec![7, 8, 9, 10, 11, 12, 13, 14, 17]));
        assert!(paths.contains(&vec![7, 8, 9, 10, 11, 12, 13, 17]));
    }
}
