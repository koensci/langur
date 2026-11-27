//! Expected GC content distribution per annotated transcript (`GCBiasParams(50, 0.05)`),
//! matching Whippet2's 20-bin (5% width) GC model in O(n) sliding-window time.

pub const GC_BINS: usize = 20;
pub const GC_WINDOW_LEN: usize = 50;

#[inline]
fn is_gc(b: u8) -> u32 {
    matches!(b, b'G' | b'g' | b'C' | b'c') as u32
}

/// Compute normalized 20-bin expected GC distribution over sliding windows of `window_len` bp.
/// Matches `Whippet2.jl/src/graph.jl::ExpectedGC` (`div(gc_content, 0.05)` clamped to `0..GC_BINS`).
pub fn expected_gc(seq: &[u8], window_len: usize) -> [f32; GC_BINS] {
    let mut bins = [0.0f32; GC_BINS];
    if window_len == 0 || seq.len() < window_len {
        return bins;
    }

    let mut gc_count: u32 = 0;
    for &b in &seq[..window_len] {
        gc_count += is_gc(b);
    }

    let n_windows = seq.len() - window_len + 1;
    let inv_win = 1.0f64 / (window_len as f64);

    for i in 0..n_windows {
        let frac = (gc_count as f64) * inv_win;
        // Bin index: floor(frac / 0.05), clamped to [0, GC_BINS - 1]
        let bin = ((frac / 0.05).floor() as usize).min(GC_BINS - 1);
        bins[bin] += 1.0;

        if i + window_len < seq.len() {
            gc_count -= is_gc(seq[i]);
            gc_count += is_gc(seq[i + window_len]);
        }
    }

    let sum: f32 = bins.iter().sum();
    if sum > 0.0 {
        let inv_sum = 1.0 / sum;
        for b in &mut bins {
            *b *= inv_sum;
        }
    }
    bins
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expected_gc_matches_whippet2() {
        // From Whippet2.jl/test/runtests.jl L103-106: ("AAAAAGCGCG" ^ 5) -> 50 bp, 50% GC -> bin 10 (0-based)
        let seq = b"AAAAAGCGCG".repeat(5);
        let egc = expected_gc(&seq, GC_WINDOW_LEN);
        assert!((egc[10] - 1.0).abs() < 1e-6);
        assert!((egc.iter().sum::<f32>() - 1.0).abs() < 1e-6);
    }
}
