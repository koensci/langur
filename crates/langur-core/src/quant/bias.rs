//! Unified 5'-hexamer sequence bias (Hansen et al., NAR 2010) and 20-bin
//! fragment GC-content bias model feeding into BOTH Stage 1 transcript TPMs
//! (via 2-pass VBEM) and Stage 2 local splice-graph PSI quantification.

use crate::graph::{GC_BINS, GC_WINDOW_LEN};
use crate::index::MmapIndex;
use crate::seq::ASCII_TO_2BIT;
use smallvec::SmallVec;

pub const HEXAMER_K: usize = 6;
pub const HEXAMER_SPACE: usize = 1 << (2 * HEXAMER_K); // 4096

/// Inline compact bias observation summary for a single node, edge, or equivalence class.
///
/// - During Wave 0 (`primer_weights == None`), accumulates exact `(hexamer, weight)` pairs in
///   `wave0_hexamers` so small datasets (`<= 1` wave) match full-sample primer weighting exactly,
///   then immediately collapses and drops `wave0_hexamers` at the end of Wave 0.
/// - During Waves 1..N (`primer_weights == Some(pw)`), performs O(1) scalar accumulation directly
///   into `(primer_weighted, primer_den)` plus inline `[f32; GC_BINS]` GC bins (zero `Box`,
///   zero `SmallVec` scans, and 22 float adds per shard merge).
#[derive(Clone, Debug)]
pub struct FeatureBias {
    /// Sum of `weight * primer_weight[h]` accumulated across reads supporting this feature.
    pub primer_weighted: f64,
    /// Sum of `weight` for reads with >= 1 valid 5'-hexamer supporting this feature.
    pub primer_den: f64,
    /// 20-bin observed fragment GC distribution for reads supporting this feature.
    pub gc_bins: [f32; GC_BINS],
    /// Exact `(hexamer_idx, weight)` list used ONLY during Wave 0; collapsed and set to `None`
    /// immediately after Wave 0 finishes.
    pub wave0_hexamers: Option<Box<Vec<(u16, f32)>>>,
}

impl Default for FeatureBias {
    fn default() -> Self {
        Self {
            primer_weighted: 0.0,
            primer_den: 0.0,
            gc_bins: [0.0f32; GC_BINS],
            wave0_hexamers: None,
        }
    }
}

impl FeatureBias {
    /// Record a Wave 0 observation with exact per-hexamer weights (collapsed immediately after Wave 0).
    #[inline]
    pub fn record_wave0(&mut self, hexamers: &[u16], gc_bin: Option<usize>, weight: f32) {
        if !hexamers.is_empty() {
            let list = self.wave0_hexamers.get_or_insert_with(|| Box::new(Vec::with_capacity(4)));
            let per_h = weight / (hexamers.len() as f32);
            for &h in hexamers {
                if let Some(entry) = list.iter_mut().find(|e| e.0 == h) {
                    entry.1 += per_h;
                } else {
                    list.push((h, per_h));
                }
            }
        }
        if let Some(b) = gc_bin {
            self.gc_bins[b] += weight;
        }
    }

    /// Record a Wave 1..N observation in O(1) time and zero heap allocations using locked primer weights.
    #[inline]
    pub fn record_fast(&mut self, primer_w: Option<f64>, gc_bin: Option<usize>, weight: f64) {
        if let Some(pw) = primer_w {
            self.primer_weighted += weight * pw;
            self.primer_den += weight;
        }
        if let Some(b) = gc_bin {
            self.gc_bins[b] += weight as f32;
        }
    }

    /// Collapse Wave 0 exact hexamer counts into `(primer_weighted, primer_den)` and free heap memory.
    #[inline]
    pub fn collapse_wave0(&mut self, primer_weight: &[f64]) {
        if let Some(hexes) = self.wave0_hexamers.take() {
            for (h, w) in *hexes {
                let wf = w as f64;
                self.primer_weighted += wf * primer_weight[h as usize];
                self.primer_den += wf;
            }
        }
    }

    #[inline]
    pub fn merge_from(&mut self, other: &Self) {
        self.primer_weighted += other.primer_weighted;
        self.primer_den += other.primer_den;
        for i in 0..GC_BINS {
            self.gc_bins[i] += other.gc_bins[i];
        }
        if let Some(ref oh) = other.wave0_hexamers {
            let list = self
                .wave0_hexamers
                .get_or_insert_with(|| Box::new(Vec::with_capacity(oh.len())));
            for &(h, w) in oh.iter() {
                if let Some(entry) = list.iter_mut().find(|e| e.0 == h) {
                    entry.1 += w;
                } else {
                    list.push((h, w));
                }
            }
        }
    }

    #[inline]
    pub fn merge_owned(&mut self, other: Self) {
        self.primer_weighted += other.primer_weighted;
        self.primer_den += other.primer_den;
        for i in 0..GC_BINS {
            self.gc_bins[i] += other.gc_bins[i];
        }
        if let Some(oh) = other.wave0_hexamers {
            if let Some(ref mut list) = self.wave0_hexamers {
                for (h, w) in *oh {
                    if let Some(entry) = list.iter_mut().find(|e| e.0 == h) {
                        entry.1 += w;
                    } else {
                        list.push((h, w));
                    }
                }
            } else {
                self.wave0_hexamers = Some(oh);
            }
        }
    }

    /// Compute primer-only adjusted count multiplier (relative to raw count).
    pub fn primer_factor(&self, model: &BiasModel) -> f64 {
        let mut num = self.primer_weighted;
        let mut den = self.primer_den;
        if let Some(ref hexes) = self.wave0_hexamers {
            for &(h, w) in hexes.iter() {
                let wf = w as f64;
                num += wf * model.primer_weight[h as usize];
                den += wf;
            }
        }
        if den > 0.0 {
            (num / den).clamp(0.1, 10.0)
        } else {
            1.0
        }
    }

    /// Compute GC-only adjusted count multiplier (relative to raw count).
    pub fn gc_factor(&self, model: &BiasModel) -> f64 {
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for b in 0..GC_BINS {
            let w = self.gc_bins[b] as f64;
            if w > 0.0 {
                num += w * model.gc_weight[b];
                den += w;
            }
        }
        if den > 0.0 {
            (num / den).clamp(0.1, 10.0)
        } else {
            1.0
        }
    }

    /// Combined primer * GC multiplier for unified bias correction.
    #[inline]
    pub fn joint_factor(&self, model: &BiasModel) -> f64 {
        (self.primer_factor(model) * self.gc_factor(model)).clamp(0.05, 20.0)
    }
}

/// Thread-local / global accumulator for 5'-hexamer foreground/background and observed GC bins.
#[derive(Clone, Debug)]
pub struct BiasAccum {
    pub primer_fore: Vec<u64>,
    pub primer_back: Vec<u64>,
    pub gc_fore: [f64; GC_BINS],
}

impl Default for BiasAccum {
    fn default() -> Self {
        Self {
            primer_fore: vec![0u64; HEXAMER_SPACE],
            primer_back: vec![0u64; HEXAMER_SPACE],
            gc_fore: [0.0f64; GC_BINS],
        }
    }
}

impl BiasAccum {
    /// Lightweight accumulator for Waves 1..N when `primer_weight` is already locked (zero 64 KB hexamer arrays).
    pub fn new_gc_only() -> Self {
        Self {
            primer_fore: Vec::new(),
            primer_back: Vec::new(),
            gc_fore: [0.0f64; GC_BINS],
        }
    }

    /// Extract 5'-hexamers (positions 0..=1) and GC bin (50 bp window from offset 4) from a read,
    /// updating global foreground/background accumulators (used during Wave 0).
    pub fn observe_read(&mut self, seq: &[u8]) -> (SmallVec<[u16; 2]>, Option<usize>) {
        if self.primer_fore.is_empty() {
            self.primer_fore.resize(HEXAMER_SPACE, 0);
            self.primer_back.resize(HEXAMER_SPACE, 0);
        }
        let mut fore_hex = SmallVec::<[u16; 2]>::new();
        let len = seq.len();

        // Foreground 5'-hexamers at 0-based positions 0 and 1 (matches Whippet2 foreoff=1, foren=2)
        for start in 0..2 {
            if start + HEXAMER_K <= len {
                if let Some(h) = encode_hexamer(&seq[start..start + HEXAMER_K]) {
                    self.primer_fore[h as usize] += 1;
                    fore_hex.push(h);
                }
            }
        }

        // Background interior hexamers (0-based positions 22..32, or middle of short reads)
        let back_start = if len >= 22 + HEXAMER_K {
            22
        } else if len > HEXAMER_K + 2 {
            len / 2
        } else {
            len
        };
        let back_end = (back_start + 10).min(len.saturating_sub(HEXAMER_K - 1));
        for start in back_start..back_end {
            if let Some(h) = encode_hexamer(&seq[start..start + HEXAMER_K]) {
                self.primer_back[h as usize] += 1;
            }
        }

        // Fragment GC bin over a 50 bp window (or full read if shorter than 50 bp)
        let gc_bin = compute_read_gc_bin(seq);
        if let Some(b) = gc_bin {
            self.gc_fore[b] += 1.0;
        }

        (fore_hex, gc_bin)
    }

    /// Fast Wave 1..N observation using locked Wave 0 `primer_weight[0..4096]`:
    /// looks up ONLY the 2 foreground 5'-hexamers (skips 10 background hexamers) and records GC bin.
    /// Returns `(primer_weight_sum, valid_hexamer_count, gc_bin)`.
    #[inline]
    pub fn observe_read_fast(&mut self, seq: &[u8], primer_weight: &[f64]) -> (f64, u32, Option<usize>) {
        let len = seq.len();
        let mut w_sum = 0.0f64;
        let mut n_hex = 0u32;
        for start in 0..2 {
            if start + HEXAMER_K <= len {
                if let Some(h) = encode_hexamer(&seq[start..start + HEXAMER_K]) {
                    w_sum += primer_weight[h as usize];
                    n_hex += 1;
                }
            }
        }
        let gc_bin = compute_read_gc_bin(seq);
        if let Some(b) = gc_bin {
            self.gc_fore[b] += 1.0;
        }
        (w_sum, n_hex, gc_bin)
    }

    pub fn merge_from(&mut self, other: &Self) {
        if !other.primer_fore.is_empty() {
            if self.primer_fore.is_empty() {
                self.primer_fore = other.primer_fore.clone();
                self.primer_back = other.primer_back.clone();
            } else {
                for (a, &b) in self.primer_fore.iter_mut().zip(other.primer_fore.iter()) {
                    *a += b;
                }
                for (a, &b) in self.primer_back.iter_mut().zip(other.primer_back.iter()) {
                    *a += b;
                }
            }
        }
        for i in 0..GC_BINS {
            self.gc_fore[i] += other.gc_fore[i];
        }
    }

    /// Finalize Pass 1 model (5'-hexamer weights ready; GC weights initialized to 1.0 until Pass 1 VBEM finishes).
    pub fn finalize_primer(&self) -> BiasModel {
        let fore_sum: u64 = self.primer_fore.iter().sum();
        let back_sum: u64 = self.primer_back.iter().sum();
        let mut primer_weight = vec![1.0f64; HEXAMER_SPACE];

        if fore_sum > 0 && back_sum > 0 {
            // Mass-proportional symmetric pseudocount: guarantees pb / pf == 1.0 identically
            // whenever fore[h]/fore_sum == back[h]/back_sum, regardless of 2:10 fore/back
            // window count ratio or sample depth!
            let prior_frac_per_bin = 0.25 / (HEXAMER_SPACE as f64);
            let pseudo_f = (fore_sum as f64) * prior_frac_per_bin;
            let pseudo_b = (back_sum as f64) * prior_frac_per_bin;
            let f_tot = (fore_sum as f64) + pseudo_f * (HEXAMER_SPACE as f64);
            let b_tot = (back_sum as f64) + pseudo_b * (HEXAMER_SPACE as f64);
            for i in 0..HEXAMER_SPACE {
                let pf = ((self.primer_fore[i] as f64) + pseudo_f) / f_tot;
                let pb = ((self.primer_back[i] as f64) + pseudo_b) / b_tot;
                primer_weight[i] = (pb / pf).clamp(0.1, 10.0);
            }
        }

        BiasModel {
            primer_weight,
            gc_weight: [1.0f64; GC_BINS],
            gc_fore: self.gc_fore,
            gc_back: [0.0f64; GC_BINS],
        }
    }
}

/// Fitted bias correction weights applied across Stage 1 (Pass 1 & Pass 2) and Stage 2.
#[derive(Clone, Debug)]
pub struct BiasModel {
    pub primer_weight: Vec<f64>,
    pub gc_weight: [f64; GC_BINS],
    pub gc_fore: [f64; GC_BINS],
    pub gc_back: [f64; GC_BINS],
}

impl Default for BiasModel {
    fn default() -> Self {
        Self {
            primer_weight: vec![1.0f64; HEXAMER_SPACE],
            gc_weight: [1.0f64; GC_BINS],
            gc_fore: [0.0f64; GC_BINS],
            gc_back: [0.0f64; GC_BINS],
        }
    }
}

impl BiasModel {
    /// Construct a Pass 1 `BiasModel` using pre-locked Wave 0 primer weights + full-sample `gc_fore`.
    pub fn with_locked_primer(primer_weight: Vec<f64>, gc_fore: [f64; GC_BINS]) -> Self {
        Self {
            primer_weight,
            gc_weight: [1.0f64; GC_BINS],
            gc_fore,
            gc_back: [0.0f64; GC_BINS],
        }
    }

    /// Fit 20-bin GC bias weights against the expressed transcriptome's `tx_expected_gc`
    /// using Pass 1 transcript counts (`tx_counts`), ready for warm-started Stage 1 Pass 2
    /// and Stage 2 event PSI.
    pub fn fit_gc_from_tx_counts(&mut self, idx: &MmapIndex, tx_counts: &[f64]) {
        self.gc_back = [0.0f64; GC_BINS];
        for (t_idx, &cnt) in tx_counts.iter().enumerate() {
            if cnt <= 0.0 {
                continue;
            }
            let egc = idx.transcript_gc(t_idx);
            for b in 0..GC_BINS.min(egc.len()) {
                self.gc_back[b] += (egc[b] as f64) * cnt;
            }
        }

        let fore_sum: f64 = self.gc_fore.iter().sum();
        let back_sum: f64 = self.gc_back.iter().sum();
        if fore_sum <= 0.0 || back_sum <= 0.0 {
            self.gc_weight = [1.0f64; GC_BINS];
            return;
        }

        // Smooth pseudocount (0.5% mass per bin) prevents extreme weights in empty GC extremes (0% / 100% GC)
        let pseudo_f = fore_sum * 0.005;
        let pseudo_b = back_sum * 0.005;
        let f_tot = fore_sum + pseudo_f * (GC_BINS as f64);
        let b_tot = back_sum + pseudo_b * (GC_BINS as f64);

        for b in 0..GC_BINS {
            let pf = (self.gc_fore[b] + pseudo_f) / f_tot;
            let pb = (self.gc_back[b] + pseudo_b) / b_tot;
            self.gc_weight[b] = (pb / pf).clamp(0.15, 6.0);
        }
    }

    /// Compute the transcript-level GC effective-length bias factor:
    /// `sum_b tx_expected_gc[t][b] / gc_weight[b]` (expected relative capture rate of transcript `t`,
    /// where `gc_weight[b] = back[b] / fore[b]` so `1 / gc_weight[b] = fore[b] / back[b]`).
    /// Multiplying `eff_len[t]` by this factor folds GC bias directly into Stage 1 Pass 2 VBEM & TPMs!
    pub fn transcript_gc_capture_rate(&self, idx: &MmapIndex, tx_idx: usize) -> f64 {
        let egc = idx.transcript_gc(tx_idx);
        let mut rate = 0.0f64;
        let mut mass = 0.0f64;
        for b in 0..GC_BINS.min(egc.len()) {
            let e = egc[b] as f64;
            if e > 0.0 {
                rate += e / self.gc_weight[b];
                mass += e;
            }
        }
        if mass > 0.0 {
            (rate / mass).clamp(0.15, 6.0)
        } else {
            1.0
        }
    }
}

#[inline]
fn encode_hexamer(slice: &[u8]) -> Option<u16> {
    let mut acc: u16 = 0;
    for &b in slice {
        let code = ASCII_TO_2BIT[b as usize];
        if code == 0xFF {
            return None;
        }
        acc = (acc << 2) | (code as u16);
    }
    Some(acc)
}

#[inline]
pub fn compute_read_gc_bin(seq: &[u8]) -> Option<usize> {
    if seq.is_empty() {
        return None;
    }
    let start = if seq.len() >= GC_WINDOW_LEN + 4 { 4 } else { 0 };
    let win_len = GC_WINDOW_LEN.min(seq.len() - start);
    if win_len == 0 {
        return None;
    }
    let mut gc = 0usize;
    for &b in &seq[start..start + win_len] {
        if matches!(b, b'G' | b'g' | b'C' | b'c') {
            gc += 1;
        }
    }
    let frac = (gc as f64) / (win_len as f64);
    Some(((frac / 0.05).floor() as usize).min(GC_BINS - 1))
}
