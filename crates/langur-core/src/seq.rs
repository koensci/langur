//! DNA sequence encoding, reverse complement, 2-bit packing, bijective
//! zero-collision hashing (`k <= 27`), and O(n) monotonic-deque minimizer extraction.

use std::collections::VecDeque;

/// Lookup table mapping ASCII byte -> 2-bit DNA value (`A=0, C=1, G=2, T/U=3`, else `0xFF`).
pub const ASCII_TO_2BIT: [u8; 256] = {
    let mut lut = [0xFFu8; 256];
    lut[b'A' as usize] = 0;
    lut[b'a' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'c' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'g' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b't' as usize] = 3;
    lut[b'U' as usize] = 3;
    lut[b'u' as usize] = 3;
    lut
};

/// Lookup table mapping 2-bit value (0..=3) -> uppercase ASCII base.
pub const TWOBIT_TO_ASCII: [u8; 4] = *b"ACGT";

/// Lookup table mapping ASCII byte -> uppercase normalized ASCII (`A/C/G/T`, or `N` for ambiguous/other).
pub const ASCII_NORM: [u8; 256] = {
    let mut lut = [b'N'; 256];
    lut[b'A' as usize] = b'A';
    lut[b'a' as usize] = b'A';
    lut[b'C' as usize] = b'C';
    lut[b'c' as usize] = b'C';
    lut[b'G' as usize] = b'G';
    lut[b'g' as usize] = b'G';
    lut[b'T' as usize] = b'T';
    lut[b't' as usize] = b'T';
    lut[b'U' as usize] = b'T';
    lut[b'u' as usize] = b'T';
    lut
};

/// Lookup table mapping ASCII byte -> uppercase complement (IUPAC-aware, `N` fallback).
pub const ASCII_COMPLEMENT: [u8; 256] = {
    let mut lut = [b'N'; 256];
    lut[b'A' as usize] = b'T';
    lut[b'a' as usize] = b'T';
    lut[b'C' as usize] = b'G';
    lut[b'c' as usize] = b'G';
    lut[b'G' as usize] = b'C';
    lut[b'g' as usize] = b'C';
    lut[b'T' as usize] = b'A';
    lut[b't' as usize] = b'A';
    lut[b'U' as usize] = b'A';
    lut[b'u' as usize] = b'A';
    lut[b'R' as usize] = b'Y';
    lut[b'r' as usize] = b'Y';
    lut[b'Y' as usize] = b'R';
    lut[b'y' as usize] = b'R';
    lut[b'S' as usize] = b'S';
    lut[b's' as usize] = b'S';
    lut[b'W' as usize] = b'W';
    lut[b'w' as usize] = b'W';
    lut[b'K' as usize] = b'M';
    lut[b'k' as usize] = b'M';
    lut[b'M' as usize] = b'K';
    lut[b'm' as usize] = b'K';
    lut[b'B' as usize] = b'V';
    lut[b'b' as usize] = b'V';
    lut[b'V' as usize] = b'B';
    lut[b'v' as usize] = b'B';
    lut[b'D' as usize] = b'H';
    lut[b'd' as usize] = b'H';
    lut[b'H' as usize] = b'D';
    lut[b'h' as usize] = b'D';
    lut
};

#[inline]
pub fn base_to_2bit(b: u8) -> u8 {
    ASCII_TO_2BIT[b as usize]
}

#[inline]
pub fn is_ambiguous(b: u8) -> bool {
    ASCII_TO_2BIT[b as usize] == 0xFF
}

/// Append uppercase-normalized sequence (`A/C/G/T/N`) to `dst`. Returns `true` if any ambiguous base was seen.
pub fn append_normalized(seq: &[u8], dst: &mut Vec<u8>) -> bool {
    let mut has_ambig = false;
    dst.reserve(seq.len());
    for &b in seq {
        let nb = ASCII_NORM[b as usize];
        if nb == b'N' {
            has_ambig = true;
        }
        dst.push(nb);
    }
    has_ambig
}

/// Append uppercase reverse-complemented (`A/C/G/T/N`) sequence to `dst`. Returns `true` if any ambiguous base was seen.
pub fn append_revcomp_normalized(seq: &[u8], dst: &mut Vec<u8>) -> bool {
    let mut has_ambig = false;
    dst.reserve(seq.len());
    for &b in seq.iter().rev() {
        let cb = ASCII_NORM[ASCII_COMPLEMENT[b as usize] as usize];
        if cb == b'N' {
            has_ambig = true;
        }
        dst.push(cb);
    }
    has_ambig
}

/// Return uppercase reverse complement of `seq`.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    append_revcomp_normalized(seq, &mut out);
    out
}

/// Pack uppercase/mixed ASCII DNA into 2-bit bytes (4 bases per byte, LSB-first within each byte).
/// Ambiguous bases (`0xFF`) are safely mapped to `0` (`A`), never `3` (`T`).
/// Always pads the final byte to a byte boundary.
#[inline]
fn ascii_to_2bit_safe(b: u8) -> u8 {
    let v = ASCII_TO_2BIT[b as usize];
    if v == 0xFF { 0 } else { v }
}

pub fn pack_2bit_append(seq: &[u8], dst: &mut Vec<u8>) -> usize {
    let start_len = dst.len();
    let n_bytes = seq.len().div_ceil(4);
    dst.reserve(n_bytes);
    let (chunks, rem) = seq.as_chunks::<4>();
    for ch in chunks {
        let b0 = ascii_to_2bit_safe(ch[0]);
        let b1 = ascii_to_2bit_safe(ch[1]);
        let b2 = ascii_to_2bit_safe(ch[2]);
        let b3 = ascii_to_2bit_safe(ch[3]);
        dst.push(b0 | (b1 << 2) | (b2 << 4) | (b3 << 6));
    }
    if !rem.is_empty() {
        let mut byte = 0u8;
        for (i, &b) in rem.iter().enumerate() {
            byte |= ascii_to_2bit_safe(b) << (2 * i);
        }
        dst.push(byte);
    }
    dst.len() - start_len
}

/// Unpack `len` bases from a byte-aligned 2-bit packed buffer into ASCII `A/C/G/T`.
pub fn unpack_2bit(packed: &[u8], base_start: usize, len: usize, dst: &mut Vec<u8>) {
    dst.clear();
    dst.reserve(len);
    for pos in base_start..(base_start + len) {
        let byte = packed[pos >> 2];
        let code = (byte >> ((pos & 3) << 1)) & 3;
        dst.push(TWOBIT_TO_ASCII[code as usize]);
    }
}

/// Returns `true` iff `2 * k <= hash_bits`, in which case `hash_kmer` is a strict mathematical
/// bijection on `[0, 2^(2k))` with zero possible collisions.
#[inline]
pub const fn is_bijective_hash(k: usize, hash_bits: u32) -> bool {
    (2 * k) <= (hash_bits as usize)
}

/// Invertible bijective mixer on `bits`-bit integers (`bits` in `34..=64`).
/// Every step is a bijection mod `2^bits` (xorshift right and multiplication by an odd constant).
#[inline]
pub fn mix_bits(mut x: u64, bits: u32) -> u64 {
    debug_assert!((32..=64).contains(&bits));
    let mask = if bits == 64 {
        u64::MAX
    } else {
        (1u64 << bits) - 1
    };
    let s1 = bits / 2;
    let s2 = (bits / 2).saturating_sub(4).max(1);
    x &= mask;
    x ^= x >> s1;
    x = x.wrapping_mul(0x3C79_AC49_2BA7_B653) & mask;
    x ^= x >> s2;
    x = x.wrapping_mul(0x1C69_B3F7_4AC4_AE35) & mask;
    x ^= x >> s1;
    x
}

/// Full 64-bit bijective SplitMix64 finalizer (used to mix all `2k` bits when `2k > hash_bits`).
#[inline]
pub fn mix64(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^= x >> 31;
    x
}

/// Hash a canonical 2-bit k-mer into `hash_bits` bits (`bucket_bits + 32`, default `22 + 32 = 54`).
///
/// - When `2 * k <= hash_bits` (all `k <= 27` at default 54 bits), uses `mix_bits` directly,
///   which is a strict mathematical bijection with **zero collisions**.
/// - When `2 * k > hash_bits` (`k > 27`), mixes all 64 bits via `mix64` and truncates to `hash_bits`.
#[inline]
pub fn hash_kmer(canon_kmer: u64, k: usize, hash_bits: u32) -> u64 {
    if is_bijective_hash(k, hash_bits) {
        mix_bits(canon_kmer, hash_bits)
    } else {
        let mask = if hash_bits == 64 {
            u64::MAX
        } else {
            (1u64 << hash_bits) - 1
        };
        mix64(canon_kmer) & mask
    }
}

/// A single canonical minimizer extracted from a sequence.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Minimizer {
    /// `hash_bits`-bit hash (`bucket = hash >> 32`, `key_lo = hash as u32`).
    pub hash: u64,
    /// Raw canonical 2-bit k-mer (for collision detection when `k > 27`).
    pub canon_kmer: u64,
    /// 0-based start offset of the k-mer within the input sequence.
    pub pos: u32,
    /// `true` iff the canonical k-mer came from the reverse complement strand.
    pub is_rc: bool,
}

#[derive(Copy, Clone)]
struct DequeElem {
    pos: u32,
    hash: u64,
    canon: u64,
    is_rc: bool,
}

/// Extract canonical odd-`k` minimizers in `O(seq.len())` time using a monotonic deque.
///
/// - Ambiguous bases reset the rolling k-mer and window cleanly (never emit k-mers with `N`).
/// - Ties in hash value select the rightmost k-mer in the window (`minimap2` convention).
/// - Short valid runs of length `k <= run < k + w - 1` emit their single window minimum at run end,
///   guaranteeing that every valid exon of length `>= k` receives at least one minimizer.
pub fn extract_minimizers(
    seq: &[u8],
    k: usize,
    w: usize,
    hash_bits: u32,
    out: &mut Vec<Minimizer>,
) {
    debug_assert!((1..=31).contains(&k) && k % 2 == 1, "k must be odd in 1..=31");
    debug_assert!(w >= 1);
    out.clear();
    if seq.len() < k {
        return;
    }

    if w < 32 {
        extract_minimizers_stack(seq, k, w, hash_bits, out);
    } else {
        extract_minimizers_heap(seq, k, w, hash_bits, out);
    }
}

/// Zero-allocation stack ring-buffer monotonic deque for `w < 32` (eliminates 80M+ heap allocs/sample).
fn extract_minimizers_stack(
    seq: &[u8],
    k: usize,
    w: usize,
    hash_bits: u32,
    out: &mut Vec<Minimizer>,
) {
    let kmer_mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let rc_shift = 2 * (k - 1);

    let mut fwd: u64 = 0;
    let mut rc: u64 = 0;
    let mut valid_bases: usize = 0;
    let mut valid_kmers_in_run: usize = 0;
    let mut emitted_in_run = false;
    let mut last_emitted_pos: u32 = u32::MAX;

    const CAP: usize = 32;
    const MASK: usize = CAP - 1;
    let mut ring = [DequeElem {
        pos: 0,
        hash: 0,
        canon: 0,
        is_rc: false,
    }; CAP];
    let mut head: usize = 0;
    let mut tail: usize = 0;

    for (i, &base) in seq.iter().enumerate() {
        let code = ASCII_TO_2BIT[base as usize];
        if code == 0xFF {
            if !emitted_in_run && valid_kmers_in_run > 0 && head != tail {
                let front = ring[head & MASK];
                out.push(Minimizer {
                    hash: front.hash,
                    canon_kmer: front.canon,
                    pos: front.pos,
                    is_rc: front.is_rc,
                });
            }
            fwd = 0;
            rc = 0;
            valid_bases = 0;
            valid_kmers_in_run = 0;
            emitted_in_run = false;
            last_emitted_pos = u32::MAX;
            head = 0;
            tail = 0;
            continue;
        }

        fwd = ((fwd << 2) | (code as u64)) & kmer_mask;
        rc = (rc >> 2) | (((3 - code) as u64) << rc_shift);
        valid_bases += 1;

        if valid_bases >= k {
            let pos = (i + 1 - k) as u32;
            let is_rc = rc < fwd;
            let canon = if is_rc { rc } else { fwd };
            let hash = hash_kmer(canon, k, hash_bits);
            valid_kmers_in_run += 1;

            while tail != head && ring[(tail - 1) & MASK].hash >= hash {
                tail -= 1;
            }
            ring[tail & MASK] = DequeElem {
                pos,
                hash,
                canon,
                is_rc,
            };
            tail += 1;

            let min_pos = pos.saturating_sub((w - 1) as u32);
            while head != tail && ring[head & MASK].pos < min_pos {
                head += 1;
            }

            if valid_kmers_in_run >= w && head != tail {
                let front = ring[head & MASK];
                if front.pos != last_emitted_pos {
                    out.push(Minimizer {
                        hash: front.hash,
                        canon_kmer: front.canon,
                        pos: front.pos,
                        is_rc: front.is_rc,
                    });
                    last_emitted_pos = front.pos;
                    emitted_in_run = true;
                }
            }
        }
    }

    if !emitted_in_run && valid_kmers_in_run > 0 && head != tail {
        let front = ring[head & MASK];
        out.push(Minimizer {
            hash: front.hash,
            canon_kmer: front.canon,
            pos: front.pos,
            is_rc: front.is_rc,
        });
    }
}

fn extract_minimizers_heap(
    seq: &[u8],
    k: usize,
    w: usize,
    hash_bits: u32,
    out: &mut Vec<Minimizer>,
) {

    let kmer_mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let rc_shift = 2 * (k - 1);

    let mut fwd: u64 = 0;
    let mut rc: u64 = 0;
    let mut valid_bases: usize = 0;
    let mut valid_kmers_in_run: usize = 0;
    let mut emitted_in_run = false;
    let mut last_emitted_pos: u32 = u32::MAX;
    let mut deque: VecDeque<DequeElem> = VecDeque::with_capacity(w + 1);

    for (i, &base) in seq.iter().enumerate() {
        let code = ASCII_TO_2BIT[base as usize];
        if code == 0xFF {
            // Flush short valid run before resetting
            if !emitted_in_run && valid_kmers_in_run > 0 {
                if let Some(front) = deque.front() {
                    out.push(Minimizer {
                        hash: front.hash,
                        canon_kmer: front.canon,
                        pos: front.pos,
                        is_rc: front.is_rc,
                    });
                }
            }
            fwd = 0;
            rc = 0;
            valid_bases = 0;
            valid_kmers_in_run = 0;
            emitted_in_run = false;
            last_emitted_pos = u32::MAX;
            deque.clear();
            continue;
        }

        fwd = ((fwd << 2) | (code as u64)) & kmer_mask;
        rc = (rc >> 2) | (((3 - code) as u64) << rc_shift);
        valid_bases += 1;

        if valid_bases >= k {
            let pos = (i + 1 - k) as u32;
            let is_rc = rc < fwd;
            let canon = if is_rc { rc } else { fwd };
            let hash = hash_kmer(canon, k, hash_bits);
            valid_kmers_in_run += 1;

            // Pop elements from back that have hash >= current (ties pick rightmost)
            while let Some(back) = deque.back() {
                if back.hash >= hash {
                    deque.pop_back();
                } else {
                    break;
                }
            }
            deque.push_back(DequeElem {
                pos,
                hash,
                canon,
                is_rc,
            });

            // Evict elements outside current window [pos - w + 1, pos]
            let min_pos = pos.saturating_sub((w - 1) as u32);
            while let Some(front) = deque.front() {
                if front.pos < min_pos {
                    deque.pop_front();
                } else {
                    break;
                }
            }

            if valid_kmers_in_run >= w {
                if let Some(front) = deque.front() {
                    if front.pos != last_emitted_pos {
                        out.push(Minimizer {
                            hash: front.hash,
                            canon_kmer: front.canon,
                            pos: front.pos,
                            is_rc: front.is_rc,
                        });
                        last_emitted_pos = front.pos;
                        emitted_in_run = true;
                    }
                }
            }
        }
    }

    // Flush short valid run at end of sequence
    if !emitted_in_run && valid_kmers_in_run > 0 {
        if let Some(front) = deque.front() {
            out.push(Minimizer {
                hash: front.hash,
                canon_kmer: front.canon,
                pos: front.pos,
                is_rc: front.is_rc,
            });
        }
    }
}

/// Encode up to 16 bases MSB-aligned in a `u32` (for `acceptor_flank_2bit`).
/// Stops at the first ambiguous base (`N`/IUPAC) and returns `(packed_u32, valid_len)`.
pub fn encode_flank_msb(seq: &[u8]) -> (u32, u8) {
    let limit = seq.len().min(16);
    let mut acc: u32 = 0;
    let mut len: u8 = 0;
    for &b in &seq[..limit] {
        let code = ASCII_TO_2BIT[b as usize];
        if code == 0xFF {
            break;
        }
        let shift = 30 - 2 * (len as u32);
        acc |= (code as u32) << shift;
        len += 1;
    }
    (acc, len)
}

/// Encode up to 16 bases ending at the 3' end of `seq`, LSB-aligned in a `u32` (for `donor_flank_2bit`).
/// Walks backward from the end and stops at the first ambiguous base; returns `(packed_u32, valid_len)`.
pub fn encode_flank_lsb_from_end(seq: &[u8]) -> (u32, u8) {
    let limit = seq.len().min(16);
    let mut acc: u32 = 0;
    let mut len: u8 = 0;
    for &b in seq.iter().rev().take(limit) {
        let code = ASCII_TO_2BIT[b as usize];
        if code == 0xFF {
            break;
        }
        acc |= (code as u32) << (2 * (len as u32));
        len += 1;
    }
    (acc, len)
}

/// Count 2-bit base mismatches between the first `cmp_len` (`1..=16`) bases of two MSB-aligned `u32` flanks.
#[inline]
pub fn mismatch_count_msb(a: u32, b: u32, cmp_len: u8) -> u32 {
    debug_assert!((1..=16).contains(&cmp_len));
    let shift = 32 - 2 * (cmp_len as u32);
    let diff = (a ^ b) >> shift;
    ((diff | (diff >> 1)) & 0x5555_5555).count_ones()
}

/// Count 2-bit base mismatches between the last `cmp_len` (`1..=16`) bases of two LSB-aligned `u32` flanks.
#[inline]
pub fn mismatch_count_lsb(a: u32, b: u32, cmp_len: u8) -> u32 {
    debug_assert!((1..=16).contains(&cmp_len));
    let mask = if cmp_len == 16 {
        u32::MAX
    } else {
        (1u32 << (2 * (cmp_len as u32))) - 1
    };
    let diff = (a ^ b) & mask;
    ((diff | (diff >> 1)) & 0x5555_5555).count_ones()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_revcomp_and_2bit_roundtrip() {
        let s = b"ACGTacgtNNTGC";
        assert_eq!(revcomp(s), b"GCANNACGTACGT");

        let clean = b"ACGTTGCAA";
        let mut packed = Vec::new();
        let nbytes = pack_2bit_append(clean, &mut packed);
        assert_eq!(nbytes, 3);
        let mut unpacked = Vec::new();
        unpack_2bit(&packed, 0, clean.len(), &mut unpacked);
        assert_eq!(unpacked, clean);

        // Verify ambiguous 'N' maps safely to 0 ('A'), NOT 3 ('T')
        let mut packed_n = Vec::new();
        pack_2bit_append(b"NNNN", &mut packed_n);
        assert_eq!(packed_n, vec![0x00]);
        unpack_2bit(&packed_n, 0, 4, &mut unpacked);
        assert_eq!(unpacked, b"AAAA");
    }

    #[test]
    fn test_flank_mismatch_masks_padding() {
        // 3 bp exon "CCT" vs "CCG" (1 mismatch) and vs "CCT" (0 mismatches),
        // even though zero-padded bits would otherwise decode as "AAAA..."
        let (a, la) = encode_flank_msb(b"CCT");
        let (b, lb) = encode_flank_msb(b"CCTTTTT");
        assert_eq!(la, 3);
        assert_eq!(lb, 7);
        assert_eq!(mismatch_count_msb(a, b, la), 0);

        let (c, _) = encode_flank_msb(b"CCGAAAA");
        assert_eq!(mismatch_count_msb(a, c, la), 1);

        let (d, ld) = encode_flank_lsb_from_end(b"TTC");
        let (e, _) = encode_flank_lsb_from_end(b"GGGGTTC");
        assert_eq!(ld, 3);
        assert_eq!(mismatch_count_lsb(d, e, ld), 0);
        let (f, _) = encode_flank_lsb_from_end(b"GGGGTTA");
        assert_eq!(mismatch_count_lsb(d, f, ld), 1);
    }

    #[test]
    fn test_minimizer_rc_invariance_and_short_runs() {
        let seq = b"GCGGATTACATTTTTTTTTTGCATTAGAAGGGGGGGGGG";
        let rc_seq = revcomp(seq);
        let mut m1 = Vec::new();
        let mut m2 = Vec::new();
        extract_minimizers(seq, 7, 5, 54, &mut m1);
        extract_minimizers(&rc_seq, 7, 5, 54, &mut m2);
        assert!(!m1.is_empty());
        assert_eq!(m1.len(), m2.len());
        // Canonical k-mer hashes must match in reverse order
        let mut h1: Vec<u64> = m1.iter().map(|m| m.hash).collect();
        let mut h2: Vec<u64> = m2.iter().map(|m| m.hash).collect();
        h1.sort_unstable();
        h2.sort_unstable();
        assert_eq!(h1, h2);

        // Short sequence k <= len < k + w - 1 still emits 1 minimizer
        let short = b"GCGGATTACA"; // len 10, k=7, w=8 -> 4 kmers < w
        let mut ms = Vec::new();
        extract_minimizers(short, 7, 8, 54, &mut ms);
        assert_eq!(ms.len(), 1);
    }
}
