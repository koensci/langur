//! Exact Whippet2 Contiguous Splice Graph `EdgeType` and `EdgeMotif` semantics
//! (Sterne-Weiler et al., Molecular Cell 2018).

/// Boundary edge type between adjacent CSG nodes (3-bit encoding matching Whippet2).
#[repr(u8)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum EdgeType {
    /// `'SL'` = `0b000`: Transcript Start
    SL = 0,
    /// `'SR'` = `0b001`: Tandem Last Exon boundary (acceptor after prior txEnd)
    SR = 1,
    /// `'LS'` = `0b010`: Tandem First Exon boundary (donor before next txStart)
    LS = 2,
    /// `'RS'` = `0b011`: Transcript End
    RS = 3,
    /// `'LR'` = `0b100`: Intron 5' -> 3' splice-site hard boundary
    LR = 4,
    /// `'LL'` = `0b101`: Alternative 5' splice site (donor inside exon)
    LL = 5,
    /// `'RR'` = `0b110`: Alternative 3' splice site (acceptor inside exon)
    RR = 6,
}

impl EdgeType {
    #[inline]
    pub const fn from_u8(v: u8) -> Option<Self> {
        match v {
            0 => Some(Self::SL),
            1 => Some(Self::SR),
            2 => Some(Self::LS),
            3 => Some(Self::RS),
            4 => Some(Self::LR),
            5 => Some(Self::LL),
            6 => Some(Self::RR),
            _ => None,
        }
    }

    #[inline]
    pub const fn as_u8(self) -> u8 {
        self as u8
    }

    /// Invert edge type for `-` strand orientation (`LR` is self-inverse; others XOR with `0b011`).
    #[inline]
    pub const fn invert(self) -> Self {
        match self {
            Self::LR => Self::LR,
            Self::SL => Self::RS,
            Self::SR => Self::LS,
            Self::LS => Self::SR,
            Self::RS => Self::SL,
            Self::LL => Self::RR,
            Self::RR => Self::LL,
        }
    }

    #[inline]
    pub const fn is_tx_start(self) -> bool {
        matches!(self, Self::SL | Self::LS)
    }

    #[inline]
    pub const fn is_tx_stop(self) -> bool {
        matches!(self, Self::SR | Self::RS)
    }

    #[inline]
    pub const fn is_exon_start(self) -> bool {
        matches!(self, Self::LR | Self::RR | Self::SR | Self::LS | Self::SL)
    }

    #[inline]
    pub const fn is_exon_stop(self) -> bool {
        matches!(self, Self::LR | Self::LL | Self::SR | Self::LS | Self::RS)
    }

    #[inline]
    pub const fn is_hard_edge(self) -> bool {
        matches!(self, Self::LR | Self::LS | Self::SR)
    }

    #[inline]
    pub const fn is_acceptor(self) -> bool {
        matches!(self, Self::LR | Self::RR | Self::SR)
    }

    #[inline]
    pub const fn is_donor(self) -> bool {
        matches!(self, Self::LR | Self::LL | Self::LS)
    }

    /// Whippet2's `is_edge(edge, left)`: `left == true` checks donor, `left == false` checks acceptor.
    #[inline]
    pub const fn is_splice_edge(self, left_donor: bool) -> bool {
        if left_donor {
            self.is_donor()
        } else {
            self.is_acceptor()
        }
    }

    /// Determine edge type from coordinate boundary indices (`0=txen, 1=acc, 2=don, 3=txst`),
    /// matching `Whippet2.jl/src/graph.jl::get_edgetype`.
    #[inline]
    pub const fn from_coord_indices(min_idx: usize, sec_idx: usize, is_node: bool, strand: bool) -> Self {
        // Transposed tables from Whippet2.jl/src/graph.jl lines 14-22:
        const INDEX_TO_EDGETYPE: [[EdgeType; 4]; 4] = [
            [EdgeType::RS, EdgeType::SR, EdgeType::RS, EdgeType::RS], // min=txen
            [EdgeType::RR, EdgeType::RR, EdgeType::RR, EdgeType::RR], // min=acc
            [EdgeType::LL, EdgeType::LR, EdgeType::LL, EdgeType::LS], // min=don
            [EdgeType::SL, EdgeType::SL, EdgeType::SL, EdgeType::SL], // min=txst
        ];
        const INDEX_TO_EDGETYPE_NODE: [[EdgeType; 4]; 4] = [
            [EdgeType::RS; 4], // min=txen
            [EdgeType::RR; 4], // min=acc
            [EdgeType::LL; 4], // min=don
            [EdgeType::SL; 4], // min=txst
        ];

        let raw = if is_node {
            INDEX_TO_EDGETYPE_NODE[min_idx][sec_idx]
        } else {
            INDEX_TO_EDGETYPE[min_idx][sec_idx]
        };
        if strand {
            raw
        } else {
            raw.invert()
        }
    }
}

/// Whippet2 event taxonomy derived from `(left_edge, right_edge)` of a CSG node.
#[repr(u8)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum EdgeMotif {
    TS = 0, // Tandem Transcription Start Site
    TE = 1, // Tandem Alternative Polyadenylation (Tx End)
    AF = 2, // Alternative First Exon
    AL = 3, // Alternative Last Exon
    RI = 4, // Retained Intron
    CE = 5, // Core / Cassette Exon
    AD = 6, // Alternative Donor (5'SS)
    AA = 7, // Alternative Acceptor (3'SS)
    NA = 8, // None / Constitutive
}

impl EdgeMotif {
    #[inline]
    pub const fn from_edges(left: EdgeType, right: EdgeType) -> Self {
        use EdgeType::*;
        match (left, right) {
            (SL, SL) | (LS, SL) => Self::TS,
            (RS, RS) | (RS, SR) => Self::TE,
            (LS, LL) | (LS, LR) | (LS, LS) => Self::AF,
            (RR, SR) | (LR, SR) | (SR, SR) => Self::AL,
            (LL, RR) => Self::RI,
            (RR, LL)
            | (RR, LR)
            | (RR, LS)
            | (LR, LL)
            | (LR, LR)
            | (LR, LS)
            | (SR, LL)
            | (SR, LR)
            | (SR, LS) => Self::CE,
            (LL, LL) | (LL, LR) | (LL, LS) => Self::AD,
            (RR, RR) | (LR, RR) | (SR, RR) => Self::AA,
            _ => Self::NA,
        }
    }

    #[inline]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::TS => "TS",
            Self::TE => "TE",
            Self::AF => "AF",
            Self::AL => "AL",
            Self::RI => "RI",
            Self::CE => "CE",
            Self::AD => "AD",
            Self::AA => "AA",
            Self::NA => "NA",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edge_inversion_and_motifs() {
        assert_eq!(EdgeType::LR.invert(), EdgeType::LR);
        assert_eq!(EdgeType::SL.invert(), EdgeType::RS);
        assert_eq!(EdgeType::RS.invert(), EdgeType::SL);
        assert_eq!(EdgeType::LS.invert(), EdgeType::SR);
        assert_eq!(EdgeType::SR.invert(), EdgeType::LS);
        assert_eq!(EdgeType::LL.invert(), EdgeType::RR);
        assert_eq!(EdgeType::RR.invert(), EdgeType::LL);

        assert_eq!(EdgeMotif::from_edges(EdgeType::LL, EdgeType::RR), EdgeMotif::RI);
        assert_eq!(EdgeMotif::from_edges(EdgeType::LR, EdgeType::LR), EdgeMotif::CE);
        assert_eq!(EdgeMotif::from_edges(EdgeType::LR, EdgeType::RR), EdgeMotif::AA);
        assert_eq!(EdgeMotif::from_edges(EdgeType::LL, EdgeType::LR), EdgeMotif::AD);
        assert_eq!(EdgeMotif::from_edges(EdgeType::SL, EdgeType::SL), EdgeMotif::TS);
        assert_eq!(EdgeMotif::from_edges(EdgeType::RS, EdgeType::RS), EdgeMotif::TE);
        assert_eq!(EdgeMotif::from_edges(EdgeType::LS, EdgeType::LR), EdgeMotif::AF);
        assert_eq!(EdgeMotif::from_edges(EdgeType::LR, EdgeType::SR), EdgeMotif::AL);
    }
}
