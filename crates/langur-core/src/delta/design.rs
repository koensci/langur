//! DESeq2-style formula parser (`~ batch + condition`, `*`/`:` interactions,
//! continuous numeric covariates, explicit `var:levelA_vs_levelB` contrasts),
//! CSV/TSV auto-delimiter metadata loader with `NA` sample filtering, classic
//! `-a`/`-b` glob/list resolution, and per-event pivoted-QR collinearity handling.

use crate::delta::math::{matrix_rank_pivoted_qr, Matrix};
use crate::error::{Result, LangurError};
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Complete design specification for Beta-Binomial GLM differential splicing.
#[derive(Clone, Debug)]
pub struct DesignSpec {
    pub files: Vec<PathBuf>,
    pub sample_names: Vec<String>,
    pub x: Matrix,
    pub col_names: Vec<String>,
    /// 0-based column indices in `x` being tested (dropped in the null model).
    pub test_cols: Vec<usize>,
    /// `true` when comparing two explicit levels (A vs B) so `Psi_A`, `Psi_B`, `DeltaPsi = Psi_A - Psi_B` are reported.
    pub is_two_group: bool,
    pub group_a_mask: Vec<bool>,
    pub group_b_mask: Vec<bool>,
    pub group_a_label: String,
    pub group_b_label: String,
    pub test_label: String,
    /// `true` when testing a single continuous numeric covariate (reports average marginal slope).
    pub is_continuous_test: bool,
}

/// Build classic 2-group design matrix for `-a` / `-b` mode.
/// Column 0 = `(Intercept)` (group B baseline); Column 1 = `1.0` for Group A, `0.0` for Group B.
/// Consequently `LogOddsRatio = beta[1] > 0` iff `DeltaPsi = Psi_A - Psi_B > 0`.
pub fn build_two_group_design(
    files_a: Vec<PathBuf>,
    files_b: Vec<PathBuf>,
    label_a: &str,
    label_b: &str,
) -> Result<DesignSpec> {
    if files_a.is_empty() || files_b.is_empty() {
        return Err(LangurError::InvalidParam(
            "2-group mode requires at least 1 file in both -a and -b".into(),
        ));
    }
    let na = files_a.len();
    let nb = files_b.len();
    let ns = na + nb;
    let mut files = Vec::with_capacity(ns);
    files.extend(files_a);
    files.extend(files_b);

    let sample_names: Vec<String> = files
        .iter()
        .map(|p| {
            p.file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| p.display().to_string())
        })
        .collect();

    let mut x = Matrix::zeros(ns, 2);
    let mut mask_a = vec![false; ns];
    let mut mask_b = vec![false; ns];
    for i in 0..ns {
        x.set(i, 0, 1.0);
        if i < na {
            x.set(i, 1, 1.0);
            mask_a[i] = true;
        } else {
            mask_b[i] = true;
        }
    }

    let tlabel = format!("{label_a}_vs_{label_b}");
    Ok(DesignSpec {
        files,
        sample_names,
        x,
        col_names: vec!["(Intercept)".into(), tlabel.clone()],
        test_cols: vec![1],
        is_two_group: true,
        group_a_mask: mask_a,
        group_b_mask: mask_b,
        group_a_label: label_a.to_owned(),
        group_b_label: label_b.to_owned(),
        test_label: tlabel,
        is_continuous_test: false,
    })
}

/// Resolve a `-a` / `-b` argument string (either a comma-delimited list of file paths/basenames
/// or a substring pattern matching `*<pattern>*.psi*` inside `dir`).
pub fn resolve_glob_or_list(pattern: &str, dir: &Path) -> Result<Vec<PathBuf>> {
    let trimmed = pattern.trim();
    if trimmed.is_empty() {
        return Ok(Vec::new());
    }
    if trimmed.contains(',') {
        let mut out = Vec::new();
        for tok in trimmed.split(',') {
            let s = tok.trim();
            if s.is_empty() {
                continue;
            }
            let p = PathBuf::from(s);
            let resolved = if p.is_absolute() { p } else { dir.join(p) };
            if !resolved.is_file() {
                return Err(LangurError::InvalidParam(format!(
                    "PSI file not found: {}",
                    resolved.display()
                )));
            }
            out.push(resolved);
        }
        return Ok(out);
    }

    // Check if it's a direct single file path first
    let direct = {
        let p = PathBuf::from(trimmed);
        if p.is_absolute() { p } else { dir.join(p) }
    };
    if direct.is_file() {
        return Ok(vec![direct]);
    }

    // Otherwise scan `dir` matching Whippet2's `*<pattern>*.psi*` rule (sorted deterministically)
    let mut matches = Vec::new();
    let entries = std::fs::read_dir(dir).map_err(|source| LangurError::IoPath {
        path: dir.to_path_buf(),
        source,
    })?;
    for entry_res in entries {
        let entry = entry_res?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Some(fname) = path.file_name().and_then(|s| s.to_str()) else {
            continue;
        };
        if fname.contains(trimmed) && fname.contains(".psi") {
            matches.push(path);
        }
    }
    matches.sort();
    Ok(matches)
}

/// Parse a DESeq2-style formula like `~ batch + donor + condition` or `~ batch * condition`
/// into an ordered, deduplicated list of additive terms (each term is a list of interacting variables).
pub fn parse_formula_terms(formula: &str) -> Result<Vec<Vec<String>>> {
    let mut s = formula.trim();
    if let Some(rest) = s.strip_prefix('~') {
        s = rest.trim();
    }
    if s.is_empty() {
        return Err(LangurError::InvalidParam(format!(
            "empty model formula: '{formula}'"
        )));
    }

    let mut terms: Vec<Vec<String>> = Vec::new();
    let mut seen: FxHashSet<Vec<String>> = FxHashSet::default();

    let mut add_term = |vars: Vec<String>| {
        let mut canon = vars.clone();
        canon.sort();
        canon.dedup();
        if seen.insert(canon) {
            terms.push(vars);
        }
    };

    for raw_part in s.split('+') {
        let part = raw_part.trim();
        if part.is_empty() || part == "1" {
            continue;
        }
        if part.contains('*') {
            let factors: Vec<String> = part
                .split('*')
                .map(|x| x.trim().to_owned())
                .filter(|x| !x.is_empty())
                .collect();
            let k = factors.len();
            if k > 8 {
                return Err(LangurError::InvalidParam(format!(
                    "interaction order too high ({k}) in '{part}'"
                )));
            }
            for mask in 1usize..((1usize << k)) {
                let mut sub = Vec::new();
                for bit in 0..k {
                    if ((mask >> bit) & 1) == 1 {
                        sub.push(factors[bit].clone());
                    }
                }
                add_term(sub);
            }
        } else if part.contains(':') {
            let factors: Vec<String> = part
                .split(':')
                .map(|x| x.trim().to_owned())
                .filter(|x| !x.is_empty())
                .collect();
            if !factors.is_empty() {
                add_term(factors);
            }
        } else {
            add_term(vec![part.to_owned()]);
        }
    }

    if terms.is_empty() {
        return Err(LangurError::InvalidParam(format!(
            "formula '{formula}' contains no predictor terms"
        )));
    }
    Ok(terms)
}

#[inline]
fn is_missing_token(s: &str) -> bool {
    let t = s.trim();
    t.is_empty()
        || t.eq_ignore_ascii_case("NA")
        || t.eq_ignore_ascii_case("NaN")
        || t.eq_ignore_ascii_case("null")
        || t == "."
}

/// Split one CSV or TSV line respecting double quotes.
fn split_delimited_line(line: &str, delim: char) -> Vec<String> {
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch == '"' {
            if in_quotes && chars.peek() == Some(&'"') {
                cur.push('"');
                chars.next();
            } else {
                in_quotes = !in_quotes;
            }
        } else if ch == delim && !in_quotes {
            fields.push(cur.trim().to_owned());
            cur.clear();
        } else {
            cur.push(ch);
        }
    }
    fields.push(cur.trim().to_owned());
    fields
}

/// Auto-detect delimiter (TSV if header contains `\t`, else CSV `,`).
fn detect_delimiter(header_line: &str) -> char {
    if header_line.contains('\t') {
        '\t'
    } else {
        ','
    }
}

/// Load a CSV/TSV sample metadata table and build the full `DesignSpec` from `formula` & `test_spec`.
///
/// - Auto-detects CSV/TSV delimiter.
/// - Resolves relative PSI paths against `dir`, falling back to `design_path`'s parent directory.
/// - Drops samples with missing/`NA` values in any formula-referenced covariate (logging a warning
///   to stderr rather than silently treating `"NA"` as a categorical level).
/// - Guarantees `DeltaPsi = Psi_A - Psi_B` and `LogOddsRatio` share identical sign for both
///   explicit `var:levelA_vs_levelB` contrasts and binary categorical variables.
pub fn build_design_from_table(
    design_path: &Path,
    formula: &str,
    test_spec: &str,
    dir: &Path,
    file_col: &str,
) -> Result<DesignSpec> {
    let f = File::open(design_path).map_err(|source| LangurError::IoPath {
        path: design_path.to_path_buf(),
        source,
    })?;
    let reader = BufReader::new(f);
    let mut lines = reader.lines();

    let header_line = loop {
        match lines.next() {
            Some(Ok(l)) => {
                let t = l.trim();
                if !t.is_empty() && !t.starts_with('#') {
                    break l;
                }
            }
            Some(Err(e)) => return Err(e.into()),
            None => {
                return Err(LangurError::InvalidParam(format!(
                    "design table '{}' is empty",
                    design_path.display()
                )))
            }
        }
    };

    let delim = detect_delimiter(&header_line);
    let headers = split_delimited_line(&header_line, delim);
    let mut col_idx_map: FxHashMap<String, usize> = FxHashMap::default();
    for (i, h) in headers.iter().enumerate() {
        col_idx_map.insert(h.clone(), i);
    }

    let &file_cidx = col_idx_map.get(file_col).ok_or_else(|| {
        LangurError::InvalidParam(format!(
            "design table '{}' missing file column '{}' (available: {:?})",
            design_path.display(),
            file_col,
            headers
        ))
    })?;
    let sample_cidx = col_idx_map.get("sample").copied();

    let terms = parse_formula_terms(formula)?;

    // Parse all referenced variable names and verify they exist in the table
    let mut all_vars: Vec<String> = Vec::new();
    for term in &terms {
        for v in term {
            if !all_vars.contains(v) {
                if !col_idx_map.contains_key(v) {
                    return Err(LangurError::InvalidParam(format!(
                        "formula variable '{v}' not found in design table columns: {:?}",
                        headers
                    )));
                }
                all_vars.push(v.clone());
            }
        }
    }

    // Parse test_spec: "", "var", "varA:varB", or "var:levelA_vs_levelB"
    let ts_trimmed = test_spec.trim();
    let (test_var, contrast_num, contrast_den) = if ts_trimmed.is_empty() {
        (terms.last().unwrap().join(":"), None, None)
    } else if ts_trimmed.contains("_vs_") && ts_trimmed.contains(':') {
        let (v, rhs) = ts_trimmed.split_once(':').unwrap();
        let parts: Vec<&str> = rhs.split("_vs_").collect();
        if parts.len() != 2 {
            return Err(LangurError::InvalidParam(format!(
                "invalid contrast '{test_spec}' (expected 'column:LevelA_vs_LevelB')"
            )));
        }
        (
            v.trim().to_owned(),
            Some(parts[0].trim().to_owned()),
            Some(parts[1].trim().to_owned()),
        )
    } else {
        (ts_trimmed.to_owned(), None, None)
    };

    // Read rows, dropping any sample with NA/missing values in file_col or formula variables
    let design_parent = design_path.parent().unwrap_or_else(|| Path::new("."));
    let mut kept_rows: Vec<Vec<String>> = Vec::new();
    let mut files: Vec<PathBuf> = Vec::new();
    let mut sample_names: Vec<String> = Vec::new();
    let mut dropped_na = 0usize;

    for (row_no, line_res) in lines.enumerate() {
        let line = line_res?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let fields = split_delimited_line(&line, delim);
        if fields.len() < headers.len() {
            return Err(LangurError::InvalidParam(format!(
                "design table row {} has {} columns, expected {}",
                row_no + 2,
                fields.len(),
                headers.len()
            )));
        }

        // Check for missing values in file column or any formula variable
        let mut missing_var: Option<&str> = None;
        if is_missing_token(&fields[file_cidx]) {
            missing_var = Some(file_col);
        } else {
            for v in &all_vars {
                let ci = col_idx_map[v];
                if is_missing_token(&fields[ci]) {
                    missing_var = Some(v.as_str());
                    break;
                }
            }
        }
        if let Some(mv) = missing_var {
            dropped_na += 1;
            let sname = sample_cidx
                .map(|ci| fields[ci].as_str())
                .unwrap_or(fields[file_cidx].as_str());
            eprintln!(
                "warning: dropping sample '{sname}' (row {}) due to missing/NA value in '{mv}'",
                row_no + 2
            );
            continue;
        }

        let raw_f = &fields[file_cidx];
        let p = PathBuf::from(raw_f);
        let resolved = if p.is_absolute() {
            p
        } else {
            let cand1 = dir.join(&p);
            if cand1.is_file() {
                cand1
            } else {
                let cand2 = design_parent.join(&p);
                if cand2.is_file() {
                    cand2
                } else {
                    cand1
                }
            }
        };
        if !resolved.is_file() {
            return Err(LangurError::InvalidParam(format!(
                "PSI file not found: '{}' (from design table entry '{}')",
                resolved.display(),
                raw_f
            )));
        }

        let sname = sample_cidx
            .map(|ci| fields[ci].clone())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| {
                resolved
                    .file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_else(|| raw_f.clone())
            });

        files.push(resolved);
        sample_names.push(sname);
        kept_rows.push(fields);
    }

    if dropped_na > 0 {
        eprintln!(
            "warning: dropped {dropped_na} sample(s) with missing/NA covariates; {} sample(s) retained.",
            kept_rows.len()
        );
    }
    let ns = kept_rows.len();
    if ns < 2 {
        return Err(LangurError::InvalidParam(format!(
            "design table has fewer than 2 valid samples ({ns})"
        )));
    }

    // Encode each referenced variable into numeric (1 col) or categorical dummy (K-1 cols)
    struct EncodedVar {
        cols: Vec<Vec<f64>>,
        col_names: Vec<String>,
        is_cat: bool,
        levels: Vec<String>,
    }

    let mut var_enc: FxHashMap<String, EncodedVar> = FxHashMap::default();
    for v in &all_vars {
        let ci = col_idx_map[v];
        let raw_vals: Vec<&str> = kept_rows.iter().map(|r| r[ci].as_str()).collect();
        let force_cat = *v == test_var && contrast_den.is_some();

        let numeric_vals: Option<Vec<f64>> = if force_cat {
            None
        } else {
            raw_vals
                .iter()
                .map(|s| s.parse::<f64>().ok().filter(|x| x.is_finite()))
                .collect()
        };

        if let Some(nums) = numeric_vals {
            var_enc.insert(
                v.clone(),
                EncodedVar {
                    cols: vec![nums],
                    col_names: vec![v.clone()],
                    is_cat: false,
                    levels: Vec::new(),
                },
            );
        } else {
            let mut levels: Vec<String> = raw_vals.iter().map(|s| (*s).to_owned()).collect();
            levels.sort();
            levels.dedup();
            if levels.len() < 2 {
                return Err(LangurError::InvalidParam(format!(
                    "categorical variable '{v}' has fewer than 2 levels: {:?}",
                    levels
                )));
            }

            let ordered_levels = if *v == test_var {
                if let Some(ref ref_lv) = contrast_den {
                    if !levels.contains(ref_lv) {
                        return Err(LangurError::InvalidParam(format!(
                            "reference level '{ref_lv}' not found in '{v}' (available: {:?})",
                            levels
                        )));
                    }
                    let mut ord = Vec::with_capacity(levels.len());
                    ord.push(ref_lv.clone());
                    for lv in &levels {
                        if lv != ref_lv {
                            ord.push(lv.clone());
                        }
                    }
                    ord
                } else {
                    levels
                }
            } else {
                levels
            };

            let mut cols = Vec::with_capacity(ordered_levels.len() - 1);
            let mut cnames = Vec::with_capacity(ordered_levels.len() - 1);
            for lv in &ordered_levels[1..] {
                cnames.push(format!("{v}:{lv}"));
                let col: Vec<f64> = raw_vals
                    .iter()
                    .map(|&s| if s == lv.as_str() { 1.0 } else { 0.0 })
                    .collect();
                cols.push(col);
            }
            var_enc.insert(
                v.clone(),
                EncodedVar {
                    cols,
                    col_names: cnames,
                    is_cat: true,
                    levels: ordered_levels,
                },
            );
        }
    }

    // Assemble full design matrix X with intercept column 0
    let mut all_cols: Vec<Vec<f64>> = vec![vec![1.0; ns]];
    let mut col_names: Vec<String> = vec!["(Intercept)".into()];
    let mut term_col_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();

    for term in &terms {
        let tname = term.join(":");
        let (term_cols, term_cnames) = if term.len() == 1 {
            let ev = &var_enc[&term[0]];
            (ev.cols.clone(), ev.col_names.clone())
        } else {
            let mut cur_cols = var_enc[&term[0]].cols.clone();
            let mut cur_names = var_enc[&term[0]].col_names.clone();
            for v in &term[1..] {
                let ev2 = &var_enc[v];
                let mut next_cols = Vec::with_capacity(cur_cols.len() * ev2.cols.len());
                let mut next_names = Vec::with_capacity(cur_cols.len() * ev2.cols.len());
                for (c1, n1) in cur_cols.iter().zip(cur_names.iter()) {
                    for (c2, n2) in ev2.cols.iter().zip(ev2.col_names.iter()) {
                        let prod: Vec<f64> =
                            c1.iter().zip(c2.iter()).map(|(&a, &b)| a * b).collect();
                        next_cols.push(prod);
                        next_names.push(format!("{n1}:{n2}"));
                    }
                }
                cur_cols = next_cols;
                cur_names = next_names;
            }
            (cur_cols, cur_names)
        };

        let start_col = all_cols.len();
        let ncols = term_cols.len();
        all_cols.extend(term_cols);
        col_names.extend(term_cnames);
        term_col_map.insert(tname, (start_col..(start_col + ncols)).collect());
    }

    let x = Matrix::from_cols(ns, &all_cols);

    let mut test_cols = term_col_map.get(&test_var).cloned().ok_or_else(|| {
        LangurError::InvalidParam(format!(
            "test term '{test_var}' not found in formula terms: {:?}",
            term_col_map.keys().collect::<Vec<_>>()
        ))
    })?;

    let mut is_two_group = false;
    let mut mask_a = vec![false; ns];
    let mut mask_b = vec![false; ns];
    let mut label_a = "A".to_owned();
    let mut label_b = "B".to_owned();
    let mut is_continuous_test = false;
    let mut test_label = test_var.clone();

    if let Some(ev) = var_enc.get(&test_var) {
        if !ev.is_cat {
            is_continuous_test = true;
        } else {
            let ci = col_idx_map[&test_var];
            let raw_vals: Vec<&str> = kept_rows.iter().map(|r| r[ci].as_str()).collect();
            if let (Some(ref c_num), Some(ref c_den)) = (contrast_num, contrast_den) {
                if !ev.levels.contains(c_num) {
                    return Err(LangurError::InvalidParam(format!(
                        "contrast level '{c_num}' not found in '{test_var}' (available: {:?})",
                        ev.levels
                    )));
                }
                let target_col = format!("{test_var}:{c_num}");
                let cidx = col_names
                    .iter()
                    .position(|n| n == &target_col)
                    .ok_or_else(|| {
                        LangurError::InvalidParam(format!(
                            "could not locate contrast column '{target_col}'"
                        ))
                    })?;
                test_cols = vec![cidx];
                is_two_group = true;
                for i in 0..ns {
                    mask_a[i] = raw_vals[i] == c_num.as_str();
                    mask_b[i] = raw_vals[i] == c_den.as_str();
                }
                label_a = c_num.clone();
                label_b = c_den.clone();
                test_label = format!("{test_var}:{c_num}_vs_{c_den}");
            } else if ev.levels.len() == 2 {
                // Binary categorical: levels[0] = reference (B, 0), levels[1] = tested (A, 1)
                is_two_group = true;
                for i in 0..ns {
                    mask_a[i] = raw_vals[i] == ev.levels[1].as_str();
                    mask_b[i] = raw_vals[i] == ev.levels[0].as_str();
                }
                label_a = ev.levels[1].clone();
                label_b = ev.levels[0].clone();
                test_label = format!("{test_var}:{}_vs_{}", ev.levels[1], ev.levels[0]);
            }
        }
    }

    Ok(DesignSpec {
        files,
        sample_names,
        x,
        col_names,
        test_cols,
        is_two_group,
        group_a_mask: mask_a,
        group_b_mask: mask_b,
        group_a_label: label_a,
        group_b_label: label_b,
        test_label,
        is_continuous_test,
    })
}

/// Full-rank event sub-matrices `(x_full, x_null, p0, keep_test)` after subsetting `design.x`
/// to `valid_idx` and dropping collinear nuisance columns via column-pivoted Householder QR.
/// Always keeps `(Intercept)` (column 0) first; returns `None` if the tested contrast becomes unidentifiable.
pub fn build_event_matrices(
    x: &Matrix,
    valid_idx: &[usize],
    test_cols: &[usize],
) -> Option<(Matrix, Matrix, usize, Vec<usize>)> {
    let m = valid_idx.len();
    let p_all = x.cols;
    if m == 0 || p_all == 0 {
        return None;
    }

    let test_set: FxHashSet<usize> = test_cols.iter().copied().collect();
    let nuis_cols: Vec<usize> = (0..p_all).filter(|c| !test_set.contains(c)).collect();

    let extract_col = |c: usize| -> Vec<f64> { valid_idx.iter().map(|&r| x.get(r, c)).collect() };

    let adds_rank = |base_cols: &[Vec<f64>], cand_col: &[f64]| -> bool {
        let mut all = Vec::with_capacity(base_cols.len() + 1);
        all.extend_from_slice(base_cols);
        all.push(cand_col.to_vec());
        let mat = Matrix::from_cols(m, &all);
        matrix_rank_pivoted_qr(&mat) > base_cols.len()
    };

    // 1. Linearly independent nuisance columns (intercept first)
    let mut x0_cols: Vec<Vec<f64>> = Vec::new();
    for nc in nuis_cols {
        let col = extract_col(nc);
        if adds_rank(&x0_cols, &col) {
            x0_cols.push(col);
        }
    }
    let p0 = x0_cols.len();

    // 2. Tested columns adding independent rank beyond X0
    let mut cur_cols = x0_cols.clone();
    let mut keep_test = Vec::new();
    for &tc in test_cols {
        let col = extract_col(tc);
        if adds_rank(&cur_cols, &col) {
            keep_test.push(tc);
            cur_cols.push(col);
        }
    }

    if keep_test.is_empty() || m < cur_cols.len() {
        return None;
    }

    let x_null = Matrix::from_cols(m, &x0_cols);
    let x_full = Matrix::from_cols(m, &cur_cols);
    Some((x_full, x_null, p0, keep_test))
}
