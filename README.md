[![DOI (paper)](https://zenodo.org/badge/DOI/10.5281/zenodo.19601866.svg)](https://doi.org/10.5281/zenodo.19601866)
[![DOI (code)](https://zenodo.org/badge/DOI/10.5281/zenodo.19601792.svg)](https://doi.org/10.5281/zenodo.19601792)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Verification Scripts for "Nonexistence of 3×3 Magic Squares of Squares"

This repository contains the SymPy-based verification scripts accompanying the paper:

**Hiroki Fukui**, *Nonexistence of 3×3 Magic Squares of Squares* (2026).

The paper proves that no 3×3 magic square with all nine entries distinct perfect squares exists. The proof is entirely algebraic; these scripts provide independent computational verification of specific lemmas, propositions, and theorems.

## Requirements

- Python 3.9+
- SymPy (install: `pip install sympy`)

All scripts run in under one minute on standard hardware (total runtime ≈ 15 seconds).

## Scripts

The verification suite consists of **17 scripts** organised into three groups, mirroring the structure of the proof.

### Group A: Pillar A — §11 (case k = 1, split primes)

Located in `scripts/phase_A_pillar_A_verification/`.

| Script | Verifies | Paper Reference |
|--------|----------|-----------------|
| `verify_P2_identity.py` | Algebraic identity Δ_l = 2 p^(2l) · J_(e−l) for Step P2 (Stage 1) | §11 (Pillar A) |
| `verify_P2_nu_p.py` | Structural proof that ν_p(J_k) = 0 for k ≥ 1 (Stage 3) | §11 (Pillar A) |
| `verify_P2_nu_p_symbolic.py` | Symbolic structural proof of ν_p(J_k) = 0 (Stage 5) | §11 (Pillar A) |
| `verify_P2_row_sum_reduction.py` | Row-sum reduction to ε₁Δ_{l₁} + ε₂Δ_{l₂} + ε₃Δ_{l₃} = 0 (Stage 6) | §11 (Pillar A) |
| `verify_c_value.py` | Constant c = 2^(4m−1) · (−i) in the Im(π̄^(4m)) congruence | §11 (Pillar A) |
| `verify_rep_count_v4.py` | Representation count #reps = e + 1 for N = 2p^(2e) | §11 (Pillar A) |

### Group B: Phase C — Gaussian shadow and dichotomy (§4)

Located in `scripts/phase_C_supplementary/`.

| Script | Verifies | Paper Reference |
|--------|----------|-----------------|
| `verify_c12_vacuousness.py` | Q-linear independence of the specialised Gaussian shadow basis {Φ_S(ρ, t)}_{|S| odd}; exhaustive Case (D) vacuousness for all 448 non-complementary k=3 triples; verified for k ∈ {3, 4, 5, 6, 7} | Proposition `prop:c12-vacuousness` (§4.3) |
| `verify_branch1_dichotomy.py` | Theorem `thm:perfectsquare` identity and the both-roots-annihilate conclusion of Lemma `lem:Res-zero-dichotomy` on each irreducible factor of √Res for an explicit k=3 triple | Lemma `lem:Res-zero-dichotomy` (§4) |
| `verify_k3_cardinality.py` | Hamming-stratified cardinality (24 triples) and orbit decomposition (2 equivalence classes) for the all-one triples at k=3 | Lemmas `lem:k3-cardinality`, `lem:k3-orbit` (Appendix D) |

### Group C: Main case-analysis scripts

Located at the repository root.

| Script | Verifies | Paper Reference |
|--------|----------|-----------------|
| `verify_k3base.py` | Irreducibility of z₊ for k=3; enumeration of 24 all-one triples into 2 equivalence classes | Lemma 7.7 |
| `verify_birational.py` | Birational map between the k=2 quartic and the twist of `48.a3` | Remark B.1 |
| `verify_cor64.py` | Trace-back of E(3) rational points; no admissible t₁ > 1 | Corollary 6.4 |
| `verify_prop81.py` | Pell descent Q₂ → Q₁ | Proposition 8.1 |
| `verify_caseIV.py` | Degree argument for Case IV (non-complementary) | Proposition 5.5 |
| `verify_k3_via_thm75.py` | All 168 k=3 triples eliminated by Theorem 7.5 alone | §4.5 / §9 |
| `generate_k3_table.py` | LaTeX table for the k=3 enumeration appendix | Table 1 / Appendix D |
| `verify_appendixA.py` | 8 line-sums reduce to 4 symmetric-pair identities | Remark A (Step 2) |

## Usage

Run any script individually:

```
python3 verify_k3base.py
```

Each script prints a verification report and a final verdict (PASS/FAIL).

To run the full 17-script suite at once, use the included harness:

```
bash run_all.sh
```

This executes all 17 scripts in order and reports the cumulative PASS/FAIL count.

## Reproducibility

All scripts use `sympy` only (no external data or random seeds). Results are deterministic and fully reproducible.

## Citation

If you use these scripts, please cite both the paper and the code:

**Paper:**
> Fukui, H. (2026). *Nonexistence of 3×3 Magic Squares of Squares.* Zenodo. [https://doi.org/10.5281/zenodo.19601866](https://doi.org/10.5281/zenodo.19601866)

**Code (this repository):**
> Fukui, H. (2026). *hirokifukui/magic-square-proof: v1.1.0 — Triple-validation release.* Zenodo. [https://doi.org/10.5281/zenodo.19601792](https://doi.org/10.5281/zenodo.19601792)

### BibTeX

```bibtex
@misc{fukui2026magicsquares,
  author       = {Fukui, Hiroki},
  title        = {Nonexistence of 3{\texttimes}3 Magic Squares of Squares},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19601866},
  url          = {https://doi.org/10.5281/zenodo.19601866},
  note         = {Preprint}
}

@software{fukui2026magicsquaresproof,
  author       = {Fukui, Hiroki},
  title        = {{hirokifukui/magic-square-proof: v1.1.0 --- Triple-validation release}},
  year         = {2026},
  publisher    = {Zenodo},
  version      = {v1.1.0},
  doi          = {10.5281/zenodo.19601792},
  url          = {https://doi.org/10.5281/zenodo.19601792}
}
```

## Version history

- **v1.1.0** (2026-05) — Triple-validation release. Adds 6 Pillar A scripts (case k=1, split primes) and 3 Phase C supplementary scripts (Gaussian shadow, dichotomy, k=3 orbit decomposition). All 17 scripts PASS. Companion paper has undergone external adversarial review, self-adversarial re-reading, and machine verification (Triple Validation Cycle).
- **v1.0.0** (2026-04) — Initial release. 8 main case-analysis scripts.

## License

MIT License (see LICENSE file).

## Contact

Hiroki Fukui (ORCID: 0009-0008-7122-522X)
Research Institute of Criminal Psychiatry / Sex Offender Medical Center;
Department of Neuropsychiatry, Kyoto University
