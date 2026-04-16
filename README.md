[![DOI (paper)](https://zenodo.org/badge/DOI/10.5281/zenodo.19601866.svg)](https://doi.org/10.5281/zenodo.19601866)
[![DOI (code)](https://zenodo.org/badge/DOI/10.5281/zenodo.19601792.svg)](https://doi.org/10.5281/zenodo.19601792)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Verification Scripts for "Nonexistence of 3×3 Magic Squares of Squares"

This repository contains the SymPy-based verification scripts accompanying the paper:

**Hiroki Fukui**, *Nonexistence of 3×3 Magic Squares of Squares* (2026).

The paper proves that no 3×3 magic square with all nine entries distinct perfect squares exists. The proof is entirely algebraic; these scripts provide independent computational verification of specific lemmas and propositions.

## Requirements

- Python 3.9+
- SymPy (install: `pip install sympy`)

All scripts run in under one minute on standard hardware.

## Scripts

| Script | Verifies | Paper Reference |
|--------|----------|-----------------|
| `verify_k3base.py` | Irreducibility of z₊ for k=3; enumeration of 24 all-one triples into 2 equivalence classes | Lemma 7.7 |
| `verify_birational.py` | Birational map between the k=2 quartic and the twist of 48.a3 | Remark B.1 |
| `verify_cor64.py` | Trace-back of E(3) rational points; no admissible t₁ > 1 | Corollary 6.4 |
| `verify_prop81.py` | Pell descent Q₂ → Q₁ | Proposition 8.1 |
| `verify_caseIV.py` | Degree argument for Case IV (non-complementary) | Proposition 5.5 |
| `verify_k3_via_thm75.py` | All 168 k=3 triples are eliminated by Theorem 7.5 alone | §4.5 / §9 |
| `generate_k3_table.py` | Generates the LaTeX table in Appendix D | Table 1 / Appendix D |
| `verify_appendixA.py` | 8 line-sums reduce to 4 symmetric-pair identities | Remark A (Step 2) |

## Usage

Run any script individually:

```
python3 verify_k3base.py
```

Each script prints a verification report and a final verdict (PASS/FAIL).

## Reproducibility

All scripts use `sympy` only (no external data or random seeds). Results are deterministic and fully reproducible.

## Citation

If you use these scripts, please cite both the paper and the code:

**Paper:**
> Fukui, H. (2026). *Nonexistence of 3×3 Magic Squares of Squares.* Zenodo. [https://doi.org/10.5281/zenodo.19601866](https://doi.org/10.5281/zenodo.19601866)

**Code (this repository):**
> Fukui, H. (2026). *hirokifukui/magic-square-proof: v1.0.0 — Initial release.* Zenodo. [https://doi.org/10.5281/zenodo.19601792](https://doi.org/10.5281/zenodo.19601792)

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
  title        = {{hirokifukui/magic-square-proof: v1.0.0 --- Initial release}},
  year         = {2026},
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.19601792},
  url          = {https://doi.org/10.5281/zenodo.19601792}
}
```

## License

MIT License (see LICENSE file).

## Contact

Hiroki Fukui (ORCID: 0009-0008-7122-522X)
Research Institute of Criminal Psychiatry / Sex Offender Medical Center;
Department of Neuropsychiatry, Kyoto University
