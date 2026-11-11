# Sharing a placebo arm across seasons: platform trials for atopic dermatitis, a pre-specified simulation

Platform trials for atopic dermatitis: sharing a placebo arm across seasons, a pre-specified simulation from a Singularity v3 forecast.

**Paper:** [paper.pdf](paper.pdf) · **Page:** [singularitysd.com/research/platform.html](https://singularitysd.com/research/platform.html) · Preprint, not yet peer reviewed.

## Abstract

Platform trials test several treatments against one shared control arm, and are now common in oncology and infectious disease but rare in dermatology. Atopic dermatitis raises a specific concern: disease activity and placebo response vary with the seasons, and a platform that compares a late-entering treatment with controls recruited months earlier (non-concurrent controls) may be biased. We simulate a three-arm platform with staggered entry (months 0, 6 and 12), 150 patients per treatment and a binary response (15% on placebo), under four time-trend scenarios: none; seasons of about $\pm5$ percentage points; seasons plus a secular drift; and seasons acting on the placebo arm only. For the two late-entering arms we compare separate two-arm trials, the platform with concurrent controls, the platform with all controls pooled, and a logistic model adjusting for period. In 2000 simulated platforms per scenario, the false-positive rate stayed between 3.4% and 5.4% for every analysis in every scenario, and power for a response of 35% was 98–100%. The platform needed 828 patients instead of 900 and 378 on placebo instead of 450 ($-16\%$); pooling non-concurrent controls doubled the controls available to each comparison (328 vs 157) without measurable bias. In this setting the seasonal argument against borrowing controls is not supported, and the pre-specified effect was large enough to put power near its ceiling, which limits the comparison of efficiency.

## How it was made

This study started from the Singularity forecast **T3v3-2026-0046**: a method predicted to enter a new field. It was then tested once, with the plan fixed before looking at the results.

## What is in this repository

| file | what it is |
|---|---|
| `experiment.py` | reproduces the study end to end |
| `fig_platform.png` | figure |
| `main.tex` | LaTeX source of the paper |
| `paper.pdf` | the paper |
| `results.json` | every number reported in the paper |

## Reproduce

```bash
pip install matplotlib numpy
python experiment.py
```

It writes `results.json` and the figures next to the script.

## Cite

```bibtex
@misc{singularity_platform_trial_shared_placebo_2026,
  title  = {Sharing a placebo arm across seasons: platform trials for atopic dermatitis, a pre-specified simulation},
  author = {Singularity},
  year   = {2026},
  url    = {https://github.com/singularity-sd/platform-trial-shared-placebo}
}
```

## License

Code: MIT (`LICENSE`). Paper, text and figures: CC BY 4.0 (`LICENSE-paper.md`). Third-party data keep their own licence where stated.

---

Part of [Singularity](https://singularitysd.com): forecasting where science goes next, then testing it.
