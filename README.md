# Diffusion Model Parameters Predict Long-Term Dementia Risk: Evidence from the EPIC-Norfolk Study

Code accompanying a study that applies a modified one-boundary diffusion model (DM) to Visual Sensitivity Test (VST) trial-level data from the EPIC-Norfolk cohort, and tests whether the resulting cognitive parameters add predictive value — over demographics and raw reaction times — for long-term dementia risk. The work is currently under peer review at *Psychology and Aging*.

## Overview

Participants (N = 7,084 after exclusions) completed two versions of the VST — a simple detection task and a more demanding complex/dynamic-dot version — during the third EPIC-Norfolk health check. Trial-level reaction times were modeled with a one-boundary diffusion process fit via simulation-based inference (BayesFlow), yielding per-participant estimates of drift rate (simple VST: *v*; complex VST: within-trial increase in drift, *Δv*), boundary separation (*a*), and non-decision time (*t0*). These parameters, along with demographic variables and raw RT quantiles, were used in (1) Cox proportional hazards models predicting incident dementia over up to 18 years of follow-up, and (2) supervised classification models (logistic regression, RBF SVM, random forest, Gaussian Naive Bayes, XGBoost) comparing predictive performance across feature sets.

## Repository structure

```
.
├── one_boundary_DM/
│   ├── simulator_simple.py        # Numba-jitted simulator + priors for the simple-VST model (v, a, t0)
│   ├── simulator_complex.py       # Numba-jitted simulator + priors for the complex-VST model (constant-drift and Δv variants)
│   ├── simplevst_modeling.ipynb   # BayesFlow workflow: trains the network, runs diagnostics, fits real HC3 simple-VST data
│   └── complexVST_model.ipynb     # BayesFlow workflow for the complex-VST (time-varying drift, Δv) model
├── Demo_and_survival/
│   ├── Functions.R                    # Shared helpers: outlier detection, demographic recoding, survival data prep, RT plotting
│   ├── demo_summary.R                 # Data cleaning, outlier elimination, merges DM parameters back in, builds demographic summary tables
│   ├── logistic_survivalanalyses.R    # Main Cox PH models (mean RT / RT quantiles / cognitive parameters) for simple and complex VST
│   └── subgroup_survival.R            # Re-runs the Cox PH models within age subgroups (<75 vs. ≥75 years)
└── ML_models/
    ├── ML_models_demo.ipynb       # Classification using demographic variables only
    ├── ML_models_rt.ipynb         # Classification using demographics + RT quantiles
    └── ML_models_cogmodel.ipynb   # Classification using demographics + DM cognitive parameters
```

## Analysis workflow

The pipeline moves back and forth between R and Python, since raw trial data is cleaned in R, modeled in Python, and the resulting parameters are merged back into R for survival analysis (and read into Python again for ML). Rough order of operations:

1. **`Demo_and_survival/demo_summary.R`** — loads the raw HC3 extracts, applies the hard RT cutoffs (250–2,500 ms simple VST; 300–10,000 ms complex VST), excludes participants with too few valid trials, and writes per-trial RT files (`HC3_simple.csv`, `HC3_complex.csv`) for the diffusion modeling step.
2. **`one_boundary_DM/simplevst_modeling.ipynb`** and **`complexVST_model.ipynb`** — define the BayesFlow simulator/adapter/workflow (DeepSet summary network + coupling flow), train the amortized network, run parameter-recovery and calibration diagnostics, then fit the trained network to the real per-participant RT data and export posterior summaries (`postsumHC3_simple.csv`, `postsumHC3_complex_inc.csv`).
3. **`Demo_and_survival/demo_summary.R`** (second pass) — merges the posterior parameter summaries and RT quantiles back into the demographic dataset, producing the master analysis file (`demog_HC3_full_py.csv`) used by every downstream script.
4. **`Demo_and_survival/logistic_survivalanalyses.R`** and **`subgroup_survival.R`** — fit and tabulate the Cox proportional hazards models (overall and age-stratified).
5. **`ML_models/ML_models_demo.ipynb`**, **`ML_models_rt.ipynb`**, **`ML_models_cogmodel.ipynb`** — run the five-fold cross-validated classification pipelines (with random oversampling for class imbalance) for each feature set and export performance metrics.

## Requirements

**Python (3.12.12)**
- bayesflow == 2.0.7
- keras (3.x, multi-backend)
- numpy, numba, pandas, scipy, seaborn, matplotlib
- scikit-learn == 1.7
- imbalanced-learn == 0.14.0
- xgboost == 3.1.1
- joblib (for pipeline caching)

**R (4.4.0)**
- survival (3.8-3), rms, survRM2, performance
- dplyr, tidyr, forcats, here, ggplot2, ggsurvfit, ggeffects, sjPlot, ggcorrplot
- caret, brms, pscl, yardstick, pROC, gee, broom, rsample

> The notebooks add `../Bayesflow/` to `sys.path` rather than importing an installed package. This reflects the secure, offline compute environment the original analyses were run in (BayesFlow was vendored into a sibling directory rather than installed via pip). If you're running this elsewhere, a standard `pip install bayesflow==2.0.7` and removing that `sys.path.append` line should work.

## Data availability

Raw EPIC-Norfolk data (the `data/` and `datasets/` directories referenced by these scripts, including trial-level RTs, demographic extracts, and dementia outcome variables) are **not included** in this repository and cannot be redistributed. Access is governed by a Data Transfer Agreement with the EPIC Management Committee; researchers can request access by emailing epic-norfolk@mrc-epid.cam.ac.uk. Model checkpoints, intermediate pickle files, and result tables produced by the notebooks/scripts (`checkpoints/`, `plots/`, `results_new/`) are likewise not tracked here.

## Reproducibility notes

- `demo_summary.R` and `logistic_survivalanalyses.R` both start with a hardcoded `setwd()` call pointing to the original secure-environment path — update this to your own working directory before running.
- Scripts assume the folder layout above (`data/`, `datasets/`, `checkpoints/`, `plots/`, `results_new/`) exists alongside the code; create these locally as needed.
- Random seeds are fixed where used (e.g., `RNG = np.random.default_rng(2024)` in the simulators, `random_state=42` in most sklearn models) for reproducibility of the modeling and ML steps, though the underlying data cannot be shared.

## Citation

If you use this code, please cite the associated manuscript (citation details to be finalized pending publication).


## Acknowledgments

Thanks to Tom Bishop (MRC Epidemiology Unit data management) for technical support throughout the project, and to Valentin Pratz for help setting up the Apptainer image used to run BayesFlow in the secure compute environment. And Stefan T. Radev for providing insights in ML models and his feedback on earlier versions of this manuscript.

## Contact

Tuba Hato — Psychology Institute, Heidelberg University
