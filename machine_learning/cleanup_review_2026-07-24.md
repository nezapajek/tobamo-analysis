# Machine Learning Folder — Review Notes (2026-07-24)

Review done while `train_model_pipeline.py --stage select --outdir model_selection`
was running in the background (started 2026-07-23, PID 498123 under the
`tobamo-model` conda env). At review time it was on iteration 3 of 5, fold 4
(last checkpoint 13:46) — on pace to finish iteration 4 the evening of
2026-07-24.

## Scripts present
All scripts referenced in `README.md` exist: `sample_refs.py`,
`getorfs_pairwise_aln.py`, `preprocess.py`, `train_model_pipeline.py`,
`predict_query_contigs.py`.

`analyze_feature_importance.py` also exists and is clean/self-contained, but
is not mentioned anywhere in the README's table of contents — a doc gap,
not a functional problem.

## Redundant code — mostly in `scripts/utils.py` (1865 lines)

- **Two name collisions where the second definition silently shadows the
  first** (Python keeps the last `def` with a given name):
  - `pairwise_refs` — defined at line 96 (fast Biopython `PairwiseAligner`)
    and again at line 1211 (MAFFT + multiprocessing). Only the second is
    ever reachable.
  - `make_distance_matrix` — defined at line 163 (has an empty-df guard and
    zeroes the diagonal) and again at line 290 (doesn't have either
    safeguard). Only the second is ever reachable.
  - `notebooks/reference_selection.ipynb` imports and calls both of these
    names, so it always silently gets the second, less-defensive version.
    The first versions (plus `calculate_identity_simple`, only used by the
    shadowed `pairwise_refs`) are dead code that looks live. **Worth fixing
    before relying on that notebook again** — pick one implementation per
    name and delete the other.
- **Copy-pasted triplet**: `prepare_bin_df`, `prepare_bin_df_refs`,
  `prepare_bin_df_refs_hist_test` are ~95% identical bodies — collapse to
  one function with a flag.
- **Near-duplicate pairs**: `make_dendrograms` / `make_dendrograms_genus`,
  `make_cluster_df` / `make_cluster_df_genus`, `predict_contigs` /
  `predict_contigs_binned_test`.
- **~85 lines of commented-out dead functions** at the bottom of the file;
  one (`train_lr_and_predict`) is a verbatim commented copy of the live
  function above it. Just delete these.
- **18 unused imports** (checked against every other script too, not just
  utils.py itself): `warnings`, `mpatches`, `sns`, `interp1d`, `AlignInfo`,
  `LeaveOneOut`, `KMeans`, `LinearRegression`, `mean_squared_error`,
  `adjusted_rand_score`, `classification_report`, `precision_recall_curve`,
  `roc_curve`, `ParameterGrid`, `cross_val_predict`, `train_test_split`,
  `LabelEncoder`, `calibration_curve`.
  - Note: some other imports that *look* unused inside `utils.py`
    (`LogisticRegression`, `SVC`, `GaussianNB`, etc.) are actually
    load-bearing — `train_model_pipeline.py` does `from utils import *`
    and relies on utils.py having imported them. Fragile pattern, but
    currently working — don't strip these without also fixing the import
    style in `train_model_pipeline.py`.
- **Latent bug**: `add_info_basic` does `missing_keys[0:10]` on a `set`
  (line 999) — would raise `TypeError` if that error path is ever hit,
  since sets aren't subscriptable. Should be `list(missing_keys)[0:10]`.
- `train_model_pipeline.py` carries a dead `model_selection_old()` (~135
  lines, starts line 108) that's never called — only `model_selection()`
  (line 244, wired to `--stage select`) is used, and it does correctly set
  `random_state` on every classifier, so the currently-running job is using
  the fixed function.

## Environments
Set up and working — proven by the background job actively running under
it. `tobamo-model` conda env exists; `mafft` and `orfipy` binaries are
present inside it.

One loose end: `analysis_environment.yml`'s pip section pins
`numpy==1.26.4` while `requirements.txt` pins `numpy==2.1.3`, and what's
actually installed is `numpy 2.2.6` / `biopython 1.84` (vs
`requirements.txt`'s `1.87`). Harmless today, but rebuilding the env from
either manifest right now wouldn't reproduce what's currently running —
worth reconciling the two files against `pip freeze` / `conda list` once
things are stable.

## Should everything be working?
Yes for the active pipeline. Two things worth fixing before trusting the
docs again:

1. README's "Notes on this migration" section says
   `results/model_selection/` is "empty pending the reproducibility-bug fix
   re-run" — stale now that the re-run is 3/5 iterations in. Update once
   tonight's run finishes.
2. Current `train_model_pipeline.py --stage final` saves the LR model as
   `lr_binned_{bin_num}_model.joblib`, but the checked-in artifact in
   `results/final_model/` is named `lr_histogram_10_model.joblib` — an
   older naming convention. `predict_query_contigs.py` only works against
   it because of an explicit fallback chain someone added for this exact
   drift (see `predict_query_contigs.py` around line 56-71). Works today,
   but it's a sign `final_model` predates the current code — regenerate it
   once the pipeline changes settle, or the fallback chain will need a
   fourth alias someday.

## Suggested cleanup order (not yet done)
1. Resolve the two shadowed-function name collisions in `utils.py` first —
   correctness risk, not just tidiness.
2. Delete dead code: `model_selection_old`, commented-out blocks, unused
   imports.
3. Collapse the `prepare_bin_df*` triplet and the `_genus` near-duplicates.
4. Fix the `missing_keys[0:10]` set-slicing bug.
5. Reconcile `analysis_environment.yml` / `requirements.txt` with the
   actually-installed package versions.
6. Regenerate `results/final_model/` under current naming, or rename the
   fallback chain's primary target.
7. Update the README migration notes once tonight's `model_selection` rerun
   completes.
