# Method (processing, analysis, and parameters)

This document describes what the MATLAB code in this folder computes from the stored prediction outputs.

The attached figure you shared (movie → RGC simulation → CNN/LSTM → predicted (x,y) trajectory) matches the **conceptual model pipeline**.
In this repo, the MATLAB code focuses on the **post-hoc analysis** of outputs: trajectories, errors, velocity, and statistics.

---

## 1) Data model: what is in each prediction `.mat`

Most analyses load `.mat` files from `Results/Mats/` containing:

### 1.1 Core arrays

- `all_paths` (true trajectory)
  - Shape: `nTrials x (2*seqLen)`
  - Format: interleaved columns `[x1 y1 x2 y2 ... xT yT]`
  - Units: **normalized coordinates** (scaled later)

- `all_paths_pred` (predicted trajectory)
  - Stored shape varies; scripts use `squeeze(all_paths_pred)`
  - Expected after squeeze: `nTrials x seqLen x 2`

- `all_paths_bg` (background motion trajectory)
  - Shape: `nTrials x (2*seqLen)`
  - Used for velocity checks in `MultiExperimentComparison.m` and other diagnostics

- `all_path_cm` (center-of-mass trajectory)
  - Shape: typically `nTrials x (2*seqLen)`
  - Interpretation: a “center of mass” estimate derived from simulated RGC responses

### 1.2 Metadata per trial

- `all_bg_file` (cell array of strings)
  - Used to categorize trials, e.g.:
    - `contains(all_bg_file, 'gray_image')` → “simple contrast” / “gray” background

- `all_id_numbers`
  - Object/image identifiers (used to look up acceptance-zone radius)

- `all_scaling_factors`
  - Per-trial scaling values over time (used to map object size/distance changes)

### 1.3 Loss arrays

- `training_losses`
  - Training loss over epochs (often 200 epochs)

- `test_losses`
  - One value per test sample/trial (often 100 samples)

---

## 2) Coordinate reshaping and scaling

### 2.1 Deinterleaving true paths

True/background/CM paths are stored as interleaved x/y columns.
Scripts convert these to a 3D array using:

- `reshapeAllPaths(all_paths)` → `nTrials x seqLen x 2`

Implementation detail (see `reshapeAllPaths.m`):
- x is `all_paths(:, 1:2:end)`
- y is `all_paths(:, 2:2:end)`

### 2.2 Physical scaling (normalized → micrometers)

A common conversion used across scripts:

- `pix_to_um = 4.375`
- `sf_scale = 0.54`
- `real_dim = [120 90] * pix_to_um / sf_scale`

The code typically treats `real_dim` as the scale factor applied to normalized x/y:

- `true_scaled = true_path .* [real_dim_x real_dim_y]`

Interpretation:
- The normalized coordinate system is multiplied by the half-field dimensions in physical units.

### 2.3 Visual degrees conversion

Many plots use “degrees” by dividing microns by:

- `vis_deg_to_um = 32.5`

So a distance in micrometers becomes degrees via:

- `deg = um / 32.5`

In summary plots (e.g., `PlotSimulationResults_case.m`), this is applied via:

- `unit_factor = 1/32.5` when `is_degree == true`

---

## 3) Error metrics computed

### 3.1 Pointwise error distance (no shift)

For a single trial, the instantaneous Euclidean error is:

$$e_t = \sqrt{(x_t^{true}-x_t^{pred})^2 + (y_t^{true}-y_t^{pred})^2}$$

After scaling into micrometers, it becomes an error trace in physical units.

### 3.2 Fixed-shift RMS error trace

The primary error trace used in the “error distance” plots is a **fixed temporal shift** RMS error.

- Implemented in `calculateFixedShiftRMSError.m`
- Parameter:
  - `fixed_shift` is often set to `-9`

Interpretation of sign:
- `shift_val < 0` → predicted signal **lags**; prediction is shifted forward (or truth backward)

After shifting, a time-aligned error is computed for all valid time points.
The output is a vector `rms_error(t)` (one value per valid time step after shift).

### 3.3 Center-of-mass (CM) RMS error

Some scripts compute a CM-based error trace as well.
Example usage in `PlotSimulationResults_case.m`:

- True path is scaled with `real_dim`
- CM path is scaled with `pix_to_um/sf_scale` and compared with `real_dim = ones(1,2)`

(That asymmetry reflects that `all_path_cm` is stored in a different coordinate convention than `all_paths`.)

### 3.4 Velocity traces

Several scripts compute the velocity of:
- the cricket (true path)
- the background motion

Velocity is typically converted to deg/s by:
- multiply per-step displacement (deg) by sampling rate (100 Hz)

The 100 Hz assumption comes from repeated use of:
- `t = (0:len-1)/100`

---

## 4) Acceptance-zone / “correct zone” correction

Many summary scripts optionally apply a tolerance so that small errors within the object region are not penalized.

### 4.1 Building per-object radius: `processed_cover_radius.mat`

The per-object “radius” is created as follows:

1) `GetCricketCoords_2.m`
- Manually select a center point `(Xc, Yc)` for each object PNG
- Writes `selected_points_summary_body.mat` containing `summary_struct` with fields including:
  - `image_id`, `X`, `Y`, `Height`, `Width`

2) `CricketZoneProcess.m`
- Loads each `processed_<id>.png` alpha channel
- Defines object pixels as `alpha > 254`
- Computes the minimal radius around `(Xc,Yc)` that covers a target fraction:
  - `target_cover_perc = 0.85`
- Uses `compute_reach_zone_radius.m` (bisection search)
- Saves `processed_cover_radius.mat`:
  - `file_index_list` (object IDs)
  - `processed_cover_radius` (radius in pixels)

### 4.2 Mapping radius to a time-varying cutoff

During analysis, `acceptance_zone_radius.m` converts the stored pixel radius into a per-time-step cutoff:

- Looks up `radius_px` from `cover_radius` using the trial’s `cid = all_id_numbers(i)`
- Applies the same `fixed_shift` alignment to the **scaling factors**
- Converts to microns using:
  - `cut_off_um = radius_px * scale_shifted * 4.375 / 0.54`

### 4.3 Applying the correction

If `is_correct_object_zone = true`, scripts apply:

$$e_t^{corr} = \max(0, e_t - cut\_off_t)$$

This is done for both fixed-shift RMS and (when available) CM RMS.

Practical effect:
- Errors smaller than the object-size threshold at that time are treated as “within zone”.

---

## 5) Background categorization

Many scripts split trials by background type using:

- `is_simple_contrast = contains(all_bg_file, 'gray_image')`

Then they compute summary metrics separately for:
- “gray/simple” backgrounds
- “image/grass/blend” backgrounds (everything else)

This impacts which subset of trials are used in `bg_type = 'grass'|'simple'|'blend'` modes.

---

## 6) Statistics and significance testing

`PlotSimulationResults_case.m` calls:

- `runFDR_ANOVA_nonparam(DataP_v(plot_line_ids, :, :))`

Where `DataP_v` is typically:

- `N_days x N_levels x n_samples`
  - days ≈ conditions/experiments
  - levels ≈ noise levels
  - samples ≈ trials

`runFDR_ANOVA_nonparam.m` behavior:

- Omnibus test **per noise level** (paired design across conditions):
  - If `N_days == 2`: paired Wilcoxon signed-rank
  - If `N_days >= 3`: Friedman test
- Multiple-comparisons control:
  - Benjamini–Hochberg FDR across noise levels
- Post-hoc:
  - Paired sign-rank between day pairs for FDR-significant noise levels

Outputs include:
- raw p-values
- BH-adjusted q-values
- which noise levels are significant
- effect summaries (Kendall’s W for Friedman, median differences for signrank)

---

## 7) Parameter catalog (what to report in papers)

These are the parameters that most strongly affect results and must be reported:

### 7.1 Temporal alignment

- `fixed_shift` (commonly `-9` frames)
- Sample rate implied: `100 Hz` (from `t = (0:len-1)/100`)

### 7.2 Unit conversion

- `pix_to_um = 4.375`
- `sf_scale = 0.54`
- `real_dim = [120 90] * pix_to_um / sf_scale`
- `vis_deg_to_um = 32.5` (used to convert um → degrees)

### 7.3 Acceptance-zone correction

- `target_cover_perc = 0.85`
- `alpha > 254` defines the object mask
- `is_correct_object_zone` enabled/disabled

### 7.4 Subsetting by background

- definition of “simple contrast”: `contains(all_bg_file, 'gray_image')`
- `bg_type` logic in summary scripts:
  - `'blend'`: include all
  - `'grass'`: exclude simple-contrast
  - `'simple'`: only simple-contrast

### 7.5 Experiment set definitions

`PlotSimulationResults_case.m` uses a large `switch exp_id` that defines:
- experiment IDs (`Dates`)
- noise levels
- filename patterns
- condition labels (`LSTM_layer_n`)

For any figure produced with that script, you should record:
- `exp_id`
- `Dates`, `Noise_level`, `plot_line_ids`
- `bg_type` and `is_correct_object_zone`

---

## 8) `PlotSimulationResults_case.m`: comparing experimental conditions (core results path)

Most of the paper-style summary results (including the figure you attached) are derived from `PlotSimulationResults_case.m`.
This script is the “batch summarizer” that compares **multiple experimental conditions** across **multiple noise levels** using a consistent error metric.

### 8.1 Inputs and naming conventions

The script loops over:

- `Dates` (called `N_days` in code): each entry is one experimental condition/run (e.g., ON vs OFF, temporal filter variants, surround strength, interocular distance, etc.).
- `Noise_level` (called `N_levels` in code): the x-axis levels.

The file to load for each (Date, Noise) cell is constructed from:

- `Folder_Name = Results/Mats` (hard-coded network path)
- `Tag = '_cricket_location_prediction_200_prediction_error_with_path'`
- `fname_pattern` selected by `exp_id` (e.g., includes optional background prefix and/or disparity fields)

and loaded variables include:

- losses: `training_losses`, `test_losses`
- paths: `all_paths`, `all_paths_pred`, `all_path_cm`
- metadata: `all_id_numbers`, `all_scaling_factors`, `all_bg_file`

### 8.2 Condition definitions (`exp_id`)

The `switch exp_id` block is effectively the “figure configuration”:

- Defines which experiments are grouped as conditions: `Dates`
- Defines the set of noise levels: `Noise_level`
- Defines human-readable labels: `LSTM_layer_n`
- Chooses `plot_line_ids` (which conditions are shown in the final figure)
- Chooses `fname_pattern` (how to build filenames for that set)

For reproducibility, a figure generated by this script should always report:

- `exp_id`
- the exact `Dates` list
- the exact `Noise_level` list
- `plot_line_ids`
- `bg_type` selection (see below)
- `fixed_shift`
- whether acceptance-zone correction was used (`is_correct_object_zone`)

### 8.3 Per-trial processing pipeline

For each loaded `.mat`:

1) **Reshape paths**
   - `all_paths_r = reshapeAllPaths(all_paths)` → `nTrials x seqLen x 2`
   - `all_paths_pred_r = squeeze(all_paths_pred)`
   - `all_path_cm` is reshaped similarly (with a fallback to the regular path if missing)

2) **Determine background category per trial**
   - `is_simple_contrast = contains(all_bg_file, 'gray_image')`

3) **Compute fixed-shift RMS error trace per trial**
   - For each trial `ii`:
     - `fixed_rms = calculateFixedShiftRMSError(true_path_trial, pred_path_trial, fixed_shift, real_dim)`
     - Optionally also compute a CM-based error trace (`fixed_cm_rms`) with the script’s CM scaling.

4) **(Optional) acceptance-zone correction**
   When `is_correct_object_zone` is enabled, the script subtracts a time-varying cutoff derived from:
   - `processed_cover_radius.mat` (object radius per `all_id_numbers(ii)`)
   - `all_scaling_factors(ii, 50:end)` aligned with the same `fixed_shift`
   - conversion `4.375/0.54`
   and then applies `max(0, error - cutoff)`.

5) **Unit conversion**
   When `is_degree` is enabled, error traces are converted from micrometers to degrees using:
   - `unit_factor = 1/32.5`

### 8.4 Aggregation (how each summary number is formed)

The script produces several summary tensors; the key ones for condition comparisons are:

- `Data_m(i,j)` / `Data_s(i,j)`
  - mean and SEM of `test_losses` for condition `i` at noise level `j`.

- `DataP_m(i,j)` / `DataP_s(i,j)`
  - **mean fixed-shift RMS error**, aggregated across time and across trials, for condition `i` at noise `j`.
  - The exact trial subset used depends on `bg_type`:
    - `'blend'`: all trials
    - `'grass'`: only trials where `~is_simple_contrast`
    - `'simple'`: only trials where `is_simple_contrast`

- `DataP_v(i,j,:)`
  - per-trial summary values used for statistics.
  - The script stores (per condition, per noise level) the **mean error over time** for each trial.
  - This is the quantity passed to `runFDR_ANOVA_nonparam`.

- `DataCM_*` equivalents
  - analogous aggregates for the center-of-mass error traces.

### 8.5 Plots produced (what the figure panels represent)

The default figure created by `PlotSimulationResults_case.m` (at the end of the main loop) contains:

- Training loss vs epoch
  - uses `training_losses(1:n_epoch)`; typically `n_epoch=200`.

- Test loss vs noise
  - `Data_m` with SEM `Data_s`.

- Error distance vs noise
  - `DataP_m` (and optionally `DataCM_m`) with SEM.
  - The y-axis label is degrees if `is_degree=true`.

### 8.6 Statistical comparison across conditions

At the end, the script runs:

- `runFDR_ANOVA_nonparam(DataP_v(plot_line_ids,:,:))`

This compares the selected conditions (`plot_line_ids`) at each noise level using paired nonparametric tests,
then applies BH-FDR across noise levels (see Section 6).

Interpretation guidance:
- The unit being tested is **trial-averaged (over time) fixed-shift RMS**.
- Pairing assumes trials are aligned across conditions in the third dimension.
  - This is appropriate only if the same trial IDs/objects/backgrounds appear in the same order across conditions.

---

## 9) `MultiExperimentComparison.m`: example trace plots (simulation vs prediction)

In addition to the population/summary plots in `PlotSimulationResults_case.m`, the figure also uses **example-trace panels**
generated by `MultiExperimentComparison.m`. This script is designed to:

- Load multiple conditions (or multiple disparities) at a fixed noise level
- Plot an **example trial** showing the ground-truth trajectory and (optionally) prediction traces
- Plot the corresponding **fixed-shift RMS error trace** and **velocity trace**
- Verify that trials are comparable across conditions (velocity and path consistency checks)

### 9.1 Comparison modes: experiments vs disparity sweeps

The script supports two mutually exclusive modes:

1) **Experiment comparison mode** (default)
   - Set `disparity_sets = {}`
   - Provide `exp_names = {'YYYYMMDDNN', ...}`
   - Each entry in `exp_names` becomes one “condition” in the comparison.

2) **Disparity sweep mode**
   - Set `disparity_sets = {'0.0','3.0','6.0','12.0',...}`
   - Provide `exp_names = {'BASE_EXPERIMENT'}` (only the first is used)
   - Each entry in `disparity_sets` becomes one “condition”, and the loader uses filenames of the form:
     - `..._disp<disp>_noise<noise>_...`

Internally, both modes call `loadDataset(...)` which standardizes output shapes using `reshapeAllPaths` and `squeeze`.

### 9.2 Key parameters that control the trace-example figure

These are the parameters most relevant to the trace panels:

- `trial_id`
  - which trial index is shown.

- `noise_level` and `bg_type`
  - select which `.mat` files to load.

- `fixed_shift` (commonly `-9`)
  - used consistently for both RMS error alignment and velocity alignment.

- `is_plot_ground_truth`
  - if true, plots ground truth trajectory from a reference condition.

- `is_plot_pred_trace`
  - if true, plots the prediction trace(s) for the condition(s) listed in `disp_trajectory_id`.

- `disp_trajectory_id`
  - list of which conditions’ predictions are overlaid.

- `is_y_axis_flip`
  - when true, `set(gca,'YDir','reverse')` for movie-style coordinate conventions.

### 9.3 Units and scaling used in the example-trace panels

The trajectory, error, and velocity are plotted in **visual degrees** when `is_visual_degree = 1`.

The script uses:

- `real_dim = [120 90] * 4.375 / 0.54`
  - converts normalized path coordinates into micrometers.

- `vis_deg_to_um = 32.5`
- `vis_scale = 1/vis_deg_to_um` (when visual-degree mode is on)

So distances in micrometers become degrees via `deg = um * vis_scale`.
Velocities are plotted as deg/s by additionally multiplying by the sampling rate:

- Sampling rate is assumed to be **100 Hz** (from `t = (0:len-1)/100`).
- Velocity magnitude (deg/s) is computed as `velocity_um_per_step * 100 * vis_scale`.

### 9.4 How trajectories are drawn (the “example trace” panel)

The script renders trajectories as line segments with a **light→dark gradient** to indicate time.

- Controlled by `alpha_start` and `alpha_end`.
- For each segment $j$ (from point $j$ to $j+1$), the segment color is interpolated between white and the condition’s base color.
- It also marks the start and end points.

The plot includes:

- a boundary rectangle based on the full `real_dim` extent
- optional 10-degree scale bars (when visual-degree mode is on)

### 9.5 Which “prediction trace” is being shown

There are two distinct ways the script can define the displayed “prediction” trajectory:

1) **Direct model output** (recommended for “model prediction” plots)
   - Set `disparity_shift_trace = 'none'`
   - The script uses `exp.all_paths_pred_r(trial_id,:,:)` (from the `.mat`).

2) **Disparity-shift reconstructed trace** (a visualization/geometry transform)
   - Set `disparity_shift_trace` to `'left'` or `'right'`
   - The script calls `reconstruct_top_img_positions_shifted(...)` using:
     - the **true path** (`exp.all_paths_r`) and
     - `scaling_factors` and
     - `interocular_dist` (in cm)
   - This produces two shifted trajectories intended to approximate left/right eye projections based on disparity.

Conceptually, the disparity reconstruction:

- maps scaling factors to distance (assumed linear between 21 cm and 4 cm),
- computes binocular disparity $2\arctan((IOD/2)/distance)$ in degrees,
- converts degrees to pixels and then to normalized units,
- shifts the x-coordinate by ±(disparity/2).

Important interpretation note:
- This reconstruction is **not** the trained model’s prediction; it is a disparity-based re-projection derived from the ground truth + scaling factors.
- Use it only when the intended figure panel is explicitly about left/right eye disparity geometry.

### 9.6 Fixed-shift RMS error trace (example trial)

For the chosen `trial_id`, the script computes:

- `fixed_rms = calculateFixedShiftRMSError(true_path, pred_path, fixed_shift, real_dim)`

Then plots:

- `t = (0:rms_len-1)/100` seconds
- `fixed_rms * vis_scale` in degrees

### 9.7 Velocity trace (example trial) and consistency validation

The script computes per-step velocity magnitudes using `calculateVelocity(true_path, bg_path, fixed_shift, real_dim)`.
It then:

- plots cricket and background velocity traces (deg/s)
- validates that cricket velocities (and background velocities) are identical across conditions when comparing experiments intended to share the same underlying stimulus

If velocities differ, the script warns/errors because cross-condition comparisons can be invalid if the underlying trajectories differ.

### 9.8 Across-trial summaries inside `MultiExperimentComparison.m`

Beyond the example trace panels, the script also computes across-trial summaries:

- **Fixed RMS and CM RMS vs time** (mean ± SEM across trials)
  - For each condition, compute `fixed_rms(t)` per trial and aggregate across all trials.

- **Velocity vs time** (mean ± SEM across trials)
  - Computed on up to the first 50 trials for efficiency.

- **Path consistency checks**
  - Overlays the first few true paths and background paths across conditions.
  - Quantifies average path differences for a subset of trials and reports whether they are effectively identical.

These checks are critical when the interpretation relies on “same stimulus / different model condition”.

---

## 10) Mapping to the figure you shared (practical guide)

Based on the scripts in this folder:

- “Trajectory panel (ground truth + predicted traces)” and “velocity panel”:
  - produced by `MultiExperimentComparison.m` (especially when `is_plot_ground_truth=true` and `is_plot_pred_trace=true`)

- “Error distance vs time” panels:
  - fixed-shift RMS traces from `calculateFixedShiftRMSError.m` (used across `PathPredictDistExtraction.m` and `MultiExperimentComparison.m`)

- “Error distance vs noise level” panels:
  - summary curves from `PlotSimulationResults_case.m`

If you want, tell me **which exact exp_id** you used for that figure, and I can add a dedicated “Figure recipe” section
that lists the exact script + parameter block to run to regenerate each panel.
