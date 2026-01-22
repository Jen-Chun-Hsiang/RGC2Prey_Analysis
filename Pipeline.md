# Pipeline (how to run this repo end-to-end)

This repo is primarily **MATLAB analysis + plotting** for cricket location prediction experiments/simulations.
The core inputs are per-run `.mat` files (typically in `Results/Mats/`) that contain prediction outputs
and training/test losses.

> Important: The **model training / inference code that produces the `.mat` files is not in this folder**.
> The MATLAB scripts here assume those `.mat` files already exist.

---

## 0) Requirements

- MATLAB (tested scripts assume relatively modern MATLAB; toolboxes are not explicitly declared).
- Access to the network paths referenced in the scripts (many scripts hard-code `\\storage1.ris.wustl.edu\...`).

Recommended:
- MATLAB Image Processing Toolbox (used by scripts that read PNG alpha channels and plot circles).

---

## 1) Folder conventions

These paths are used repeatedly:

- **Model output mats**: `Results/Mats/`
  - Expected to contain `.mat` files with names like:
    - `YYYYMMDDNN_cricket_<bg>_noise<noise>_cricket_location_prediction_200_prediction_error_with_path.mat`
    - `YYYYMMDDNN_cricket_<bg>disp<disp>_noise<noise>_cricket_location_prediction_200_prediction_error_with_path.mat`
  - (Exact pattern depends on the plotting script; see below.)

- **Figures out**: `Summary/Illustrator/`
  - Many scripts write `.eps` and `.png` here.

- **Acceptance-zone metadata** (optional but often enabled):
  - `selected_points_summary_body.mat`
  - `processed_cover_radius.mat`

---

## 2) (Optional but recommended) Build the “acceptance zone” correction file

Several analysis scripts can subtract an object-dependent tolerance radius from error traces.
This requires a mapping from **object ID → radius** stored in `processed_cover_radius.mat`.

### 2.1 Create `selected_points_summary_body.mat` (manual centers)

Run **GetCricketCoords_2.m**:

- You will be prompted to select a folder of PNGs.
- For each PNG, click the object center and confirm.
- The script writes a MAT file to the hard-coded save folder:
  - `selected_points_summary_body.mat`

Notes:
- The code extracts an integer `image_id` from filenames (first digit sequence in the filename).
- It also stores image width/height and converts to centered coordinates.

### 2.2 Create `processed_cover_radius.mat` (radius that covers X% of the object)

Run **CricketZoneProcess.m**:

- Reads alpha channels from images matching `processed_*.png`.
- Uses centers from `selected_points_summary_body.mat`.
- Computes `reach_zone_radius` for each object using **compute_reach_zone_radius.m** with:
  - `target_cover_perc = 0.85` (85% of alpha-defined object pixels)
- Saves:
  - `processed_cover_radius.mat` containing:
    - `file_index_list` (object IDs)
    - `processed_cover_radius` (radius per ID, in pixels)

Important:
- `CricketZoneProcess.m` includes a `keyboard` statement (intentional breakpoint).
  - Remove/continue past it if you want unattended batch processing.

---

## 3) Confirm your model-output `.mat` files exist

All main analysis scripts load `.mat` files from:

- `Results/Mats/`

and expect these variables at minimum:

- `test_losses` (vector, typically length 100)
- `training_losses` (vector, typically length 200)
- `all_paths` (numeric, `nTrials x (2*seqLen)` interleaved x/y)
- `all_paths_pred` (numeric, predicted paths; squeezes to `nTrials x seqLen x 2`)
- `all_paths_bg` (numeric, background paths; `nTrials x (2*seqLen)` interleaved)
- `all_path_cm` (numeric, “center of mass” path; `nTrials x (2*seqLen)` interleaved)
- `all_bg_file` (cell array of strings; used to classify background type)
- `all_id_numbers` (object/image IDs)
- `all_scaling_factors` (per-trial scaling factors across time)

If these fields are missing, most plotting scripts will error.

---

## 4) Run analysis/plot scripts (common workflows)

### Workflow A: Single dataset inspection (one run)

Run **PathPredictDistExtraction.m**:

- Set:
  - `exp_name` (e.g., `2025091802`)
  - `bg_type` (`'blend'` or `'grass'`)
  - `noise_level` (`'0.0'`, `'0.016'`, `'0.256'`, etc.)
  - `trial_id` for the single-trial view
- Outputs:
  - Single-trial path plot (true vs predicted vs center-of-mass path)
  - Time-series traces (error distance, fixed-shift RMS, velocity)
  - Mean±SEM error traces across trials

Use this when you want to verify one model run looks correct before doing comparisons.

### Workflow B: Compare multiple experiments at a fixed noise level

Run **MultiExperimentComparison.m**:

- Configure:
  - `exp_names` (cell array of experiment IDs)
  - `noise_level` (string)
  - `bg_type`
  - `trial_id`
- Outputs:
  - Figure with:
    - Ground-truth trajectory (optionally)
    - Prediction traces (optional)
    - Fixed RMS error trace
    - Velocity traces
  - Fixed RMS (and optional CM RMS) mean±SEM across time
  - Consistency checks for velocity (and later in script: path consistency)

This is the best script to reproduce “panel-style” plots with trajectory + error-vs-time + velocity.

### Workflow C: Summary plots across many experiments and many noise levels

Run **PlotSimulationResults_case.m**:

- Set `exp_id` to choose a pre-configured analysis set.
  - Each `exp_id` defines:
    - `Dates` (experiment IDs)
    - `Noise_level` list
    - `LSTM_layer_n` labels
    - `fname_pattern` (how to build MAT filenames)
    - `plot_line_ids` (which lines/conditions to plot)
- Common outputs:
  - Training loss vs epoch
  - Test loss vs noise
  - “Error dist.” (fixed-shift RMS, optionally acceptance-zone corrected) vs noise
  - Saves summary figure to `Summary/Illustrator/` as `.eps` and `.png`

This is the script that generates the **summary error-vs-noise** plots seen in many paper-style figures.

### Workflow D: 2D sweep (disparity × noise)

Run **PlotSimulationResults_case_disparity.m**:

- This script sweeps:
  - `Second_level` (e.g., disparity)
  - `Noise_level`
- Produces the same training/test/error summary plots but for a 2D parameter grid.

---

## 5) Notes on reproducibility

- Many scripts hard-code:
  - `fixed_shift = -9` (frames)
  - sampling rate implied by `t = (0:len-1)/100` → **100 Hz** (0.01 s per step)
  - pixel/degree scaling constants

- Background grouping:
  - Many scripts classify “simple contrast” trials as those where `all_bg_file` contains `gray_image`.

- “Correct zone” option:
  - When enabled (`is_correct_object_zone = true/1`), scripts subtract a per-time-step cutoff radius
    computed from `processed_cover_radius.mat` and `all_scaling_factors`.

---

## 6) Quick reference: which script does what

- **GetCricketCoords_2.m**: manual object center picking, saves `selected_points_summary_body.mat`.
- **CricketZoneProcess.m**: builds `processed_cover_radius.mat` from PNG alpha channels.
- **PathPredictDistExtraction.m**: single-run inspection + modular comparisons.
- **MultiExperimentComparison.m**: across-experiment comparisons, trajectory/error/velocity panels.
- **PlotSimulationResults_case.m**: curated “experiment sets” and summary curves vs noise.
- **PlotSimulationResults_case_disparity.m**: disparity×noise summary.
- **DetailsPathObjectAnalysis.m**: deeper diagnostics (loss vs background types; CM vs model; etc.).

