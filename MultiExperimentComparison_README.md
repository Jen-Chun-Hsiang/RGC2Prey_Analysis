# Multi-Experiment Comparison Script Documentation

## Overview

The `MultiExperimentComparison.m` script is designed to compare cricket path prediction data across multiple experiments. It extends the functionality of `PathPredictDistExtraction.m` to provide comprehensive multi-experiment analysis.

## Features

### 1. **Multi-Experiment Loading and Configuration**
- Configurable list of experiment names
- Automatic data loading for all specified experiments
- Consistent parameter settings across experiments

### 2. **Single Trial Visualization Across Experiments**
- **Path Trajectories**: Compare true and predicted paths for a specific trial across all experiments
- **Fixed RMS Error**: Time-series comparison of prediction accuracy
- **Fixed CM RMS Error**: Center-of-mass based error comparison
- **Velocity Comparison**: Cricket velocity across experiments
- Different colors for each experiment for easy identification

### 3. **Fixed RMS Error Time Series Analysis**
- Mean ± SEM plots for fixed RMS error across time
- Separate analysis for standard RMS and center-of-mass RMS
- Statistical comparison across all trials in each experiment

### 4. **Velocity and Path Consistency Verification**
- **True Cricket Velocity**: Verify consistency across experiments
- **Background Velocity**: Check if background motion is consistent
- **Path Consistency Check**: Visual and quantitative verification that true paths are identical
- **Background Path Verification**: Ensure background paths match across experiments

## Usage

### Basic Configuration

```matlab
% Configure experiments to compare
exp_names = {'2025091502', '2025091506', '2025091802'}; % Modify as needed

% Set trial ID for single trial visualization
trial_id = 96;

% Configure other parameters
bg_type = 'blend'; % or 'grass'
noise_level = '0.0'; % Fixed noise level for comparison
```

### Key Parameters

- **`exp_names`**: Cell array of experiment names to compare
- **`trial_id`**: Specific trial ID to visualize across experiments
- **`bg_type`**: Background type ('blend' or 'grass')
- **`noise_level`**: Noise level for comparison (string)
- **`fixed_shift`**: Fixed shift value for correlation analysis (-9)

## Output Figures

### Figure 1: Single Trial Comparison
- **Subplot 1**: Path trajectories (true vs predicted for each experiment)
- **Subplot 2**: Fixed RMS error comparison
- **Subplot 3**: Fixed CM RMS error comparison  
- **Subplot 4**: Velocity comparison

### Figure 2: Fixed RMS Error Comparison
- **Subplot 1**: Mean fixed RMS error across all trials with SEM
- **Subplot 2**: Mean fixed CM RMS error across all trials with SEM

### Figure 3: Velocity and Path Consistency Verification
- **Subplot 1**: True cricket velocity across experiments
- **Subplot 2**: Background velocity across experiments
- **Subplot 3**: True paths consistency verification (first 5 trials)
- **Subplot 4**: Background paths consistency verification (first 5 trials)

## Consistency Verification

The script performs quantitative verification to ensure proper comparison:

### Path Consistency Checks
- Compares sequence lengths across experiments
- Calculates average path differences for true paths
- Calculates average path differences for background paths
- Reports whether paths are identical or different

### Expected Results
For proper pairwise comparison:
- ✓ True paths should be IDENTICAL across experiments
- ✓ Background paths should be IDENTICAL across experiments
- Sequence lengths should match

If paths differ, this indicates different experimental conditions that may affect comparison validity.

## Dependencies

The script requires the following external functions:
- `loadDataset()` (included in script)
- `calculateFixedShiftRMSError.m`
- `reshapeAllPaths.m`
- `acceptance_zone_radius.m` (for zone correction, if used)
- `processed_cover_radius.mat` (coverage data file)

## Color Scheme

The script uses a predefined color scheme for up to 8 experiments:
1. Red (1.0, 0.2, 0.2)
2. Green (0.2, 0.8, 0.2)
3. Blue (0.2, 0.2, 1.0)
4. Orange (1.0, 0.5, 0.0)
5. Magenta (0.8, 0.2, 0.8)
6. Cyan (0.0, 0.8, 0.8)
7. Gray (0.5, 0.5, 0.5)
8. Yellow (0.8, 0.8, 0.2)

For more than 8 experiments, additional colors are automatically generated.

## Example Output

The script will print progress information and consistency verification results:

```
Loading data for 3 experiments...
Loading experiment 1/3: 2025091502
  Loading: 2025091502_cricket_blend_noise0.0_cricket_location_prediction_200_prediction_error_with_path.mat
  - Loaded 200 trials with sequence length 100

=== PATH CONSISTENCY VERIFICATION ===
Comparing 2025091502 vs 2025091506:
  Sequence lengths: 100 vs 100
  Avg true path difference (first 10 trials): 0.000000
  Avg background path difference (first 10 trials): 0.000000
  ✓ True paths are IDENTICAL
  ✓ Background paths are IDENTICAL
```

## Customization

### Adding More Experiments
Simply add experiment names to the `exp_names` cell array:

```matlab
exp_names = {'2025091502', '2025091506', '2025091802', 'new_experiment'};
```

### Changing Analysis Parameters
Modify the configuration section:

```matlab
trial_id = 50; % Different trial for visualization
fixed_shift = -5; % Different shift value
bg_type = 'grass'; % Different background
noise_level = '0.016'; % Different noise level
```

### Adjusting Plot Limits
Modify the `ylim()` calls throughout the script:

```matlab
ylim([0 20]); % For different RMS error ranges
ylim([0 200]); % For different velocity ranges
```

The script provides a comprehensive framework for comparing multiple experiments while maintaining the core functionality and analysis methods from the original `PathPredictDistExtraction.m` script.