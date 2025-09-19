# Your Project Name

Brief description of your MATLAB project.

## Requirements
- MATLAB R20XX or later
- Required toolboxes (list them)

## Setup
1. Clone this repository
2. Open MATLAB
3. Navigate to the project folder
4. Run the main script

## Usage
Describe how to use your scripts
 
## compareTraceVelocity utility
A helper `compareTraceVelocity(traceA, traceB, ...)` was added to compare the
moving velocity between two 2D traces. Each trace should be an N-by-2 array
of `[x y]` positions sampled at equal time intervals (or provide a `Time`
vector). The function returns a struct containing velocity vectors, speed
time-series, angular differences, cosine similarity of velocity vectors,
and cross-correlation of speeds.

Quick example:

```matlab
res = compareTraceVelocity(traceA, traceB, 'Fs', 30, 'Smooth', 5);
fprintf('Speed correlation: %.3f\n', res.speedCorr);
fprintf('Mean angular diff (deg): %.2f\n', rad2deg(res.angDiffMean));
fprintf('Mean cosine similarity: %.3f\n', res.cosSimMean);
```

Interpretation tips:
- `speedCorr` shows how similar the magnitudes of velocity are over time.
- `angDiff` (and `angDiffMean`) describe the signed turning difference from
	`traceA` to `traceB` at each sample (in radians); convert to degrees for
	readability. Values near 0 mean the two traces move in the same direction.
- `cosSim` is the cosine of the angle between the instantaneous velocity
	vectors; values near 1 mean they point the same way, near -1 mean opposite.
- Use the `crossCorr` output to find lags where speed time-series align.

There's a demo script `test_compareTraceVelocity.m` that builds synthetic
examples and produces a quick figure `compareTraceVelocity_demo.png`.
