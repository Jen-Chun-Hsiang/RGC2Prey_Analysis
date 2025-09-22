function [reach_zone_radius, covered_frac] = compute_reach_zone_radius(alpha, Xc, Yc, target_cover_perc, opts)
% compute_reach_zone_radius  Find radius around center that covers target fraction of object
% 
% Usage:
%   r = compute_reach_zone_radius(alpha, Xc, Yc, target_cover_perc)
%   [r, frac] = compute_reach_zone_radius(..., opts)
%
% Inputs:
%   alpha               - 2D alpha channel image (numeric, e.g. uint8)
%   Xc, Yc              - center coordinates in image pixel space (x=cols, y=rows)
%   target_cover_perc   - fraction in [0,1] of object pixels to cover
%   opts (optional)     - struct with fields:
%                          .tol (default 0.5)    : radius tolerance in pixels
%                          .maxIter (default 50) : maximum bisection iterations
%                          .verbose (default false)
%
% Outputs:
%   reach_zone_radius   - radius (in pixels) such that at least target_cover_perc of
%                         object pixels lie within distance <= radius from (Xc,Yc)
%   covered_frac        - achieved covered fraction for returned radius
%
% Note: object pixels are defined as alpha > 254.

if nargin < 5, opts = struct(); end
if ~isfield(opts,'tol'), opts.tol = 0.5; end
if ~isfield(opts,'maxIter'), opts.maxIter = 50; end
if ~isfield(opts,'verbose'), opts.verbose = false; end

assert(~isempty(alpha) && ndims(alpha)==2, 'alpha must be a 2D image');
assert(~isnumeric(target_cover_perc) || (target_cover_perc>=0 && target_cover_perc<=1), 'target_cover_perc in [0,1]');

alpha = double(alpha);  % work in double for comparisons
objMask = alpha > 254;
total_obj = nnz(objMask);
if total_obj == 0
    error('No object pixels found (alpha > 254).');
end

[rows, cols] = size(alpha);
% build distance matrix from center (Xc,Yc) to every pixel center
[xGrid, yGrid] = meshgrid(1:cols, 1:rows);
dist = hypot(xGrid - Xc, yGrid - Yc);

% compute maximum distance to any object pixel -> upper bound
objDists = dist(objMask);
r_low = 0;
r_high = max(objDists);

% trivial cases
if target_cover_perc <= 0
    reach_zone_radius = 0;
    covered_frac = 0;
    error('No object pixels found (alpha > 254).');
    return;
end
if target_cover_perc >= 1
    reach_zone_radius = r_high;
    covered_frac = 1;
    error('No object pixels found (alpha > 254).');
    return;
end

iter = 0;
while (r_high - r_low) > opts.tol && iter < opts.maxIter
    iter = iter + 1;
    r_mid = (r_low + r_high) / 2;
    covered = nnz(objMask & (dist <= r_mid));
    frac = covered / total_obj;
    if opts.verbose
        fprintf('iter %d: r_low=%.3f r_mid=%.3f r_high=%.3f frac=%.4f\n', iter, r_low, r_mid, r_high, frac);
    end
    if frac >= target_cover_perc
        % can shrink radius
        r_high = r_mid;
    else
        % need larger radius
        r_low = r_mid;
    end
end

reach_zone_radius = r_high;
covered_frac = nnz(objMask & (dist <= reach_zone_radius)) / total_obj;

end