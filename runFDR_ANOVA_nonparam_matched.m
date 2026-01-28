function res = runFDR_ANOVA_nonparam_matched(Data_v, Key_v, alpha)
% runFDR_ANOVA_nonparam_matched  Paired nonparametric per-level omnibus test
% with key-based matching across conditions ("days"), BH-FDR across levels,
% and paired Wilcoxon (signrank) post-hoc for significant levels.
%
% This is a safer variant of runFDR_ANOVA_nonparam for research use when the
% third-dimension sample ordering may not be identical across conditions.
% It matches subjects within each level using a stable key (e.g., prey id +
% background file) and runs paired tests on the matched set.
%
% Inputs:
%   Data_v : N_days x N_levels x num_sample numeric (can contain NaN)
%   Key_v  : cell array N_days x N_levels, where Key_v{i,j} is a vector of
%            subject identifiers (length n_i_j). Keys can be numeric, string,
%            char, or cellstr. Must uniquely identify a subject within a
%            given (day, level).
%   alpha  : FDR level (optional, default 0.05)
%
% Output:
%   res    : struct with fields similar to runFDR_ANOVA_nonparam, plus:
%            - matching_used (true)
%            - n_matched_per_level
%            - n_unique_per_day_level
%
% Notes:
%   - Per-level matching: for each level, the matched set is the intersection
%     of keys across all days.
%   - If your samples are independent across days, paired tests are not
%     appropriate; use Kruskal-Wallis instead.

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end

if ndims(Data_v) ~= 3
    error('Data_v must be N_days x N_levels x num_sample');
end

[N_days, N_levels, num_sample] = size(Data_v);

if ~iscell(Key_v) || ~isequal(size(Key_v), [N_days, N_levels])
    error('Key_v must be a cell array of size N_days x N_levels.');
end

% Top-level metadata
res = struct();
res.alpha = alpha;
res.N_days = N_days;
res.N_levels = N_levels;
res.num_sample = num_sample;
res.paired_design = true;
res.matching_used = true;
res.note = 'Per-level key-matched paired omnibus (signrank if N_days==2, friedman if N_days>=3), BH-FDR across levels, paired signrank post-hoc for FDR-significant levels.';

% Preallocate
pvals = nan(1, N_levels);
qvals = nan(1, N_levels);
signif_levels = false(1, N_levels);

omnibus_tbls = cell(1, N_levels);
omnibus_stats = cell(1, N_levels);
omnibus_test = cell(1, N_levels);
descriptives = cell(1, N_levels);
effect_info = cell(1, N_levels);

n_unique_per_day_level = nan(N_days, N_levels);
n_matched_per_level = nan(1, N_levels);

% Paper-friendly summaries
n_total_per_level = nan(1, N_levels);
n_complete_per_level = nan(1, N_levels);
n_excluded_per_level = nan(1, N_levels);
df_per_level = nan(1, N_levels);
chi2_per_level = nan(1, N_levels);
kendallW_per_level = nan(1, N_levels);
median_diff_per_level = nan(1, N_levels);
rank_biserial_per_level = nan(1, N_levels);
n_eff_per_level = nan(1, N_levels);

Xv_per_level = cell(1, N_levels);

for j = 1:N_levels
    % Build per-day key/value tables for this level
    keys_by_day = cell(1, N_days);
    vals_by_day = cell(1, N_days);

    for i = 1:N_days
        vals = squeeze(Data_v(i, j, :));
        keys = Key_v{i, j};

        % Normalize keys to string vector
        keys = i_keysToString(keys);

        % Truncate to min length (defensive)
        n_use = min(numel(vals), numel(keys));
        vals = vals(1:n_use);
        keys = keys(1:n_use);

        % Drop NaN values and empty keys
        valid = ~isnan(vals) & keys ~= "";
        vals = vals(valid);
        keys = keys(valid);

        % Make keys unique within this day/level to avoid ambiguous matching
        keys = i_makeUniqueKeys(keys);

        keys_by_day{i} = keys;
        vals_by_day{i} = vals;
        n_unique_per_day_level(i, j) = numel(unique(keys));
    end

    % Intersection across all days
    common = keys_by_day{1};
    for i = 2:N_days
        common = intersect(common, keys_by_day{i}, 'stable');
    end

    n_matched = numel(common);
    n_matched_per_level(j) = n_matched;

    % Total counts (for reporting)
    n_total = 0;
    for i = 1:N_days
        n_total = max(n_total, numel(keys_by_day{i}));
    end
    n_total_per_level(j) = n_total;

    % Build matched matrix subjects x days
    if n_matched == 0
        Xv = nan(0, N_days);
    else
        Xv = nan(n_matched, N_days);
        for i = 1:N_days
            [tf, loc] = ismember(common, keys_by_day{i});
            if ~all(tf)
                % Should not happen because common is intersection
                Xv(:, i) = NaN;
            else
                Xv(:, i) = vals_by_day{i}(loc);
            end
        end
    end

    Xv_per_level{j} = Xv;
    n_complete = size(Xv, 1);
    n_complete_per_level(j) = n_complete;
    n_excluded_per_level(j) = n_total_per_level(j) - n_complete;

    if isempty(Xv)
        descriptives{j} = struct('n_total', n_total_per_level(j), 'n_complete', n_complete, 'n_excluded', n_excluded_per_level(j), ...
            'median_per_day', nan(1, N_days), 'mean_per_day', nan(1, N_days), 'std_per_day', nan(1, N_days));
    else
        descriptives{j} = struct('n_total', n_total_per_level(j), 'n_complete', n_complete, 'n_excluded', n_excluded_per_level(j), ...
            'median_per_day', median(Xv, 1), 'mean_per_day', mean(Xv, 1), 'std_per_day', std(Xv, 0, 1));
    end

    if size(Xv, 1) < 2 || N_days < 2
        pvals(j) = NaN;
        omnibus_tbls{j} = [];
        omnibus_stats{j} = [];
        omnibus_test{j} = '';
        effect_info{j} = struct('kendallW', NaN, 'chi2', NaN, 'df', NaN, ...
            'median_diff', NaN, 'rank_biserial', NaN, 'n_total', n_total_per_level(j), 'n_complete', n_complete, 'n_eff', NaN);
        continue;
    end

    if N_days == 2
        % Omnibus: signrank (paired)
        try
            [p_pair, ~, stats_sr] = signrank(Xv(:,1), Xv(:,2));
            p = p_pair;
            meddiff = median(Xv(:,1) - Xv(:,2), 'omitnan');
        catch ME
            warning('signrank (omnibus) failed for level %d: %s', j, ME.message);
            p = NaN; stats_sr = []; meddiff = NaN;
        end

        omnibus_test{j} = 'signrank';
        pvals(j) = p;
        omnibus_tbls{j} = [];
        omnibus_stats{j} = stats_sr;

        % Effect: matched-pairs rank-biserial
        [r_rb, n_eff] = i_rankBiserialFromSignrankStats(Xv(:,1), Xv(:,2), stats_sr);

        effect_info{j} = struct('kendallW', NaN, 'chi2', NaN, 'df', 1, ...
            'median_diff', meddiff, 'rank_biserial', r_rb, 'n_total', n_total_per_level(j), 'n_complete', n_complete, 'n_eff', n_eff);

        df_per_level(j) = 1;
        chi2_per_level(j) = NaN;
        kendallW_per_level(j) = NaN;
        median_diff_per_level(j) = meddiff;
        rank_biserial_per_level(j) = r_rb;
        n_eff_per_level(j) = n_eff;
    else
        % Omnibus: Friedman (paired, >=3 conditions)
        try
            [p_f, tbl_f, stats_f] = friedman(Xv, 1, 'off');
            p = p_f; tbl = tbl_f; stats = stats_f;
        catch ME
            warning('friedman failed for level %d: %s', j, ME.message);
            p = NaN; tbl = []; stats = [];
        end

        omnibus_test{j} = 'friedman';
        pvals(j) = p;
        omnibus_tbls{j} = tbl;
        omnibus_stats{j} = stats;

        % Effect: Kendall''s W from chi2 if available
        [W, chi2, df] = i_kendallWFromFriedmanTbl(tbl, size(Xv,1), N_days);
        effect_info{j} = struct('kendallW', W, 'chi2', chi2, 'df', df, ...
            'median_diff', NaN, 'rank_biserial', NaN, 'n_total', n_total_per_level(j), 'n_complete', n_complete, 'n_eff', NaN);

        df_per_level(j) = df;
        chi2_per_level(j) = chi2;
        kendallW_per_level(j) = W;
        median_diff_per_level(j) = NaN;
        rank_biserial_per_level(j) = NaN;
        n_eff_per_level(j) = NaN;
    end
end

% BH-FDR across non-NaN p-values
valid_mask = ~isnan(pvals);
sp = pvals(valid_mask);
if any(valid_mask)
    sp = sp(:);
    [sp_sorted, sort_idx_rel] = sort(sp, 'ascend');
    m_valid = numel(sp_sorted);
    ranks = (1:m_valid)';

    raw_q = sp_sorted .* (m_valid ./ ranks);
    adj_q = raw_q;
    for t = m_valid-1:-1:1
        adj_q(t) = min(adj_q(t), adj_q(t+1));
    end
    adj_q = min(adj_q, 1);

    valid_idx = find(valid_mask);
    sort_idx_global = valid_idx(sort_idx_rel);
    qvals(sort_idx_global) = adj_q;

    thresh = (ranks / m_valid) * alpha;
    kf = find(sp_sorted <= thresh, 1, 'last');
    if ~isempty(kf)
        sig_global_sorted = sort_idx_global(1:kf);
        signif_levels(sig_global_sorted) = true;
    end
end

% Build output
res.pvals = pvals;
res.qvals = qvals;
res.signif_levels = signif_levels;
res.omnibus_test = omnibus_test;
res.n_total_per_level = n_total_per_level;
res.n_complete_per_level = n_complete_per_level;
res.n_excluded_per_level = n_excluded_per_level;
res.n_unique_per_day_level = n_unique_per_day_level;
res.n_matched_per_level = n_matched_per_level;
res.df_per_level = df_per_level;
res.chi2_per_level = chi2_per_level;
res.kendallW_per_level = kendallW_per_level;
res.median_diff_per_level = median_diff_per_level;
res.rank_biserial_per_level = rank_biserial_per_level;
res.n_eff_per_level = n_eff_per_level;

if N_days == 2
    res.effect_size_label = 'rank_biserial';
    res.effect_size_per_level = rank_biserial_per_level;
else
    res.effect_size_label = 'kendallW';
    res.effect_size_per_level = kendallW_per_level;
end

% Pairwise post-hoc for FDR-significant levels
res.level = struct([]);
for j = 1:N_levels
    res.level(j).omnibus_p = pvals(j);
    res.level(j).omnibus_test = omnibus_test{j};
    res.level(j).omnibus_stats = omnibus_stats{j};
    res.level(j).omnibus_tbl = omnibus_tbls{j};
    res.level(j).effect_info = effect_info{j};
    res.level(j).descriptives = descriptives{j};
    res.level(j).is_significant = signif_levels(j);
    res.level(j).n_matched = n_matched_per_level(j);
    res.level(j).pairwise = [];

    if ~res.level(j).is_significant
        continue;
    end

    Xv = Xv_per_level{j};
    n_common = size(Xv, 1);

    if N_days == 2
        [p_pair, h, stats_sr] = signrank(Xv(:,1), Xv(:,2));
        meddiff = median(Xv(:,1) - Xv(:,2), 'omitnan');
        [r_rb, n_eff] = i_rankBiserialFromSignrankStats(Xv(:,1), Xv(:,2), stats_sr);

        res.level(j).pairwise(1).day1 = 1;
        res.level(j).pairwise(1).day2 = 2;
        res.level(j).pairwise(1).p_signrank = p_pair;
        res.level(j).pairwise(1).h = h;
        res.level(j).pairwise(1).stats = stats_sr;
        res.level(j).pairwise(1).n_common = n_common;
        res.level(j).pairwise(1).median_diff = meddiff;
        res.level(j).pairwise(1).n_eff = n_eff;
        res.level(j).pairwise(1).rank_biserial = r_rb;
    else
        pair_counter = 0;
        for di = 1:N_days-1
            xi = Xv(:, di);
            for dk = di+1:N_days
                xk = Xv(:, dk);
                pair_counter = pair_counter + 1;
                [p_pair, h, stats_sr] = signrank(xi, xk);
                meddiff = median(xi - xk, 'omitnan');
                [r_rb, n_eff] = i_rankBiserialFromSignrankStats(xi, xk, stats_sr);

                res.level(j).pairwise(pair_counter).day1 = di;
                res.level(j).pairwise(pair_counter).day2 = dk;
                res.level(j).pairwise(pair_counter).p_signrank = p_pair;
                res.level(j).pairwise(pair_counter).h = h;
                res.level(j).pairwise(pair_counter).stats = stats_sr;
                res.level(j).pairwise(pair_counter).n_common = n_common;
                res.level(j).pairwise(pair_counter).median_diff = meddiff;
                res.level(j).pairwise(pair_counter).n_eff = n_eff;
                res.level(j).pairwise(pair_counter).rank_biserial = r_rb;
            end
        end
    end
end

end

function keys = i_keysToString(keys)
% Normalize keys to string vector
if isstring(keys)
    keys = keys(:);
    return;
end
if iscell(keys)
    try
        keys = string(keys(:));
    catch
        keys = string(cellfun(@(x) char(x), keys(:), 'UniformOutput', false));
    end
    return;
end
if isnumeric(keys)
    keys = string(keys(:));
    return;
end
if ischar(keys)
    % single char array -> 1 key
    keys = string({keys});
    return;
end
error('Unsupported key type: %s', class(keys));
end

function keys_u = i_makeUniqueKeys(keys)
% Ensure uniqueness by appending occurrence counts for duplicates.
keys = keys(:);
[uk, ~, ic] = unique(keys, 'stable');
counts = accumarray(ic, 1);
if all(counts == 1)
    keys_u = keys;
    return;
end

% Create suffixes for duplicates
idx_in_group = zeros(size(ic));
for k = 1:numel(uk)
    ids = find(ic == k);
    if numel(ids) > 1
        idx_in_group(ids) = 1:numel(ids);
    end
end

keys_u = keys;
dup = idx_in_group > 0;
keys_u(dup) = keys(dup) + "#" + string(idx_in_group(dup));
end

function [r_rb, n_eff] = i_rankBiserialFromSignrankStats(x1, x2, stats_sr)
% Compute matched-pairs rank-biserial correlation from signrank output.
r_rb = NaN;
n_eff = NaN;
try
    d = x1 - x2;
    d = d(~isnan(d));
    d = d(d ~= 0);
    n_eff = numel(d);
    if n_eff == 0
        r_rb = 0;
        return;
    end
    if ~isempty(stats_sr) && isfield(stats_sr, 'signedrank')
        totalRank = n_eff * (n_eff + 1) / 2;
        r_rb = (2 * double(stats_sr.signedrank) / totalRank) - 1;
    end
catch
    % keep defaults
end
end

function [W, chi2, df] = i_kendallWFromFriedmanTbl(tbl, n, k)
W = NaN; chi2 = NaN; df = k - 1;
try
    if isempty(tbl) || size(tbl, 1) < 2
        return;
    end
    hdr = tbl(1, :);
    chi2_col = find(strcmpi(hdr, 'Chi-sq') | strcmpi(hdr, 'Chi-Sq') | strcmpi(hdr, 'Chi^2') | strcmpi(hdr, 'Chi2'), 1, 'first');
    df_col = find(strcmpi(hdr, 'df'), 1, 'first');
    if ~isempty(chi2_col)
        chi2 = tbl{2, chi2_col};
    end
    if ~isempty(df_col)
        df = tbl{2, df_col};
    end
    if isnumeric(chi2) && ~isempty(chi2) && n > 0 && k > 1
        W = chi2 / (n * (k - 1));
    end
catch
    W = NaN; chi2 = NaN;
end
end
