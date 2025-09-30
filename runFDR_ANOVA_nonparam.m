function res = runFDR_ANOVA_nonparam(Data_v, alpha)
% runFDR_ANOVA_nonparam  Paired nonparametric per-level omnibus test,
% BH-FDR across levels, paired Wilcoxon (signrank) post-hoc for significant levels.
%
% Behavior change:
%  - If N_days == 2: use paired signrank as the omnibus test per level.
%  - If N_days >= 3: use Friedman as before.
%
% Inputs:
%   Data_v : N_days x N_levels x num_sample (numeric, can contain NaN)
%   alpha  : FDR level (optional, default 0.05)
%
% Notes:
%   - Paired design: third-dimension index is subjects repeated across N_days.
%     Only complete subjects (no NaNs across days) are used for omnibus and
%     pairwise tests for each level.
%   - If samples are independent across days, use Kruskal-Wallis instead.

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end

if ndims(Data_v) ~= 3
    error('Data_v must be N_days x N_levels x num_sample');
end

[N_days, N_levels, num_sample] = size(Data_v);

% Preallocate
pvals = nan(1, N_levels);
effect_info = cell(1, N_levels);    % store extra info (kendallW, meddiff, stats)
Xv_per_level = cell(1, N_levels);   % store complete-pairs matrix per level
omnibus_tbls = cell(1, N_levels);
omnibus_stats = cell(1, N_levels);

for j = 1:N_levels
    % Build subjects x days matrix for level j (paired design)
    X = nan(num_sample, N_days);
    for i = 1:N_days
        X(:, i) = squeeze(Data_v(i, j, :));
    end

    % Keep only rows with no NaN across any day (complete subjects)
    valid_rows = all(~isnan(X), 2);
    Xv = X(valid_rows, :);
    Xv_per_level{j} = Xv;   % even if empty

    % Need at least 2 complete subjects and >=2 days
    if size(Xv, 1) < 2 || N_days < 2
        p = NaN; tbl = []; stats = [];
        pvals(j) = p;
        omnibus_tbls{j} = tbl;
        omnibus_stats{j} = stats;
        effect_info{j} = struct('kendallW', NaN, 'median_diff', NaN);
        continue;
    end

    % Choose omnibus based on number of days
    if N_days == 2
        % Use paired Wilcoxon signed-rank as omnibus (equivalent to pairwise)
        try
            [p_pair, h_pair, stats_sr] = signrank(Xv(:,1), Xv(:,2));
            p = p_pair;
            tbl = []; stats = stats_sr;
            meddiff = median(Xv(:,1) - Xv(:,2), 'omitnan');
        catch ME
            warning('signrank (omnibus) failed for level %d: %s', j, ME.message);
            p = NaN; tbl = []; stats = [];
            meddiff = NaN;
        end
        pvals(j) = p;
        omnibus_tbls{j} = tbl;
        omnibus_stats{j} = stats;
        effect_info{j} = struct('kendallW', NaN, 'median_diff', meddiff);
    else
        % N_days >= 3: Paired, rank-based omnibus across days (Friedman)
        try
            [p_f, tbl_f, stats_f] = friedman(Xv, 1, 'off');
            p = p_f;
            tbl = tbl_f;
            stats = stats_f;
        catch ME
            warning('friedman failed for level %d: %s', j, ME.message);
            p = NaN; tbl = []; stats = [];
        end
        pvals(j) = p;
        omnibus_tbls{j} = tbl;
        omnibus_stats{j} = stats;

        % Effect size: Kendall's W = chi2 / (n * (k-1)) if chi2 present
        W = NaN;
        try
            if ~isempty(tbl) && size(tbl,1) >= 2
                chi2 = tbl{2,5};    % 'Chi-sq' cell from friedman table
                n = size(Xv, 1);    % subjects/blocks
                k = N_days;         % treatments/days
                if isnumeric(chi2) && ~isempty(chi2) && n > 0 && k > 1
                    W = chi2 / (n * (k - 1));
                end
            end
        catch
            W = NaN;
        end
        effect_info{j} = struct('kendallW', W, 'median_diff', NaN);
    end
end

% Benjamini-Hochberg FDR across non-NaN p-values
qvals = nan(size(pvals));
signif_levels = false(size(pvals));
valid_mask = ~isnan(pvals);
sp = pvals(valid_mask);
if any(valid_mask)
    sp = sp(:);  % ensure column
    [sp_sorted, sort_idx_rel] = sort(sp, 'ascend');
    m_valid = numel(sp_sorted);
    ranks = (1:m_valid)';

    % BH raw q and monotone adjustment
    raw_q = sp_sorted .* (m_valid ./ ranks);
    adj_q = raw_q;
    for t = m_valid-1:-1:1
        adj_q(t) = min(adj_q(t), adj_q(t+1));
    end
    adj_q = min(adj_q, 1);

    % map back to original positions among valid_mask
    valid_idx = find(valid_mask);
    sort_idx_global = valid_idx(sort_idx_rel);
    qvals(sort_idx_global) = adj_q;

    % determine which are significant under BH at level alpha
    thresh = (ranks / m_valid) * alpha;
    kf = find(sp_sorted <= thresh, 1, 'last');
    if ~isempty(kf)
        sig_global_sorted = sort_idx_global(1:kf);
        signif_levels(sig_global_sorted) = true;
    end
end

% build output and run pairwise signrank only for levels significant after FDR
res.pvals = pvals;
res.qvals = qvals;
res.signif_levels = signif_levels;
for j = 1:N_levels
    res.level(j).omnibus_p = pvals(j);
    res.level(j).omnibus_stats = omnibus_stats{j};
    res.level(j).omnibus_tbl = omnibus_tbls{j};
    res.level(j).effect_info = effect_info{j};
    res.level(j).is_significant = signif_levels(j);
    res.level(j).pairwise = [];

    if ~res.level(j).is_significant
        continue;
    end

    Xv = Xv_per_level{j};   % complete subjects x days for this level
    n_common = size(Xv, 1);

    if N_days == 2
        % Only one pair to report
        if n_common < 2
            p_pair = NaN; h = 0; stats_sr = struct(); meddiff = NaN;
        else
            try
                [p_pair, h, stats_sr] = signrank(Xv(:,1), Xv(:,2));
                meddiff = median(Xv(:,1) - Xv(:,2), 'omitnan');
            catch ME
                warning('signrank failed for level %d pair 1-2: %s', j, ME.message);
                p_pair = NaN; h = 0; stats_sr = struct(); meddiff = NaN;
            end
        end
        res.level(j).pairwise(1).day1 = 1;
        res.level(j).pairwise(1).day2 = 2;
        res.level(j).pairwise(1).p_signrank = p_pair;
        res.level(j).pairwise(1).h = h;
        res.level(j).pairwise(1).stats = stats_sr;
        res.level(j).pairwise(1).n_common = n_common;
        res.level(j).pairwise(1).median_diff = meddiff;
    else
        % Multiple days: run paired signrank for each day pair using complete-subject set
        pair_counter = 0;
        for di = 1:N_days-1
            xi = Xv(:, di);
            for dk = di+1:N_days
                xk = Xv(:, dk);
                pair_counter = pair_counter + 1;
                if numel(xi) < 2
                    p_pair = NaN; h = 0; stats_sr = struct(); meddiff = NaN;
                else
                    try
                        [p_pair, h, stats_sr] = signrank(xi, xk);
                        meddiff = median(xi - xk, 'omitnan');
                    catch ME
                        warning('signrank failed for level %d pair %d-%d: %s', j, di, dk, ME.message);
                        p_pair = NaN; h = 0; stats_sr = struct(); meddiff = NaN;
                    end
                end
                res.level(j).pairwise(pair_counter).day1 = di;
                res.level(j).pairwise(pair_counter).day2 = dk;
                res.level(j).pairwise(pair_counter).p_signrank = p_pair;
                res.level(j).pairwise(pair_counter).h = h;
                res.level(j).pairwise(pair_counter).stats = stats_sr;
                res.level(j).pairwise(pair_counter).n_common = n_common;
                res.level(j).pairwise(pair_counter).median_diff = meddiff;
            end
        end
    end
end

end