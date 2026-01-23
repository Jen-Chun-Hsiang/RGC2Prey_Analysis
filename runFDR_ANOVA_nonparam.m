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

% Top-level metadata for paper reporting / provenance
res = struct();
res.alpha = alpha;
res.N_days = N_days;
res.N_levels = N_levels;
res.num_sample = num_sample;
res.paired_design = true;
res.note = 'Per-level paired omnibus (signrank if N_days==2, friedman if N_days>=3), BH-FDR across levels, paired signrank post-hoc for FDR-significant levels.';

% Preallocate
pvals = nan(1, N_levels);
effect_info = cell(1, N_levels);    % store extra info (kendallW, meddiff, stats)
Xv_per_level = cell(1, N_levels);   % store complete-pairs matrix per level
omnibus_tbls = cell(1, N_levels);
omnibus_stats = cell(1, N_levels);
omnibus_test = cell(1, N_levels);
descriptives = cell(1, N_levels);

% Paper-friendly per-level summaries (aligned with pvals/qvals order)
n_total_per_level = nan(1, N_levels);
n_complete_per_level = nan(1, N_levels);
n_excluded_per_level = nan(1, N_levels);
df_per_level = nan(1, N_levels);
chi2_per_level = nan(1, N_levels);
kendallW_per_level = nan(1, N_levels);
median_diff_per_level = nan(1, N_levels);
rank_biserial_per_level = nan(1, N_levels);
n_eff_per_level = nan(1, N_levels); % signrank effective N after dropping zeros

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

    n_total = size(X, 1);
    n_complete = size(Xv, 1);
    n_excluded = n_total - n_complete;

    n_total_per_level(j) = n_total;
    n_complete_per_level(j) = n_complete;
    n_excluded_per_level(j) = n_excluded;

    % Per-day descriptives (computed on complete paired rows)
    if isempty(Xv)
        descriptives{j} = struct('n_total', n_total, 'n_complete', n_complete, 'n_excluded', n_excluded, ...
            'median_per_day', nan(1, N_days), 'mean_per_day', nan(1, N_days), 'std_per_day', nan(1, N_days));
    else
        descriptives{j} = struct('n_total', n_total, 'n_complete', n_complete, 'n_excluded', n_excluded, ...
            'median_per_day', median(Xv, 1), 'mean_per_day', mean(Xv, 1), 'std_per_day', std(Xv, 0, 1));
    end

    % Need at least 2 complete subjects and >=2 days
    if size(Xv, 1) < 2 || N_days < 2
        p = NaN; tbl = []; stats = [];
        pvals(j) = p;
        omnibus_tbls{j} = tbl;
        omnibus_stats{j} = stats;
        omnibus_test{j} = '';
        effect_info{j} = struct('kendallW', NaN, 'chi2', NaN, 'df', NaN, ...
            'median_diff', NaN, 'rank_biserial', NaN, 'n_total', n_total, 'n_complete', n_complete, 'n_eff', NaN);
        continue;
    end

    % Choose omnibus based on number of days
    if N_days == 2
        % Use paired Wilcoxon signed-rank as omnibus (equivalent to pairwise)
        try
            [p_pair, ~, stats_sr] = signrank(Xv(:,1), Xv(:,2));
            p = p_pair;
            tbl = []; stats = stats_sr;
            meddiff = median(Xv(:,1) - Xv(:,2), 'omitnan');
        catch ME
            warning('signrank (omnibus) failed for level %d: %s', j, ME.message);
            p = NaN; tbl = []; stats = [];
            meddiff = NaN;
        end
        omnibus_test{j} = 'signrank';
        pvals(j) = p;
        omnibus_tbls{j} = tbl;
        omnibus_stats{j} = stats;

        % Effect size: matched-pairs rank-biserial correlation (from signed-rank W)
        % r_rb = (Wpos - Wneg) / (Wpos + Wneg) = 2*Wpos/totalRank - 1
        r_rb = NaN;
        n_eff = NaN;
        try
            d = Xv(:,1) - Xv(:,2);
            d = d(~isnan(d));
            d = d(d ~= 0); % signrank drops zeros
            n_eff = numel(d);
            if n_eff == 0
                % All differences are zero => no rank information; treat effect as 0
                r_rb = 0;
            elseif ~isempty(stats_sr) && isfield(stats_sr, 'signedrank') && n_eff >= 1
                totalRank = n_eff * (n_eff + 1) / 2;
                r_rb = (2 * double(stats_sr.signedrank) / totalRank) - 1;
            end
        catch
            % keep defaults
        end

        effect_info{j} = struct('kendallW', NaN, 'chi2', NaN, 'df', 1, ...
            'median_diff', meddiff, 'rank_biserial', r_rb, 'n_total', n_total, 'n_complete', n_complete, 'n_eff', n_eff);

        df_per_level(j) = 1;
        chi2_per_level(j) = NaN;
        kendallW_per_level(j) = NaN;
        median_diff_per_level(j) = meddiff;
        rank_biserial_per_level(j) = r_rb;
        n_eff_per_level(j) = n_eff;
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
        omnibus_test{j} = 'friedman';
        pvals(j) = p;
        omnibus_tbls{j} = tbl;
        omnibus_stats{j} = stats;

        % Effect size: Kendall's W = chi2 / (n * (k-1)) if chi2 present
        W = NaN;
        chi2 = NaN;
        df = N_days - 1;
        try
            if ~isempty(tbl) && size(tbl,1) >= 2 && size(tbl,2) >= 3
                hdr = tbl(1,:);
                chi2_col = find(strcmpi(hdr, 'Chi-sq') | strcmpi(hdr, 'Chi-Sq') | strcmpi(hdr, 'Chi^2') | strcmpi(hdr, 'Chi2'), 1, 'first');
                df_col = find(strcmpi(hdr, 'df'), 1, 'first');
                if ~isempty(chi2_col)
                    chi2 = tbl{2, chi2_col};
                end
                if ~isempty(df_col)
                    df = tbl{2, df_col};
                end
                n = size(Xv, 1);    % subjects/blocks
                k = N_days;         % treatments/days
                if isnumeric(chi2) && ~isempty(chi2) && n > 0 && k > 1
                    W = chi2 / (n * (k - 1));
                end
            end
        catch
            W = NaN;
            chi2 = NaN;
        end
        effect_info{j} = struct('kendallW', W, 'chi2', chi2, 'df', df, ...
            'median_diff', NaN, 'rank_biserial', NaN, 'n_total', n_total, 'n_complete', n_complete, 'n_eff', NaN);

        df_per_level(j) = df;
        chi2_per_level(j) = chi2;
        kendallW_per_level(j) = W;
        median_diff_per_level(j) = NaN;
        rank_biserial_per_level(j) = NaN;
        n_eff_per_level(j) = NaN;
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

% Top-level, paper-friendly per-level summaries
res.omnibus_test = omnibus_test;
res.n_total_per_level = n_total_per_level;
res.n_complete_per_level = n_complete_per_level;
res.n_excluded_per_level = n_excluded_per_level;
res.df_per_level = df_per_level;
res.chi2_per_level = chi2_per_level;
res.kendallW_per_level = kendallW_per_level;
res.median_diff_per_level = median_diff_per_level;
res.rank_biserial_per_level = rank_biserial_per_level;
res.n_eff_per_level = n_eff_per_level;

% Unified effect size vector (for convenience in paper tables)
if N_days == 2
    res.effect_size_label = 'rank_biserial';
    res.effect_size_per_level = rank_biserial_per_level;
else
    res.effect_size_label = 'kendallW';
    res.effect_size_per_level = kendallW_per_level;
end
for j = 1:N_levels
    res.level(j).omnibus_p = pvals(j);
    res.level(j).omnibus_test = omnibus_test{j};
    res.level(j).omnibus_stats = omnibus_stats{j};
    res.level(j).omnibus_tbl = omnibus_tbls{j};
    res.level(j).effect_info = effect_info{j};
    res.level(j).descriptives = descriptives{j};
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

        % Effect size: matched-pairs rank-biserial correlation
        r_rb = NaN;
        n_eff = NaN;
        try
            d = Xv(:,1) - Xv(:,2);
            d = d(~isnan(d));
            d = d(d ~= 0);
            n_eff = numel(d);
            if ~isempty(stats_sr) && isfield(stats_sr, 'signedrank') && n_eff >= 1
                totalRank = n_eff * (n_eff + 1) / 2;
                r_rb = (2 * double(stats_sr.signedrank) / totalRank) - 1;
            end
        catch
            % keep defaults
        end
        res.level(j).pairwise(1).n_eff = n_eff;
        res.level(j).pairwise(1).rank_biserial = r_rb;
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

                % Effect size: matched-pairs rank-biserial correlation
                r_rb = NaN;
                n_eff = NaN;
                try
                    d = xi - xk;
                    d = d(~isnan(d));
                    d = d(d ~= 0);
                    n_eff = numel(d);
                    if ~isempty(stats_sr) && isfield(stats_sr, 'signedrank') && n_eff >= 1
                        totalRank = n_eff * (n_eff + 1) / 2;
                        r_rb = (2 * double(stats_sr.signedrank) / totalRank) - 1;
                    end
                catch
                    % keep defaults
                end
                res.level(j).pairwise(pair_counter).n_eff = n_eff;
                res.level(j).pairwise(pair_counter).rank_biserial = r_rb;
            end
        end
    end
end

end