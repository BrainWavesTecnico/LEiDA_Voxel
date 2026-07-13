function stats = bootstrap_within_permutation_ttest2(data, design, niter, nboot, pthr, htest)
% BOOTSTRAP_WITHIN_PERMUTATION_TTEST2
% Permutation test with within-permutation bootstrap variance estimation.
%
% Tests whether the means of two independent groups differ significantly,
% using a permutation null distribution where variance is estimated via
% bootstrap resampling within each permuted sample.
%
% INPUT:
%   data    - 1 x N vector (single feature, N subjects)
%   design  - 1 x N vector of group labels (1 or 2)
%   niter   - number of permutations (e.g. 1000)
%   nboot   - bootstrap samples per permutation (e.g. 20)
%   pthr    - significance threshold (not used internally, kept for compatibility)
%   htest   - 'ttest' (equal variances) or 'welchtest' (unequal variances)
%
% OUTPUT:
%   stats.tvals       - observed test statistic
%   stats.diffs       - observed difference in group means (group1 - group2)
%   stats.pvals       - [p_right, p_left] one-sided p-values
%   stats.pvals_2sided- two-sided p-value
%   stats.eff         - Hedge's g effect size
%
% Authors: Joana Cabral, Miguel Farinha
% Revised: 2025

%% Input validation
N = length(data);
if length(design) ~= N
    error('data and design must have the same length.');
end

g1 = find(design == 1);
g2 = find(design == 2);

if isempty(g1) || isempty(g2)
    error('Both groups must be non-empty. Group 1: %d subjects, Group 2: %d subjects.', length(g1), length(g2));
end
if ~ismember(htest, {'ttest', 'welchtest'})
    error('htest must be ''ttest'' or ''welchtest''.');
end

n1 = length(g1);
n2 = length(g2);

%% Observed test statistic
obs_stat = compute_stat(data(g1), data(g2), htest);
stats.tvals = obs_stat;
stats.diffs = mean(data(g1)) - mean(data(g2));

%% Permutation null distribution
null_dist = zeros(niter, 1);

for iter = 1:niter
    perm     = randperm(N);
    perm_g1  = data(perm(1:n1));
    perm_g2  = data(perm(n1+1:end));

    % Bootstrap variance estimation within permuted groups
    var_g1 = bootstrap_variance(perm_g1, nboot);
    var_g2 = bootstrap_variance(perm_g2, nboot);

    diff_perm = mean(perm_g1) - mean(perm_g2);
    se = compute_se(var_g1, var_g2, n1, n2, htest);

    if se > 0
        null_dist(iter) = diff_perm / se;
    end
end

%% P-values via kernel density estimate of null distribution
NCDF = max(200, round(200 * niter / 5000));
[fi, xi] = ksdensity(null_dist, 'function', 'cdf', 'npoints', NCDF);

xi_ext = [atanh(-1 + eps), xi, atanh(1 - eps)];
fi_ext = [0, fi, 1];

pval_left  = interp1(xi_ext, fi_ext, obs_stat, 'linear', 'extrap');
pval_left  = max(0, min(1, pval_left));
pval_right = 1 - pval_left;

stats.pvals        = [pval_right, pval_left];
stats.pvals_2sided = 2 * min(pval_right, pval_left);

%% Hedge's g effect size
s_pooled = sqrt(((n1 - 1) * var(data(g1)) + (n2 - 1) * var(data(g2))) / (n1 + n2 - 2));
if s_pooled > 0
    stats.eff = abs(mean(data(g1)) - mean(data(g2))) / s_pooled;
else
    stats.eff = 0;
end

end

%% Helper: compute test statistic from two data vectors
function t = compute_stat(x, y, htest)
    nx = length(x);
    ny = length(y);
    diff = mean(x) - mean(y);
    switch htest
        case 'ttest'
            sp = sqrt(((nx-1)*var(x) + (ny-1)*var(y)) / (nx+ny-2));
            se = sp * sqrt(1/nx + 1/ny);
        case 'welchtest'
            se = sqrt(var(x)/nx + var(y)/ny);
    end
    if se > 0
        t = diff / se;
    else
        t = 0;
    end
end

%% Helper: bootstrap variance estimate for a group
function v = bootstrap_variance(x, nboot)
    if nboot == 0
        v = var(x);   % direct estimate, fine for large samples
        return
    end
    n = length(x);
    boot_vars = zeros(nboot, 1);
    for b = 1:nboot
        idx = randi(n, n, 1);
        boot_vars(b) = var(x(idx));
    end
    v = mean(boot_vars);
end

%% Helper: standard error from variances
function se = compute_se(v1, v2, n1, n2, htest)
    switch htest
        case 'ttest'
            sp = sqrt(((n1-1)*v1 + (n2-1)*v2) / (n1+n2-2));
            se = sp * sqrt(1/n1 + 1/n2);
        case 'welchtest'
            se = sqrt(v1/n1 + v2/n2);
    end
end

