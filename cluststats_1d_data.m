function [sig_flag, mean_data1, mean_data2, clust_vector, apes] = cluststats__1d_data(data1, data2, varargin)

    %
    % Function performs a cluster based permutation test on 1d data.
    %
    % Inputs:
    % ------------------
    %
    % data1 & data2: 1d data as dimord (unit, dimension1).
    %
    %
    % Optional inputs:
    % ------------------
    %
    % tail:                      Can be 0 (default), 1 or -1. Corresponding to two-sided, left and right testing.
    % pval_voxel:                Thresholding p value (default 0.01).
    % pval_cluster:              Significance level for inference statistic (default 0.05).
    % n_perms:                   Number of permutations for creating the null-hypothesis distribution (default 1000).
    %
    % Returns:
    % ------------------
    %
    % sig_flag:                  The flag of significance. 1 if there are significant differences, 0 if not.          
    % mean_data1 & mean_data2:   The datasets averaged across unit dimension.
    % clust_vector:              Vector of identified clusters.
    % apes:                      The effect sizes (Mordkoff, 2019).
    %
    %
    % Stefan Arnau, June 2020
    % Email: arnau@ifado.de
    %

    % Set defaults
    default_tail         = 0;
    default_pval_voxel   = 0.01;
    default_pval_cluster = 0.05;
    default_n_perms      = 1000;
    
    % Init input parser
    p = inputParser;

    % Parse inputs and set defaults
    p.FunctionName  = mfilename;
    p.CaseSensitive = false;
    p.addRequired('data1', @isnumeric);
    p.addRequired('data2', @isnumeric);
    p.addParamValue('tail', default_tail, @isnumeric);
    p.addParamValue('pval_voxel', default_pval_voxel, @isnumeric);
    p.addParamValue('pval_cluster', default_pval_cluster, @isnumeric);
    p.addParamValue('n_perms', default_n_perms, @isnumeric);
    parse(p, data1, data2, varargin{:});

    % Adjust cluster-p-value for 2-sided testing
    if p.Results.tail
        cpv = p.Results.pval_cluster;
    else
        cpv = p.Results.pval_cluster / 2;
    end

    % Catch some of the possible errors
	if size(p.Results.data2) ~= size(p.Results.data1)
		error('Datasets must be of equal shape... :)');
		return;
	end
    
    % Init matrices
    permuted_t = zeros(p.Results.n_perms, size(p.Results.data1, 2));
    max_clust = zeros(p.Results.n_perms, 2);
    desmat = [zeros(size(p.Results.data1, 1), 1), ones(size(p.Results.data1, 1), 1)];
    for perm = 1 : p.Results.n_perms

        % Talk
        fprintf('\nPermutation %i/%i\n', perm, p.Results.n_perms);

        % Permute
        toflip = find(round(rand(size(p.Results.data1, 1), 1)));
        d1_perm = p.Results.data1;
        d1_perm(toflip, :) = p.Results.data2(toflip, :);
        d2_perm = p.Results.data2;
        d2_perm(toflip, :) = p.Results.data1(toflip, :);

        % Calculate and save t values
        tnum = squeeze(mean(d1_perm - d2_perm, 1));
        tdenum = squeeze(std(d1_perm - d2_perm, 0, 1)) / sqrt(size(p.Results.data1, 1));
        fake_t = tnum ./ tdenum;
        permuted_t(perm, :) = fake_t;

        % Threshold t values
        fake_t(abs(fake_t) < tinv(1 - p.Results.pval_voxel, size(p.Results.data1, 1) - 1)) = 0;

        % Identify clusters
        [clust_labels, n_clusts] = bwlabel(fake_t);

        % Determine min and mux sum of t in clusters
        sum_t = [];
        for clu = 1 : n_clusts
            sum_t(end + 1) = sum(fake_t(clust_labels == clu));
        end
        max_clust(perm, 1) = min([0, sum_t]);
        max_clust(perm, 2) = max([0, sum_t]);      
    end

    % T-test the real thing
    tnum = squeeze(mean(p.Results.data1 - p.Results.data2, 1));
    tdenum = squeeze(std(p.Results.data1 - p.Results.data2, 0, 1)) / sqrt(size(p.Results.data1, 1));
    tmat = tnum ./ tdenum;
    tvals = tmat;
    tmat(abs(tmat) < tinv(1 - p.Results.pval_voxel, size(p.Results.data1, 1) - 1)) = 0;  

    % Identify clusters
    [clust_labels, n_clusts] = bwlabel(fake_t);

    % Determine min and mux sum of t in clusters
    sum_t = [];
    for clu = 1 : n_clusts
        sum_t(end + 1) = sum(fake_t(clust_labels == clu));
    end

    % Determine upper and lower thresholds
    clust_thresh_lower = prctile(max_clust(:, 1), cpv * 100);
    clust_thresh_upper = prctile(max_clust(:, 2), 100 - cpv * 100);

    % Determine cluster to keep
    clust2keep = find(sum_t <= clust_thresh_lower | sum_t >= clust_thresh_upper);

    % Build cluster vector
    clust_vector = zeros(size(tmat));
    for clu = 1 : length(clust2keep)
        clust_vector(clust_labels == clust2keep(clu)) = 1;
    end

    % Set the flag of significance
    sig_flag = logical(sum(clust_vector(:)));

    % Calculate effect sizes
    x = tvals.^2 ./ (tvals.^2 + (size(p.Results.data1, 1) - 1));
    apes = x - (1 - x) .* (1 / (size(p.Results.data1, 1) - 1));

    % Calculate averages
    mean_data1 = squeeze(mean(p.Results.data1, 1));
    mean_data2 = squeeze(mean(p.Results.data2, 1));

end