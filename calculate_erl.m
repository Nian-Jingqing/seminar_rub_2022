
% Clear residuals
clear all;

% Path variables
PATH_EEGLAB        = '/home/plkn/eeglab2022.0/';
PATH_AUTOCLEANED   = '/mnt/data_fast/rub_seminar_2022/3_autocleaned/';

% Get vhdr file list
fl = dir([PATH_AUTOCLEANED, '*.set']);

% Initialize eeglab
addpath(PATH_EEGLAB);
eeglab;

% Load info
EEG = pop_loadset('filename', fl(1).name, 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

% An erl matrix
erl_data = zeros(numel(fl), 3, EEG.nbchan, EEG.pnts);

erl_data_ipsi = zeros(numel(fl), 3, EEG.nbchan, EEG.pnts);
erl_data_contra = zeros(numel(fl), 3, EEG.nbchan, EEG.pnts);


% =============================== COMPUTE THE EVENT_RELATED LATERALIZATIONS ========================================================================

% Some variables...
behavior = []
counter = 0;
id_list = [];

% Iterate datasets
for s = 1 : numel(fl)

    % Load data
    EEG = pop_loadset('filename', fl(s).name, 'filepath', PATH_AUTOCLEANED, 'loadmode', 'all');

    % Get id as integer
    id = str2num(fl(s).name(1 : 2));
    id_list(end + 1) = id;

    % Define matching lateralized electrodes
    chans_left   = [1,  2,  6,  7,  8, 11, 12, 16, 19, 20, 24, 25, 29];
    chans_right  = [5,  4, 10,  9, 32, 15, 14, 18, 23, 22, 28, 27, 31];
    chans_center = [3, 13, 17, 21, 26, 30, 33];

    % Iterate conditions
    for cnd = 1 : 3


        % Trialinfo columns:
        %
        % 1: condition from EEG (1=AV, 2=A, 3=V)
        % 2: stimulus side from EEG (1=left, 2=right)
        % 3: condition from eprime-file (1=AV, 2=A, 3=V)
        % 4: correct response from eprime-file (1=left, 2=right, 0=none)
        % 5: observed response from eprime-file (1=left, 2=right, 0=none)
        % 6: accuracy (1=correct, 2=omission, 0=incorrect)
        % 7: RT

        % Get indices of left and right targets. Only correct trials.
        target_left  = EEG.trialinfo(:, 2) == 1 & EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 1;
        target_right = EEG.trialinfo(:, 2) == 2 & EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 1;

        % Get ipsi and contra potentials
        ipsi_l   = squeeze(mean(EEG.data(chans_left,  :, target_left), 3));
        ipsi_r   = squeeze(mean(EEG.data(chans_right, :, target_right), 3));
        contra_l = squeeze(mean(EEG.data(chans_right, :, target_left), 3));
        contra_r = squeeze(mean(EEG.data(chans_left,  :, target_right), 3));

        % Calculate ipsi and contra activations
        erl_data_ipsi(s, cnd, chans_left,  :) = (ipsi_l + ipsi_r) / 2;
        erl_data_ipsi(s, cnd, chans_right, :) = (ipsi_l + ipsi_r) / 2;
        erl_data_contra(s, cnd, chans_left,  :) = (contra_l + contra_r) / 2;
        erl_data_contra(s, cnd, chans_right, :) = (contra_l + contra_r) / 2;

        % Calculate lateralized potentials
        erl_data(s, cnd, chans_left,  :) = ((contra_l - ipsi_l) + (contra_r - ipsi_r)) / 2;
        erl_data(s, cnd, chans_right, :) = ((contra_l - ipsi_l) + (contra_r - ipsi_r)) / 2;

        % Get number of trials 
        n_trials =  size(EEG.trialinfo(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 2) ~= 0, :), 1);
        n_trials_correct =  size(EEG.trialinfo(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 2) ~= 0 & EEG.trialinfo(:, 6) == 1, :), 1);

        % Get behavioral measures
        rt =        nanmean(EEG.trialinfo(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 1 & EEG.trialinfo(:, 2) ~= 0, 7), 1);
        acc =       size(EEG.trialinfo(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 1 & EEG.trialinfo(:, 2) ~= 0, 1), 1) / n_trials;
        incorrect = size(EEG.trialinfo(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 0 & EEG.trialinfo(:, 2) ~= 0, 1), 1) / n_trials;
        omission =  size(EEG.trialinfo(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 2 & EEG.trialinfo(:, 2) ~= 0, 1), 1) / n_trials;

        % Collect behavioral measure in a matrix
        counter = counter + 1;
        behavior(counter, :) = [id, cnd, rt, acc, incorrect, omission, n_trials_correct];

    end
end

% Identify subjects that had an accuracy of at least chance level
av_to_keep = behavior(behavior(:, 4) > 0.5 & behavior(:, 2) == 1, 1); % Auvi above chance
v_to_keep  = behavior(behavior(:, 4) > 0.5 & behavior(:, 2) == 3, 1); % Vi only above chance
id_to_keep = intersect(av_to_keep, v_to_keep);

% Get idx in file-list of those subjects
idx_to_keep = [];
for i = 1 : length(id_to_keep)
    idx_to_keep(i) = find(id_list == id_to_keep(i));
end

% Remove sub-threshold performers from eeeg-data
erl_data_ipsi = erl_data_ipsi(idx_to_keep, :, :, :);
erl_data_contra = erl_data_contra(idx_to_keep, :, :, :);
erl_data = erl_data(idx_to_keep, :, :, :);

% Select posterior electrode patch (P7/8, PO7/8, PO3/4)
idx_chan_v = [19, 24, 25];

% Average across selected electrodes
erl_posterior_ipsi_av   = squeeze(mean(mean(erl_data_ipsi(:, 1, idx_chan_v, :), 1), 3));
erl_posterior_ipsi_v    = squeeze(mean(mean(erl_data_ipsi(:, 3, idx_chan_v, :), 1), 3));
erl_posterior_contra_av = squeeze(mean(mean(erl_data_contra(:, 1, idx_chan_v, :), 1), 3));
erl_posterior_contra_v  = squeeze(mean(mean(erl_data_contra(:, 3, idx_chan_v, :), 1), 3));
erl_posterior_diff_av   = squeeze(mean(mean(erl_data(:, 1, idx_chan_v, :), 1), 3));
erl_posterior_diff_v    = squeeze(mean(mean(erl_data(:, 3, idx_chan_v, :), 1), 3));

% Define time window for topography
time_win = [230, 280];
idx_time_win = EEG.times >= time_win(1) & EEG.times <= time_win(2);

% Get values for plotting topography
[~, time_idx] = min(abs(EEG.times - 400));
topovals_av = squeeze(mean(mean(erl_data(:, 1, :, idx_time_win), 1), 4));
topovals_v  = squeeze(mean(mean(erl_data(:, 3, :, idx_time_win), 1), 4));

% Plot frontal asymmetry auvi versus visu
figure()
subplot(2, 2, 1)
plot(EEG.times, erl_posterior_ipsi_av, 'k:', 'LineWidth', 2)
hold on
plot(EEG.times, erl_posterior_ipsi_v, 'r:', 'LineWidth', 2)
plot(EEG.times, erl_posterior_contra_av, 'k-', 'LineWidth', 2)
plot(EEG.times, erl_posterior_contra_v, 'r-', 'LineWidth', 2)
legend({'auvi-ipsi', 'visu-ipsi', 'auvi-contra', 'visu-contra'})
grid on
ylim([-3.5, 2.5])
title('contra vs ipsi at [P7/8, PO7/8, PO3/4]')

subplot(2, 2, 2)
plot(EEG.times, erl_posterior_diff_av, 'k-', 'LineWidth', 2.5)
hold on
plot(EEG.times, erl_posterior_diff_v, 'r-', 'LineWidth', 2.5)
legend({'auvi', 'visu'})
grid on
ylim([-2.2, 1.5])
title('Lateralization at [P7/8, PO7/8, PO3/4]')
xline(time_win(1))
xline(time_win(2))

% Color limits for topos
clim = [-1.8, 1.8];

subplot(2, 2, 3)
topoplot(topovals_av, EEG.chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
colormap('jet');
caxis(clim);
title(['Auvi from ', num2str(time_win(1)), ' to ', num2str(time_win(2)), ' ms'])
colorbar()

subplot(2, 2, 4)
topoplot(topovals_v, EEG.chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
colormap('jet');
caxis(clim);
title(['Visu from ', num2str(time_win(1)), ' to ', num2str(time_win(2)), ' ms'])
colorbar()

% Set some parameters
pval_cluster = 0.25;
n_perms = 1000;
n_subjects = size(erl_data, 1);
pval_voxel = 0.01;

% Select a time window for analysis
idx_time_win_analysis = EEG.times >= 0 & EEG.times <= 800;

% Data for statistical analysis as dubects x time
analysis_times = EEG.times(idx_time_win_analysis);
data1 = squeeze(mean(erl_data(:, 1, idx_chan_v, idx_time_win_analysis), 3)); % Auvi
data2 = squeeze(mean(erl_data(:, 3, idx_chan_v, idx_time_win_analysis), 3)); % Visu

% Init matrices
permuted_t = zeros(n_perms, size(data1, 2));
max_tsum = zeros(n_perms, 2);
max_nvox = zeros(n_perms, 1);
desmat = [zeros(n_subjects, 1), ones(n_subjects, 1)];

% Iterate permutations
for perm = 1 : n_perms

    % Permute
    toflip = randsample(n_subjects, floor(n_subjects / 2));
    d1_perm = data1;
    d1_perm(toflip, :) = data2(toflip, :);
    d2_perm = data2;
    d2_perm(toflip, :) = data1(toflip, :);

    % Calculate and save t values
    tnum = squeeze(mean(d1_perm - d2_perm, 1));
    tdenum = squeeze(std(d1_perm - d2_perm, 0, 1)) / sqrt(n_subjects);
    fake_t = tnum ./ tdenum;
    permuted_t(perm, :) = fake_t;

    % Threshold t values
    fake_t(abs(fake_t) < tinv(1 - pval_voxel, n_subjects) - 1) = 0;
    fake_t = logical(fake_t);

    % Identify clusters
    [clust_labels, n_clusts] = bwlabel(fake_t);

    % Determine min and mux sum of t in clusters
    sum_t = [];
    sum_vox = [];
    for clu = 1 : n_clusts
        sum_t(end + 1) = sum(fake_t(clust_labels == clu));
        sum_vox(end + 1) = sum(clust_labels == clu);
    end

    % Collect min and max cluster statistics
    max_tsum(perm, 1) = min([0, sum_t]);
    max_tsum(perm, 2) = max([0, sum_t]);
    max_nvox(perm) = max([0, sum_vox]);

end

% T-test the real thing
tnum = squeeze(mean(data1 - data2, 1));
tdenum = squeeze(std(data1 - data2, 0, 1)) / sqrt(size(data1, 1));
tmat = tnum ./ tdenum;

% Save for later before thresholding
tvals = tmat;

% Threshold t values
tmat(abs(tmat) < tinv(1 - pval_voxel, size(data1, 1) - 1)) = 0;
tmat = logical(tmat);

% Identify clusters
[clust_labels, n_clusts] = bwlabel(tmat);

% Determine min and mux sum of t in clusters
sum_t = [];
sum_vox = [];
for clu = 1 : n_clusts
    sum_t(end + 1) = sum(fake_t(clust_labels == clu));
    sum_vox(end + 1) = sum(clust_labels == clu);
end

% Determine upper and lower thresholds
clust_thresh_lower = prctile(max_tsum(:, 1), pval_cluster * 100);
clust_thresh_upper = prctile(max_tsum(:, 2), 100 - pval_cluster * 100);
clust_thresh_nvox  = prctile(max_nvox, 100 - pval_cluster * 100);

% Determine cluster to keep
clust2keep = find(sum_t <= clust_thresh_lower | sum_t >= clust_thresh_upper);
%clust2keep = find(sum_vox >= clust_thresh_nvox);

% Build cluster vector
clust_vector = zeros(size(tmat));
for clu = 1 : length(clust2keep)
    clust_vector(clust_labels == clust2keep(clu)) = 1;
end

% Set the flag of significance
sig_flag = logical(sum(clust_vector(:)));

% Calculate effect sizes
x = tvals.^2 ./ (tvals.^2 + (n_subjects - 1));
apes = x - (1 - x) .* (1 / (n_subjects - 1));

% Calculate averages
mean_data1 = squeeze(mean(data1, 1));
mean_data2 = squeeze(mean(data2, 1));


figure()
plot(analysis_times, mean_data1)
hold on
plot(analysis_times, mean_data2)
plot(analysis_times, tvals)
plot(analysis_times, apes)
plot(analysis_times, tmat)
plot(analysis_times, clust_vector)
legend('erl-av', 'erl-v', 't-vals', 'eta', 'cluster-thresh', 'cluster')
















aa=bb











% Get average across subjects for both visual conditions
gga_v = squeeze(mean(mean(erl_data(:, [1, 3], 24, :), 2), 1));

% Get average across subjects for both auditory conditions
gga_a = squeeze(mean(mean(erl_data(:, [1, 2], 2, :), 2), 1));

% Get idx of n2pc time window
[~, idx_n2pc] = min(gga_v);
time_n2pc = EEG.times(idx_n2pc);
width_n2pc = 50;
idx_win_n2pc = EEG.times >= time_n2pc - (width_n2pc / 2) & EEG.times <= time_n2pc + (width_n2pc / 2);

% Get idx of n2ac time window
[~, idx_n2ac] = min(gga_a);
time_n2ac = EEG.times(idx_n2ac);
width_n2ac = 50;
idx_win_n2ac = EEG.times >= time_n2ac - (width_n2ac / 2) & EEG.times <= time_n2ac + (width_n2ac / 2);

% Plot time window for n2pc analysis
figure()
subplot(1, 2, 1)
plot(EEG.times, gga_v)
xline(time_n2pc - (width_n2pc / 2))
xline(time_n2pc + (width_n2pc / 2))
title('n2pc time window in averaged data at po7/8')

subplot(1, 2, 2)
plot(EEG.times, gga_a)
xline(time_n2ac - (width_n2ac / 2))
xline(time_n2ac + (width_n2ac / 2))
title('n2ac time window in averaged data at fc3/4')



% n2ac
idx_chan_a = [12, 7, 6];
ipsi_fc34_av   = squeeze(mean(mean(erl_data_ipsi(:, 1, idx_chan_a, :), 1), 3));
ipsi_fc34_a    = squeeze(mean(mean(erl_data_ipsi(:, 2, idx_chan_a, :), 1), 3));
ipsi_fc34_v    = squeeze(mean(mean(erl_data_ipsi(:, 3, idx_chan_a, :), 1), 3));
contra_fc34_av = squeeze(mean(mean(erl_data_contra(:, 1, idx_chan_a, :), 1), 3));
contra_fc34_a  = squeeze(mean(mean(erl_data_contra(:, 2, idx_chan_a, :), 1), 3));
contra_fc34_v  = squeeze(mean(mean(erl_data_contra(:, 3, idx_chan_a, :), 1), 3));
erl_fc34_av    = squeeze(mean(mean(erl_data(:, 1, idx_chan_a, :), 1), 3));
erl_fc34_a     = squeeze(mean(mean(erl_data(:, 2, idx_chan_a, :), 1), 3));
erl_fc34_v     = squeeze(mean(mean(erl_data(:, 3, idx_chan_a, :), 1), 3));



subplot(2, 2, 3)
plot(EEG.times, erl_po78_av, 'k-', 'LineWidth', 2.5)
hold on
plot(EEG.times, erl_po78_a, 'm-', 'LineWidth', 2.5)
plot(EEG.times, erl_po78_v, 'c-', 'LineWidth', 2.5)
legend({'av', 'a', 'v'})
grid on
ylim([-2.2, 1.5])
title('PO7/8 - ERL')

subplot(2, 2, 4)
plot(EEG.times, erl_fc34_av, 'k-', 'LineWidth', 2.5)
hold on
plot(EEG.times, erl_fc34_a, 'm-', 'LineWidth', 2.5)
plot(EEG.times, erl_fc34_v, 'c-', 'LineWidth', 2.5)
legend({'av', 'a', 'v', '270ms'})
grid on
ylim([-2.2, 2.5])
title('FC3/4 - ERL')

% Topo n2pc
av270 = squeeze(mean(mean(erl_data(:, 1, :, idx_win_n2pc), 1), 4));
a270  = squeeze(mean(mean(erl_data(:, 2, :, idx_win_n2pc), 1), 4));
v270  = squeeze(mean(mean(erl_data(:, 3, :, idx_win_n2pc), 1), 4));

% Topo at 400 ms
[~, time_idx] = min(abs(EEG.times - 400));
av400 = squeeze(mean(mean(erl_data(:, 1, :, idx_win_n2ac), 1), 4));
a400  = squeeze(mean(mean(erl_data(:, 2, :, idx_win_n2ac), 1), 4));
v400  = squeeze(mean(mean(erl_data(:, 3, :, idx_win_n2ac), 1), 4));

figure()


clim = [-1.8, 1.8];
subplot(2, 3, 1)
topoplot(av270, EEG.chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
colormap('jet');
caxis(clim);
title('av - 270ms')
subplot(2, 3, 2)
topoplot(a270, EEG.chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
colormap('jet');
caxis(clim);
title('a - 270ms')
subplot(2, 3, 3)
topoplot(v270, EEG.chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
colormap('jet');
caxis(clim);
title('v - 270ms')


clim = [-1, 1];
subplot(2, 3, 4)
topoplot(av400, EEG.chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
colormap('jet');
caxis(clim);
title('av - 400ms')
subplot(2, 3, 5)
topoplot(a400, EEG.chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
colormap('jet');
caxis(clim);
title('a - 400ms')
subplot(2, 3, 6)
topoplot(v400, EEG.chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'on');
colormap('jet');
caxis(clim);
title('v - 400ms')



