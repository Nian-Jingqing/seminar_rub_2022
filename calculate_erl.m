
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


% Iterate datasets
behavior = []
counter = 0;
for s = 1 : numel(fl)

    % Load data
    EEG = pop_loadset('filename', fl(s).name, 'filepath', PATH_AUTOCLEANED, 'loadmode', 'all');

    % Get id as integer
    id = str2num(fl(s).name(1 : 2));

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

        % Get behavioral measures
        rt =        mean(EEG.trialinfo(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 1 & EEG.trialinfo(:, 2) ~= 0, 7), 1);
        acc =       size(EEG.trialinfo(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 1 & EEG.trialinfo(:, 2) ~= 0, 1), 1);
        incorrect = size(EEG.trialinfo(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 0 & EEG.trialinfo(:, 2) ~= 0, 1), 1);
        omission =  size(EEG.trialinfo(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 2 & EEG.trialinfo(:, 2) ~= 0, 1), 1);

        counter = counter + 1;
        behavior(counter, :) = [id, cnd, rt, acc, incorrect, omission, n_trials, sum(target_left), sum(target_right)];


    end
end

mean_behav = [];
for cnd = 1 : 3

    mean_behav(cnd, :) = mean(behavior(behavior(:, 2) == cnd, :), 1);

end

figure()
labels = {'rt', 'accuracy', 'incorrect', 'omission'}
for p = 3 : 6
    subplot(2, 2, p - 2)
    plot([1, 2, 3], mean_behav(:, p))
    title(labels{p - 2})
end

% n2pc
ipsi_po78_av = squeeze(nanmean(erl_data_ipsi(:, 1, 24, :), 1));
ipsi_po78_a  = squeeze(nanmean(erl_data_ipsi(:, 2, 24, :), 1));
ipsi_po78_v  = squeeze(nanmean(erl_data_ipsi(:, 3, 24, :), 1));

contra_po78_av = squeeze(nanmean(erl_data_contra(:, 1, 24, :), 1));
contra_po78_a  = squeeze(nanmean(erl_data_contra(:, 2, 24, :), 1));
contra_po78_v  = squeeze(nanmean(erl_data_contra(:, 3, 24, :), 1));

erl_po78_av = squeeze(nanmean(erl_data(:, 1, 24, :), 1));
erl_po78_a  = squeeze(nanmean(erl_data(:, 2, 24, :), 1));
erl_po78_v  = squeeze(nanmean(erl_data(:, 3, 24, :), 1));

% n2ac
ipsi_fc34_av = squeeze(nanmean(erl_data_ipsi(:, 1, 2, :), 1));
ipsi_fc34_a  = squeeze(nanmean(erl_data_ipsi(:, 2, 2, :), 1));
ipsi_fc34_v  = squeeze(nanmean(erl_data_ipsi(:, 3, 2, :), 1));

contra_fc34_av = squeeze(nanmean(erl_data_contra(:, 1, 2, :), 1));
contra_fc34_a  = squeeze(nanmean(erl_data_contra(:, 2, 2, :), 1));
contra_fc34_v  = squeeze(nanmean(erl_data_contra(:, 3, 2, :), 1));

erl_fc34_av = squeeze(nanmean(erl_data(:, 1, 2, :), 1));
erl_fc34_a  = squeeze(nanmean(erl_data(:, 2, 2, :), 1));
erl_fc34_v  = squeeze(nanmean(erl_data(:, 3, 2, :), 1));

% Plot contra versus ipsi and ERLs
figure()
subplot(2, 2, 1)
plot(EEG.times, ipsi_po78_av, 'k-', 'LineWidth', 2)
hold on
plot(EEG.times, ipsi_po78_a, 'm-', 'LineWidth', 2)
plot(EEG.times, ipsi_po78_v, 'c-', 'LineWidth', 2)
plot(EEG.times, contra_po78_av, 'k:', 'LineWidth', 2)
plot(EEG.times, contra_po78_a, 'm:', 'LineWidth', 2)
plot(EEG.times, contra_po78_v, 'c:', 'LineWidth', 2)
legend({'av-ipsi', 'a-ipsi', 'v-ipsi', 'av-contra', 'a-contra', 'v-contra'})
grid on
ylim([-3.5, 2.5])
title('PO7/8 - contra vs ipsi')


subplot(2, 2, 2)
plot(EEG.times, ipsi_fc34_av, 'k-', 'LineWidth', 2)
hold on
plot(EEG.times, ipsi_fc34_a, 'm-', 'LineWidth', 2)
plot(EEG.times, ipsi_fc34_v, 'c-', 'LineWidth', 2)
plot(EEG.times, contra_fc34_av, 'k:', 'LineWidth', 2)
plot(EEG.times, contra_fc34_a, 'm:', 'LineWidth', 2)
plot(EEG.times, contra_fc34_v, 'c:', 'LineWidth', 2)
legend({'av-ipsi', 'a-ipsi', 'v-ipsi', 'av-contra', 'a-contra', 'v-contra'})
grid on
ylim([-3.5, 2.5])
title('FC3/4 - contra vs ipsi')

subplot(2, 2, 3)
plot(EEG.times, erl_po78_av, 'k-', 'LineWidth', 2.5)
hold on
plot(EEG.times, erl_po78_a, 'm-', 'LineWidth', 2.5)
plot(EEG.times, erl_po78_v, 'c-', 'LineWidth', 2.5)
legend({'av', 'a', 'v'})
grid on
ylim([-1.8, 1.5])
title('PO7/8 - ERL')

subplot(2, 2, 4)
plot(EEG.times, erl_fc34_av, 'k-', 'LineWidth', 2.5)
hold on
plot(EEG.times, erl_fc34_a, 'm-', 'LineWidth', 2.5)
plot(EEG.times, erl_fc34_v, 'c-', 'LineWidth', 2.5)
legend({'av', 'a', 'v', '270ms'})
grid on
ylim([-1.8, 2.5])
title('FC3/4 - ERL')







aa=bb


% Topo at 270 ms
[~, time_idx] = min(abs(EEG.times - 270));
av270 = squeeze(mean(erl_data(:, 1, :, time_idx), 1));
a270  = squeeze(mean(erl_data(:, 2, :, time_idx), 1));
v270  = squeeze(mean(erl_data(:, 3, :, time_idx), 1));

% Topo at 400 ms
[~, time_idx] = min(abs(EEG.times - 400));
av400 = squeeze(mean(erl_data(:, 1, :, time_idx), 1));
a400  = squeeze(mean(erl_data(:, 2, :, time_idx), 1));
v400  = squeeze(mean(erl_data(:, 3, :, time_idx), 1));

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
