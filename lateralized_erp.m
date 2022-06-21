
% Clear residuals
clear all;

% Path variables
PATH_EEGLAB        = '/home/plkn/eeglab2022.0/';
PATH_AUTOCLEANED   = '/mnt/data_fast/rub_seminar_2022/autocleaned/';

% Get vhdr file list
fl = dir([PATH_AUTOCLEANED, '*.set']);

% Initialize eeglab
addpath(PATH_EEGLAB);
eeglab;

% Load info
EEG = pop_loadset('filename', fl(1).name, 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

% An erl matrix
erl_data = zeros(numel(fl), 3, EEG.nbchan, EEG.pnts);

% Iterate datasets
for s = 1 : numel(fl)

    % Load data
    EEG = pop_loadset('filename', fl(s).name, 'filepath', PATH_AUTOCLEANED, 'loadmode', 'all');

    % Define matching lateralized electrodes
    chans_left   = [1,  2,  6,  7,  8, 11, 12, 16, 19, 20, 24, 25, 29];
    chans_right  = [5,  4, 10,  9, 32, 15, 14, 18, 23, 22, 28, 27, 31];
    chans_center = [3, 13, 17, 21, 26, 30, 33];

    % Iterate conditions
    for cnd = 1 : 3

        % Get indices of left and right targets
        target_left  = EEG.trialinfo(:, 2) == 1 & EEG.trialinfo(:, 1) == cnd;
        target_right = EEG.trialinfo(:, 2) == 2 & EEG.trialinfo(:, 1) == cnd;

        % Get ipsi and contra potentials
        ipsi_l   = squeeze(mean(EEG.data(chans_left,  :, target_left), 3));
        ipsi_r   = squeeze(mean(EEG.data(chans_right, :, target_right), 3));
        contra_l = squeeze(mean(EEG.data(chans_right, :, target_left), 3));
        contra_r = squeeze(mean(EEG.data(chans_left,  :, target_right), 3));

        % Calculate lateralized potentials
        erl_data(s, cnd, chans_left, :)  = ((contra_l - ipsi_l) + (contra_r - ipsi_r)) / 2;
        erl_data(s, cnd, chans_right, :) = ((contra_l - ipsi_l) + (contra_r - ipsi_r)) / 2;

    end
end


% n2pc
po78_av = squeeze(mean(erl_data(:, 1, 24, :), 1));
po78_a  = squeeze(mean(erl_data(:, 2, 24, :), 1));
po78_v  = squeeze(mean(erl_data(:, 3, 24, :), 1));

% n2pc
fc34_av = squeeze(mean(erl_data(:, 1, 2, :), 1));
fc34_a  = squeeze(mean(erl_data(:, 2, 2, :), 1));
fc34_v  = squeeze(mean(erl_data(:, 3, 2, :), 1));

figure()
subplot(2, 2, 1)
plot(EEG.times, po78_av, 'LineWidth', 2.5)
hold on
plot(EEG.times, po78_a, 'LineWidth', 2.5)
plot(EEG.times, po78_v, 'LineWidth', 2.5)
legend({'av', 'a', 'v'})
xline(270)
title('PO7/8')

subplot(2, 2, 2)
plot(EEG.times, fc34_av, 'LineWidth', 2.5)
hold on
plot(EEG.times, fc34_a, 'LineWidth', 2.5)
plot(EEG.times, fc34_v, 'LineWidth', 2.5)
legend({'av', 'a', 'v', '270ms'})
xline(400)
title('FC3/4')

% Tpo at 270 ms
[~, time_idx] = min(abs(EEG.times - 270));
av270 = squeeze(mean(erl_data(:, 1, :, time_idx), 1));
a270  = squeeze(mean(erl_data(:, 2, :, time_idx), 1));
v270  = squeeze(mean(erl_data(:, 3, :, time_idx), 1));

% Tpo at 400 ms
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
