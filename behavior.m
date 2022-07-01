clear all;

% Path variables
PATH_EEGLAB        = '/home/plkn/eeglab2022.0/';
PATH_RAW_DATA      = '/mnt/data2/rub_seminar_2022/0_raw_data/';
PATH_BEHAVIOR      = '/mnt/data2/rub_seminar_2022/0_eprime_export/';
PATH_CODED_DATA    = '/mnt/data2/rub_seminar_2022/1_coded_data/';


% Get vhdr file list
fl = dir([PATH_CODED_DATA, '*.set']);

% Initialize eeglab
addpath(PATH_EEGLAB);
eeglab;

% Iterate datasets
b = [];
c = [];
counter = 0;
for s = 1 : numel(fl)

        % Get id
        tmp = regexp(fl(s).name,'\d*','Match');
        id = str2num(tmp{1});

        % Load coded data
        EEG = pop_loadset('filename', fl(s).name, 'filepath', PATH_CODED_DATA, 'loadmode', 'info');
        aa=bb
        % Double check info
        for t = 1 : size(EEG.trialinfo, 1)

            % Check if eprime and eeg-marker info match
            if (EEG.trialinfo(t, 1) ~= EEG.trialinfo(t, 3)) | (EEG.trialinfo(t, 2) ~= EEG.trialinfo(t, 4))
                error('eeg-marker and eprime info mismatch!!!');
            end

            % Check accuracy coding
            if (EEG.trialinfo(t, 4) == EEG.trialinfo(t, 5)) & (EEG.trialinfo(t, 6) ~= 1)
                error('accuracy should be 1 here...');
            end
            if (EEG.trialinfo(t, 4) ~= EEG.trialinfo(t, 5)) & (EEG.trialinfo(t, 6) == 1)
                error('accuracy should not be 1 here...');
            end
            if (EEG.trialinfo(t, 4) == 0) & (EEG.trialinfo(t, 6) == 2)
                error('This must not be an omission error...');
            end

        end

        % Loop conditions
        for cnd = 1 : 3

            % Some idx
            n_trials_targ0 = sum(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 4) == 0);
            n_trials_targ1 = sum(EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 4) ~= 0);
            co_trial_targ0_idx  = EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 1 & EEG.trialinfo(:, 4) == 0;
            in_trial_targ0_idx  = EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 0 & EEG.trialinfo(:, 4) == 0;
            om_trial_targ0_idx  = EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 2 & EEG.trialinfo(:, 4) == 0;
            co_trial_targ1_idx  = EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 1 & EEG.trialinfo(:, 4) ~= 0;
            in_trial_targ1_idx  = EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 0 & EEG.trialinfo(:, 4) ~= 0;
            om_trial_targ1_idx  = EEG.trialinfo(:, 1) == cnd & EEG.trialinfo(:, 6) == 2 & EEG.trialinfo(:, 4) ~= 0;

            % Get vals
            acc1_targ0 = sum(co_trial_targ0_idx);
            acc0_targ0 = sum(in_trial_targ0_idx);
            acc2_targ0 = sum(om_trial_targ0_idx);
            acc1_targ1 = sum(co_trial_targ1_idx);
            acc0_targ1 = sum(in_trial_targ1_idx);
            acc2_targ1 = sum(om_trial_targ1_idx);
            acc1_targ0_perc = sum(co_trial_targ0_idx) / n_trials_targ0;
            acc0_targ0_perc = sum(in_trial_targ0_idx) / n_trials_targ0;
            acc2_targ0_perc = sum(om_trial_targ0_idx) / n_trials_targ0;
            acc1_targ1_perc = sum(co_trial_targ1_idx) / n_trials_targ1;
            acc0_targ1_perc = sum(in_trial_targ1_idx) / n_trials_targ1;
            acc2_targ1_perc = sum(om_trial_targ1_idx) / n_trials_targ1;
            rt = mean(EEG.trialinfo(co_trial_targ1_idx, 7));

            % Add to table
            counter = counter + 1;
            b(counter, :) = [id,...
                             cnd,...
                             n_trials_targ0,...
                             n_trials_targ1,...
                             acc1_targ0,...
                             acc0_targ0,...
                             acc2_targ0,...
                             acc1_targ1,...
                             acc0_targ1,...
                             acc2_targ1,...
                             acc1_targ0_perc,...
                             acc0_targ0_perc,...
                             acc2_targ0_perc,...
                             acc1_targ1_perc,...
                             acc0_targ1_perc,...
                             acc2_targ1_perc,...
                             rt,...
                            ];

            c(s, 1) = id;
            c(s, cnd + 1) = acc1_targ1_perc;
                                      
        end

        % Save matrix
        writematrix(b, [PATH_BEHAVIOR, 'auvi_behavior.csv']);

end