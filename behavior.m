clear all;

% Path variables
PATH_EEGLAB        = '/home/plkn/eeglab2022.0/';
PATH_RAW_DATA      = '/mnt/data_fast/rub_seminar_2022/0_raw_data/';
PATH_BEHAVIOR      = '/mnt/data_fast/rub_seminar_2022/0_eprime_export/';
PATH_CODED_DATA    = '/mnt/data_fast/rub_seminar_2022/1_coded_data/';


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
        id_str = tmp{1};
        id = str2num(tmp{1});

        % Open trial table
        trial_table = readtable([PATH_BEHAVIOR, 'S', num2str(id_str), '_eprime.xlsx'], 'ReadVariableNames', true);

        % Condition column
        conds = trial_table.Condition;

        % Correct response columns
        cresp_fix = trial_table.Fixation_CRESP_SubTrial_;
        cresp_a   = trial_table.TrialDisplayAonly_CRESP;
        cresp_v   = trial_table.TrialDisplay_CRESP;
        cresp_av  = trial_table.TrialDisplayAuvi_CRESP;
        
        % Response columns
        resp_fix = trial_table.Fixation_RESP_SubTrial_;
        resp_a   = trial_table.TrialDisplayAonly_RESP;
        resp_v   = trial_table.TrialDisplay_RESP;
        resp_av  = trial_table.TrialDisplayAuvi_RESP;

        % RT columns
        rt_fix = trial_table.Fixation_RT_SubTrial_;
        rt_a   = trial_table.TrialDisplayAonly_RT;
        rt_v   = trial_table.TrialDisplay_RT;
        rt_av  = trial_table.TrialDisplayAuvi_RT;

        % Iterate non practice block trials
        behav = [];
        counter = 0;
        for t = 25 : size(trial_table, 1)
    
            % Get condition
            if strcmp(conds{t}, 'Auvi')
                cnd = 1;
            elseif strcmp(conds{t}, 'Aonly')
                cnd = 2;
            elseif strcmp(conds{t}, 'Vonly')
                cnd = 3;
            end
    
            % Get correct response
            cresp = 0;
            if iscell(cresp_fix)
                if strcmp(cresp_fix{t}, 'm')
                    cresp = 2;
                elseif strcmp(cresp_fix{t}, 'n')
                    cresp = 1;
                end
            end
            if iscell(cresp_a)
                if strcmp(cresp_a{t}, 'm')
                    cresp = 2;
                elseif strcmp(cresp_a{t}, 'n')
                    cresp = 1;
                end
    
            end
            if iscell(cresp_v)
                if strcmp(cresp_v{t}, 'm')
                    cresp = 2;
                elseif strcmp(cresp_v{t}, 'n')
                    cresp = 1;
                end
    
            end
            if iscell(cresp_av)
                if strcmp(cresp_av{t}, 'm')
                    cresp = 2;
                elseif strcmp(cresp_av{t}, 'n')
                    cresp = 1;
                end
            end
    
            % Get response
            resp = 0;
            if iscell(resp_fix)
                if strcmp(resp_fix{t}, 'm')
                    resp = 2;
                elseif strcmp(resp_fix{t}, 'n')
                    resp = 1;
                end
            end
            if iscell(resp_a)
                if strcmp(resp_a{t}, 'm')
                    resp = 2;
                elseif strcmp(resp_a{t}, 'n')
                    resp = 1;
                end
    
            end
            if iscell(resp_v)
                if strcmp(resp_v{t}, 'm')
                    resp = 2;
                elseif strcmp(resp_v{t}, 'n')
                    resp = 1;
                end
    
            end
            if iscell(resp_av)
                if strcmp(resp_av{t}, 'm')
                    resp = 2;
                elseif strcmp(resp_av{t}, 'n')
                    resp = 1;
                end
            end
    
            % Get RT
            rt = NaN;
            if ~isnan(rt_fix(t)) & rt_fix(t) ~= 0
                rt = rt_fix(t) + 500;
    
            elseif ~isnan(rt_a(t)) & rt_a(t) ~= 0
                rt = rt_a(t);
    
            elseif ~isnan(rt_v(t)) & rt_v(t) ~= 0
                rt = rt_v(t);
    
            elseif ~isnan(rt_av(t)) & rt_av(t) ~= 0
                rt = rt_av(t);
            end
    
            % Get accuracy
            if cresp == resp
                acc = 1;
            elseif resp == 0
                acc = 2;
            else
                acc = 0;
            end
    
            % Fill matrix
            counter = counter + 1;
            behav(counter, :) = [cnd, cresp, resp, acc, rt];
    
        end

        % Load coded data
        EEG = pop_loadset('filename', fl(s).name, 'filepath', PATH_CODED_DATA, 'loadmode', 'info');

        for t = 1 : size(behav, 1)
            ['trial: ', num2str(t), ' - ', num2str(behav(t, 3)), ' ', resp_fix{t}, ' ', resp_a{t}, ' ', resp_v{t}, ' ', resp_av(t)]
        end



        aa=bb;

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