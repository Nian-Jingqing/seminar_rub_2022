

% Path variables
PATH_EEGLAB        = '/home/plkn/eeglab2022.0/';
PATH_RAW_DATA      = '/mnt/data_fast/rub_seminar_2022/raw_data/';
PATH_ICSET         = '/mnt/data_fast/rub_seminar_2022/ic_sets/';
PATH_BEHAVIOR      = '/mnt/data_fast/rub_seminar_2022/eprime_export/';
PATH_AUTOCLEANED   = '/mnt/data_fast/rub_seminar_2022/autocleaned/';

% Get vhdr file list
fl = dir([PATH_RAW_DATA, '*.vhdr']);

% Initialize eeglab
addpath(PATH_EEGLAB);
eeglab;

% Get channel location file
channel_location_file = which('dipplot.m');
channel_location_file = channel_location_file(1 : end - length('dipplot.m'));
channel_location_file = [channel_location_file, 'standard_BESA/standard-10-5-cap385.elp'];

% Iterate datasets
for s = 1 : numel(fl)

        % Get id
        tmp = regexp(fl(s).name,'\d*','Match');
        id = tmp{1};

        % Open trial table
        trial_table = readtable([PATH_BEHAVIOR, 'S', id, '_eprime.xlsx'], 'ReadVariableNames', true);

        % Get behavior
        resp_fix = trial_table.Fixation_RESP_SubTrial_;
        resp_a   = trial_table.TrialDisplayAonly_RESP;
        resp_v   = trial_table.TrialDisplay_RESP;
        resp_av  = trial_table.TrialDisplayAuvi_RESP;

        rt_fix = trial_table.Fixation_RT_SubTrial_;
        rt_a   = trial_table.TrialDisplayAonly_RT;
        rt_v   = trial_table.TrialDisplay_RT;
        rt_av  = trial_table.TrialDisplayAuvi_RT;

        counter = 0;
        behavior = []
        for t = 25 : size(trial_table, 1)
            
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

            counter = counter + 1;
            behavior(counter, :) = [resp, rt];

        end
        
        % Load dataset
        EEG = pop_loadbv(PATH_RAW_DATA, fl(s).name, [], []);

        % A matrix for the trials
        EEG.trialinfo = [];

        % Iterate events
        counter = 0;
        for e = 1 : length(EEG.event)

            % If stimulus event
            if strcmp(EEG.event(e).type(1 : 2), 'S ')

                % Increase counter
                counter = counter + 1;

                % Get condition
                condition_number = str2num(EEG.event(e).type(3));
                if ismember(condition_number, [1, 2, 3])
                    av_cond = 1; % AV
                elseif ismember(condition_number, [4, 5, 6])
                    av_cond = 2; % A
                elseif ismember(condition_number, [7, 8, 9])
                    av_cond = 3; % V
                end

                % Get position
                condition_number = str2num(EEG.event(e).type(3));
                if ismember(condition_number, [1, 4, 7])
                    stimpos = 1;
                elseif ismember(condition_number, [2, 5, 8])
                    stimpos = 2;
                elseif ismember(condition_number, [3, 6, 9])
                    stimpos = 0;
                end

                % Get stimulus identity
                stim_id = str2num(EEG.event(e).type(4));

                % Fill trial table
                EEG.trialinfo(counter, :) = [av_cond, stimpos];

                % Rename trigger
                EEG.event(e).type = 'stimulus';

            end
        end

        % Check number of trials
        if size(EEG.trialinfo, 1) == 648 & size(behavior, 1) == 648
            for t = 1 : size(EEG.trialinfo, 1)

                % Add response
                EEG.trialinfo(t, 3) = behavior(t, 1);

                % Add accuracy
                if EEG.trialinfo(t, 2) == EEG.trialinfo(t, 3)
                    EEG.trialinfo(t, 4) = 1;
                elseif EEG.trialinfo(t, 3) == 0
                    EEG.trialinfo(t, 4) = 2;
                else
                    EEG.trialinfo(t, 4) = 0;
                end

                % Add rt
                EEG.trialinfo(t, 5) = behavior(t, 2);

            end
        else
            error('number of trials seems to be not correct...')
        end

        % Add channel locations
        EEG = pop_chanedit(EEG, 'lookup', channel_location_file);

        % Save original channel locations (for later interpolation)
        EEG.chanlocs_original = EEG.chanlocs;

        % Filter
        EEG = pop_basicfilter(EEG, [1 : EEG.nbchan], 'Cutoff', [0.5, 20], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 4, 'RemoveDC', 'on', 'Boundary', 'boundary');
            
        % Bad channel detection
        [EEG, EEG.chans_rejected] = pop_rejchan(EEG, 'elec', [1 : EEG.nbchan], 'threshold', 5, 'norm', 'on', 'measure', 'kurt');

        % Interpolate channels
        EEG = pop_interp(EEG, EEG.chanlocs_original, 'spherical');

        % Reref common average
        EEG = pop_reref(EEG, []);

        % Determine rank of data
        dataRank = sum(eig(cov(double(EEG.data'))) > 1e-6);

        % Epoch data
        EEG = pop_epoch(EEG, {'stimulus'}, [-0.3, 1], 'newname', [id '_epoched'], 'epochinfo', 'yes');
        EEG = pop_rmbase(EEG, [-200, 0]);

        % Reject epochs
        [EEG, rejected_epochs] = pop_autorej(EEG, 'nogui', 'on');
        EEG.trialinfo(rejected_epochs, :) = [];


        % Runica & ICLabel
        EEG = pop_runica(EEG, 'extended', 1, 'interrupt', 'on', 'PCA', dataRank);
        EEG = iclabel(EEG);

        % Find nobrainer
        EEG.nobrainer = find(EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.3 | EEG.etc.ic_classification.ICLabel.classifications(:, 3) > 0.3);

        % Save IC set
        pop_saveset(EEG, 'filename', [id, '_icset.set'], 'filepath', PATH_ICSET, 'check', 'on');

        % Remove components
        EEG = pop_subcomp(EEG, EEG.nobrainer, 0);

        % Save clean data
        pop_saveset(EEG, 'filename', [id, '_cleaned.set'], 'filepath', PATH_AUTOCLEANED, 'check', 'on');


end