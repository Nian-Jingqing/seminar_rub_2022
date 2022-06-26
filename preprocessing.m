

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

        % Behavior matrix
        behavior = [];

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
            behavior(counter, :) = [cnd, cresp, resp, acc, rt];

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

        % Behavior columns: 1=condition | 2=correct response | 3=response | 4=accuracy | 5=rt

        % Check number of trials and add behavioral data to trialinfo
        if size(EEG.trialinfo, 1) == 648 & size(behavior, 1) == 648
            EEG.trialinfo = horzcat(EEG.trialinfo, behavior);
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