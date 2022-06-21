

% Path variables
PATH_EEGLAB        = '/home/plkn/eeglab2022.0/';
PATH_RAW_DATA      = '/mnt/data_fast/rub_seminar_2022/raw_data/';
PATH_ICSET         = '/mnt/data_fast/rub_seminar_2022/ic_sets/';
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
                    stimpos = 3;
                end

                % Get stimulus identity
                stim_id = str2num(EEG.event(e).type(4));

                % Fill trial table
                EEG.trialinfo(counter, :) = [av_cond, stimpos];

                % Rename trigger
                EEG.event(e).type = 'stimulus';

            end
        end

        % Add FCz as empty channel
        EEG.data(end + 1, :) = 0;
        EEG.nbchan = size(EEG.data, 1);
        EEG.chanlocs(end + 1).labels = 'FCz';

        % Add channel locations
        EEG = pop_chanedit(EEG, 'lookup', channel_location_file);

        % Save original channel locations (for later interpolation)
        EEG.chanlocs_original = EEG.chanlocs;

        % Reref to CPz, so that FCz obtains non-interpolated data
        EEG = pop_reref(EEG, 'CPz');

        % Filter
        EEG = pop_basicfilter(EEG, [1 : EEG.nbchan], 'Cutoff', [0.5, 40], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 4, 'RemoveDC', 'on', 'Boundary', 'boundary');
            
        % Bad channel detection
        [EEG, EEG.chans_rejected] = pop_rejchan(EEG, 'elec', [1 : EEG.nbchan], 'threshold', 5, 'norm', 'on', 'measure', 'kurt');

        % Interpolate channels
        EEG = pop_interp(EEG, EEG.chanlocs_original, 'spherical');

        % Reref common average
        EEG = pop_reref(EEG, []);

        % Determine rank of data
        dataRank = sum(eig(cov(double(EEG.data'))) > 1e-6);

        % Epoch data
        EEG = pop_epoch(EEG, {'stimulus'}, [-1, 1.8], 'newname', [num2str(s) '_epoched'], 'epochinfo', 'yes');
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
        pop_saveset(EEG, 'filename', [num2str(s), '_icset.set'], 'filepath', PATH_ICSET, 'check', 'on');

        % Remove components
        EEG = pop_subcomp(EEG, EEG.nobrainer, 0);

        % Save clean data
        pop_saveset(EEG, 'filename', [num2str(s), '_cleaned.set'], 'filepath', PATH_AUTOCLEANED, 'check', 'on');


end