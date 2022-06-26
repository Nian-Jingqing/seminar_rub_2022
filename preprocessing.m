

% Path variables
PATH_EEGLAB        = '/home/plkn/eeglab2022.0/';
PATH_CODED_DATA    = '/mnt/data_fast/rub_seminar_2022/1_coded_data/';
PATH_ICSET         = '/mnt/data_fast/rub_seminar_2022/2_ic_sets/';
PATH_AUTOCLEANED   = '/mnt/data_fast/rub_seminar_2022/3_autocleaned/';


% Get vhdr file list
fl = dir([PATH_CODED_DATA, '*.set']);

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

        % Load coded data
        EEG = pop_loadset('filename', fl(s).name, 'filepath', PATH_CODED_DATA, 'loadmode', 'all');

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