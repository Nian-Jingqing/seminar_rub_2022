%% Script 01 - EEG-S-2022
% CODING
% add relevant fields to EEG.event structure (i.e. translate event.type trigger
% information into variables)
% add responses, response times, accuracy (from behavioral data files)

% (c) Laura Klatt & Julian Reiser
% can be shared under CC BY-NC-SA licence (https://creativecommons.org/licenses/by-nc-sa/4.0/)
%%
clear variables; clc;

PATHIN_RAW = '/Users/klatt/Documents/Lehre/EEG_S_2022/data/EEG_RAW/'; %where raw EEG data files are stored
PATHIN_BEHAV = '/Users/klatt/Documents/Lehre/EEG_S_2022/data/EPrime/eprimeExport/'; %where behavioral data files are stored
PATHOUT = '/Users/klatt/Documents/Lehre/EEG_S_2022/data/EEG/'; %where to save coded output files

SUBJECT = {'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S10' 'S11' 'S12' 'S13' 'S15'} %all subjects
SUBJECT = {'S01' 'S02' 'S05' 'S07' 'S08' 'S10' 'S11'}; %subjects with final naming scheme, other still need to be re-named

trial_count = zeros(length(SUBJECT), 3);

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %start eeglab

for sub = 1 :length(SUBJECT)

    %Define subject path based on PATHOUT directory and subject ID of current subject
    Subject_Path = [PATHOUT SUBJECT{sub} filesep];
    
    if ~exist(Subject_Path, "dir")
    
        mkdir(Subject_Path)

    end

    %load raw EEG dataset
    EEG = pop_loadbv(PATHIN_RAW, [SUBJECT{sub} '_AuVi2.vhdr']);

   
    %%% add condition variable
    %AuVi = audio-visual trials
    %Aonly = auditory-only trials
    %Vonly = visual-only trials

    Auvi_trials = 0; Aonly_trials = 0; Vonly_trials = 0; %initialize trial counters   

    for ev = 1:length(EEG.event)
        
        if any(strcmp(EEG.event(ev).type, {'S101', 'S102', 'S103'}))

            EEG.event(ev).condition = 'practice';

        elseif any(strcmp(EEG.event(ev).type(3), {'1' '2' '3'})) && strcmp(EEG.event(ev).type(2), ' ')

            EEG.event(ev).condition = 'AuVi';
            Auvi_trials = Auvi_trials +1;

        elseif any(strcmp(EEG.event(ev).type(3),  {'4' '5' '6'})) && strcmp(EEG.event(ev).type(2), ' ')

            EEG.event(ev).condition = 'Aonly';
            Aonly_trials = Aonly_trials +1;

        elseif any(strcmp(EEG.event(ev).type(3), {'7' '8' '9'})) && strcmp(EEG.event(ev).type(2), ' ')

            EEG.event(ev).condition = 'Vonly';
            Vonly_trials = Vonly_trials +1;
            
        end

    end
    
    %save trial_count for all subjects
    trial_count(sub,1) = Auvi_trials;
    trial_count(sub,2) = Aonly_trials;
    trial_count(sub,3) = Vonly_trials;

    %%% add target-position variable (left vs. right vs. none)
    for ev = 1:length(EEG.event)

        if strcmp(EEG.event(ev).type(2), ' ') %only go through two-digit S-trigger

            if any(strcmp(EEG.event(ev).type(3), {'1', '4', '7'}))

                EEG.event(ev).targetpos = 'left';

            elseif any(strcmp(EEG.event(ev).type(3), {'2', '5', '8'}))

                EEG.event(ev).targetpos = 'right' ;

            elseif any(strcmp(EEG.event(ev).type(3), {'3', '6', '9'}))

                EEG.event(ev).targetpos = 'none';

            end
        end

    end %end for-event-loop

    %%% add target digit (1-9)
    for ev = 1:length(EEG.event)

        if strcmp(EEG.event(ev).type(2), ' ') %only go through two-digit S-trigger

            if strcmp(EEG.event(ev).type(4), '1')

                EEG.event(ev).targetdigit = '1';

            elseif strcmp(EEG.event(ev).type(4), '2')

                EEG.event(ev).targetdigit = '2';

            elseif strcmp(EEG.event(ev).type(4), '3')

                EEG.event(ev).targetdigit = '3';

            elseif strcmp(EEG.event(ev).type(4), '4')

                EEG.event(ev).targetdigit = '4';

            elseif strcmp(EEG.event(ev).type(4), '5')

                EEG.event(ev).targetdigit = '5';

            elseif strcmp(EEG.event(ev).type(4), '6')

                EEG.event(ev).targetdigit = '6';

            elseif strcmp(EEG.event(ev).type(4), '7')

                EEG.event(ev).targetdigit = '7';

            elseif strcmp(EEG.event(ev).type(4), '8')

                EEG.event(ev).targetdigit = '8';

            elseif strcmp(EEG.event(ev).type(4), '9')

                EEG.event(ev).targetdigit = '9';

            end %end if-digit-statement

        end

    end %end for-event loop

    %%%%% load behavioral data file for current subject
    file = [PATHIN_BEHAV SUBJECT{sub} '_Auvi_Behav_allTrials.csv'];
    behav = readtable(file);
    
    %%% find first trial after practice block
    for ev = 1:length(EEG.event)

        if strcmp(EEG.event(ev).condition, 'practice') == 1 && strcmp(EEG.event(ev+1).condition, 'practice') == 0

            start_trigger = ev + 1;

            break

        end

    end

    %%% find index of boundary events
    clear boundary;
    for ev = start_trigger :length(EEG.event)
        if strcmp(EEG.event(ev).type, 'boundary')
            boundary = ev;
        end
    end

    % create event-array excluding boundary events
    clear event_array;
    event_array = start_trigger:length(EEG.event);
    if exist('boundary', 'var')
        event_array(event_array == boundary) = []; %omit boundary events
    end

    %%% add response time
    RTsub = behav.RT; %NaN = trials with two responses, flagged as NaN; otherwise, if no response was given, RT = 0;
    trial = 1; 
    for ev = event_array
        EEG.event(ev).RT = RTsub(trial);
        trial = trial + 1;
    end
    
    %%% add accuracy (correct vs. incorrect)
    ACCsub = behav.ACC;
    trial = 1;
    for ev = event_array
        EEG.event(ev).ACC = ACCsub(trial);
        trial = trial + 1;
    end
    
    %%% add continuous trial counter
    TRIALsub = behav.trialnum;
    trial = 1;
    for ev = event_array
        EEG.event(ev).trialnum = TRIALsub(trial);
        trial = trial + 1;
    end
    %for practice trials set trial counter to 0
    for ev = event_array
        EEG.event(ev).trialnum = 0;
    end

    %%% letterkilla - get rid of non-numeric fields in EEG.event.type
    EEG = letterkilla(EEG);
    EEG = eeg_checkset( EEG );

    % Save as a new dataset
    EEG=pop_saveset(EEG,'filename',[SUBJECT{sub} '_coded_raw.set'],'filepath',Subject_Path);
    
    %get rid of practice trials
    EEG = pop_selectevent( EEG, 'omitcondition',{'practice'},'deleteevents','on');
    
    %%%%% re-code event triggers to include information whether trial is correct
    %%%%% or incorrect
    for ev = 1:length(EEG.event)

        if str2double(EEG.event(ev).type) < 100 && EEG.event(ev).ACC == 1

            EEG.event(ev).type = str2double(EEG.event(ev).type) + 100; %event-fields starting with 1 == correct trials

        elseif str2double(EEG.event(ev).type) < 100 && EEG.event(ev).ACC == 0

            EEG.event(ev).type = str2double(EEG.event(ev).type) + 300; %event-fields starting with 3 == false trials

        end

    end

    EEG = eeg_checkset( EEG );
    
    EEG = pop_saveset(EEG,'filename',[SUBJECT{sub} '_recoded.set'],'filepath',Subject_Path);

end