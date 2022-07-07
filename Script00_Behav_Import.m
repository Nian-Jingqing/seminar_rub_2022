%% EEG-S-2022, Script00_BehavioralData Import

PATHIN = '/Users/klatt/Documents/Lehre/EEG_S_2022/data/EPrime';
PATHOUT = '/Users/klatt/Documents/Lehre/EEG_S_2022/data/EPrime/aggregated_data';

SUBJECT = {'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S10' 'S11' 'S12' 'S13' 'S15' 'S34' 'S35' 'S36' 'S37' 'S38'}; %subject ID


%keep columns depending on condition
%Fixation_CRESP_SubTrial_ and TrialDisplay_CRESP and KeytoPress_SubTrial all indicate correct response
keepcolumn_Aonly = {'Condition', 'trialnum', 'TargetPos_SubTrial_', 'KeytoPress_SubTrial_', 'Fixation_RESP_SubTrial_', 'TrialDisplayAonly_RESP', ...
    'Fixation_ACC_SubTrial_','TrialDisplayAonly_ACC', 'Fixation_RT_SubTrial_', 'TrialDisplayAonly_RT'};
keepcolumn_Vonly = {'Condition', 'trialnum' 'TargetPos_SubTrial_', 'KeytoPress_SubTrial_', 'Fixation_RESP_SubTrial_', 'TrialDisplay_RESP', ...
      'Fixation_ACC_SubTrial_', 'TrialDisplay_ACC', 'Fixation_RT_SubTrial_', 'TrialDisplay_RT'};
keepcolumn_Auvi = {'Condition', 'trialnum', 'TargetPos_SubTrial_', 'KeytoPress_SubTrial_', 'Fixation_RESP_SubTrial_', 'TrialDisplayAuvi_RESP', ...
     'Fixation_ACC_SubTrial_', 'TrialDisplayAuvi_ACC'  'Fixation_RT_SubTrial_', 'TrialDisplayAuvi_RT'};

 
 for sub = 1:length(SUBJECT)
     
     %load table
     baseFileName = [SUBJECT{sub} '_eprime.xlsx'];
     folder = PATHIN;
     fullFileName = fullfile(folder, baseFileName);
     
     if exist(fullFileName, 'file')
         opts = detectImportOptions(fullFileName);
         opts.VariableNamesRange = "A3:KW3"; 
         %opts.VariableNames = cell(1,309);
         opts.DataRange = 'A4:KW675';
         columns = [101 203 215 227 91 193 205 217];
         for c = 1:length(columns)
             opts.VariableTypes{columns(c)} = 'double'; 
         end
         % File exists.  Read it into a table.
         %fulltable = readtable(fullFileName,'Range', 'A3:KW675');
         fulltable = readtable(fullFileName,opts, 'ReadVariableNames', true);
     else
         % File does not exist.  Warn the user.
         errorMessage = sprintf('Error: file not found:\n\n%s', fullFileName);
         uiwait(errordlg(errorMessage));
         return;
     end
     
     %delete practice trials
     fulltable(strcmp(fulltable.Procedure_Block_, 'PracticeProc'),:) = [];
     
     %add column with continuous trial counter to fulltable
     trialnum = ones(size(fulltable,1),1);
     trialnum(1:end) = 1:size(fulltable,1);
     fulltable = addvars(fulltable,trialnum,'Before','Condition');
     
     %%%%%%%%%%%%%%%%%%%%%
     %%% AONLY
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     %Aonly-Table
     Aonlytable = fulltable; %create copy
     Aonlytable = Aonlytable(strcmp(Aonlytable.Condition, 'Aonly'),:); %keep only trials from Aonly condition
     
     %exclude columns we don't need
     match = 0;
     clear resp_idx
     for i = 1:length(Aonlytable.Properties.VariableNames)
         if any(strcmp(Aonlytable.Properties.VariableNames(i), keepcolumn_Aonly))
             match = match +1;
             resp_idx(match) = i;
         end
     end
     Aonlytable = Aonlytable(:,resp_idx);
     
     order = keepcolumn_Aonly;
     [~, ind] = ismember(order, Aonlytable.Properties.VariableNames);
     Aonlytable = Aonlytable(:,ind);
     
     %%%%%%%%%%%%%%%%%%%%%
     %%% VONLY
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     %Vonly-Table
     Vonlytable = fulltable; %create copy
     Vonlytable = Vonlytable(strcmp(Vonlytable.Condition, 'Vonly'),:); %keep only trials from Vonly condition
     
     %exclude columns we don't need
     match = 0;
     clear resp_idx
     for i = 1:length(Vonlytable.Properties.VariableNames)
         if any(strcmp(Vonlytable.Properties.VariableNames(i), keepcolumn_Vonly))
             match = match +1;
             resp_idx(match) = i;
         end
     end
     Vonlytable = Vonlytable(:,resp_idx);
     
     order = keepcolumn_Vonly;
     [~, ind] = ismember(order, Vonlytable.Properties.VariableNames);
     Vonlytable = Vonlytable(:,ind);
     
     %%%%%%%%%%%%%%%%%%%%%
     %%% AUVI
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     %Auvi-Table
     Auvitable = fulltable; %create copy
     Auvitable = Auvitable(strcmp(Auvitable.Condition, 'Auvi'),:); %keep only trials from Auvi condition
     
     %exclude columns we don't need
     match = 0;
     clear resp_idx
     for i = 1:length(Auvitable.Properties.VariableNames)
         if any(strcmp(Auvitable.Properties.VariableNames(i), keepcolumn_Auvi))
             match = match +1;
             resp_idx(match) = i;
         end
     end
     Auvitable = Auvitable(:,resp_idx);
     
     order = keepcolumn_Auvi;
     [~, ind] = ismember(order, Auvitable.Properties.VariableNames);
     Auvitable = Auvitable(:,ind);
     
     
     %%%%%%%%%%%%%
     %%%% AGGREGATE
     %one matrix with trial-counter / condition / accuracy / RT
     
     %rename vars so that we can concacenate Aonly, Vonly, and Auvi-Table
     Aonlytable.Properties.VariableNames([6 8 10]) = {'TrialDisplay_RESP' 'TrialDisplay_ACC' 'TrialDisplay_RT'};
     Auvitable.Properties.VariableNames([6 8 10]) =  {'TrialDisplay_RESP' 'TrialDisplay_ACC' 'TrialDisplay_RT'};
     
     %concacenate
     AllTrialsTable = [Aonlytable; Auvitable; Vonlytable];
     AllTrialsTable = sortrows(AllTrialsTable,2); %sort according to trial number (order of conditions is not the same in all subjects)
     
     %rename some additional columns for better data handling
     AllTrialsTable.Properties.VariableNames([3 4 5 7 9]) =  {'TargetPos' 'KeytoPress' 'Fixation_RESP' 'Fixation_ACC' 'Fixation_RT'};
      
     
     %%%%%%%%%%%
     
     %flag trials for which two responses were recorded
     for trial = 1:size(AllTrialsTable,1)
         
         if AllTrialsTable.Fixation_RT(trial) ~= 0 & AllTrialsTable.TrialDisplay_RT(trial) ~= 0
             
             AllTrialsTable.flagdoubleRESP(trial) = 1;
             
         else
             
             AllTrialsTable.flagdoubleRESP(trial) = 0;
             
         end
     end
     
     %%% Check accuracy
     
     %add accuracy column to table
     ACC = nan(size(AllTrialsTable,1),1);
     AllTrialsTable = addvars(AllTrialsTable,ACC,'Before','Fixation_ACC', 'NewVariableNames', 'ACC');
     
     %fill ACC variable
     for trial = 1:size(AllTrialsTable,1)
         
         if AllTrialsTable.KeytoPress(trial) == "" & AllTrialsTable.Fixation_RESP(trial) == "" & AllTrialsTable.TrialDisplay_RESP(trial) == ""
             
             AllTrialsTable.ACC(trial) = 1; %correct no-go
             
         elseif AllTrialsTable.KeytoPress(trial) == "" & (AllTrialsTable.Fixation_RESP(trial) ~= "" | AllTrialsTable.TrialDisplay_RESP(trial) ~= "")
             
             AllTrialsTable.ACC(trial) = 0; %incorrect go
             
         elseif strcmp(AllTrialsTable.KeytoPress(trial), 'm') & (strcmp(AllTrialsTable.Fixation_RESP(trial), 'm') | strcmp(AllTrialsTable.TrialDisplay_RESP(trial), 'm'))
             
             AllTrialsTable.ACC(trial) = 1; %correct target-right identification
             
         elseif strcmp(AllTrialsTable.KeytoPress(trial), 'n') & (strcmp(AllTrialsTable.Fixation_RESP(trial), 'n') | strcmp(AllTrialsTable.TrialDisplay_RESP(trial), 'n'))
             
             AllTrialsTable.ACC(trial) = 1; %correct target-left identification
             
         elseif strcmp(AllTrialsTable.KeytoPress(trial), 'm') & (strcmp(AllTrialsTable.Fixation_RESP(trial), 'n') | strcmp(AllTrialsTable.TrialDisplay_RESP(trial), 'n'))
             
             AllTrialsTable.ACC(trial) = 0; %incorrect go (target-right classified as target-left)
             
         elseif strcmp(AllTrialsTable.KeytoPress(trial), 'n') & (strcmp(AllTrialsTable.Fixation_RESP(trial), 'm') | strcmp(AllTrialsTable.TrialDisplay_RESP(trial), 'm'))
             
             AllTrialsTable.ACC(trial) = 0; %incorrect go (target-left classified as target-right)
             
         elseif any(strcmp(AllTrialsTable.KeytoPress(trial), {'n', 'm'})) & (strcmp(AllTrialsTable.Fixation_RESP(trial), "") | strcmp(AllTrialsTable.TrialDisplay_RESP(trial), ""))
             
             AllTrialsTable.ACC(trial) = 0; %incorrect no-go (target-present classified as target-absent)
             
         end
     end
     
     %if trial was flagged to have two responses, overwrite ACC to 0
     for trial = 1:size(AllTrialsTable,1)
         if AllTrialsTable.flagdoubleRESP(trial) == 1
             AllTrialsTable.ACC(trial) = 0;
         end
     end
     
     %%% merge RT columns
     %RT occuring during fixation (in ISI) is returned relative to fixation
     %cross onset, not relative to stimulus onset; stimulus duration was 500
     %ms; thus we add 500 ms to all fixation-recorded RTs
     
     %add RT column
     RT_column = nan(size(AllTrialsTable,1),1);
     AllTrialsTable = addvars(AllTrialsTable,RT_column,'Before','Fixation_RT', 'NewVariableNames', 'RT');
     
     for trial = 1:size(AllTrialsTable,1)
         
         if AllTrialsTable.Fixation_RT(trial) == 0 & AllTrialsTable.TrialDisplay_RT(trial) == 0 %no-response
             
             AllTrialsTable.RT(trial) = 0;
             
         elseif AllTrialsTable.Fixation_RT(trial) > 0 & AllTrialsTable.TrialDisplay_RT(trial) > 0 %if two responses were recorded
             
             AllTrialsTable.RT(trial) = NaN;
             
         elseif AllTrialsTable.Fixation_RT(trial) > 0 & AllTrialsTable.TrialDisplay_RT(trial) == 0 %response during fixation period
             
             AllTrialsTable.RT(trial) = AllTrialsTable.Fixation_RT(trial) + 500; %add 500 ms to RTs recorded during fixation period
             
         elseif AllTrialsTable.Fixation_RT(trial) == 0 & AllTrialsTable.TrialDisplay_RT(trial) > 0 %response while trial display was still on screen
             
             AllTrialsTable.RT(trial) = AllTrialsTable.TrialDisplay_RT(trial); %copy value to RT column
             
         end
         
     end %end for-loop
     
     %save data table for subject sub
     writetable(AllTrialsTable,[PATHIN filesep SUBJECT{sub} '_Auvi_Behav_allTrials.csv']) 
      
     %% compute mean RT per condition (for correct, target-present trials)
     
     CorrectTargets_Aonly = AllTrialsTable(AllTrialsTable.ACC == 1 & strcmp(AllTrialsTable.TargetPos, 'none') == 0 ...
         & AllTrialsTable.flagdoubleRESP == 0 & strcmp(AllTrialsTable.Condition, 'Aonly') == 1, :);
     RT(sub,1) = mean(CorrectTargets_Aonly.RT);
     
     CorrectTargets_Vonly = AllTrialsTable(AllTrialsTable.ACC == 1 & strcmp(AllTrialsTable.TargetPos, 'none') == 0 ...
         & AllTrialsTable.flagdoubleRESP == 0 & strcmp(AllTrialsTable.Condition, 'Vonly') == 1, :); 
     RT(sub,2) = mean(CorrectTargets_Vonly.RT);
     
     CorrectTargets_Auvi = AllTrialsTable(AllTrialsTable.ACC == 1 & strcmp(AllTrialsTable.TargetPos, 'none') == 0 ...
         & AllTrialsTable.flagdoubleRESP == 0 & strcmp(AllTrialsTable.Condition, 'Auvi') == 1, :);
     RT(sub,3) = mean(CorrectTargets_Auvi.RT);   
     
     %% compute mean ACC per condition 
     Trials_Aonly = AllTrialsTable(strcmp(AllTrialsTable.Condition, 'Aonly') == 1, :);
     Trials_Vonly = AllTrialsTable(strcmp(AllTrialsTable.Condition, 'Vonly') == 1, :);
     Trials_Auvi = AllTrialsTable(strcmp(AllTrialsTable.Condition, 'Auvi') == 1, :);
     
     %compute accuracy for all trial types (target-left, target-right, target-absent)
     ACC_all(sub,1) = (sum(Trials_Aonly.ACC)/size(Trials_Aonly,1)) * 100;
     ACC_all(sub,2) = (sum(Trials_Vonly.ACC)/size(Trials_Vonly,1)) * 100;
     ACC_all(sub,3) = (sum(Trials_Auvi.ACC)/size(Trials_Auvi,1)) * 100;
     
     %compute accuracy separately for target-present trials vs target-absent trials
     
     %%% target-present trials
     Targets_Aonly = AllTrialsTable(strcmp(AllTrialsTable.Condition, 'Aonly') == 1 & strcmp(AllTrialsTable.TargetPos, 'none') == 0, :);
     Targets_Vonly = AllTrialsTable(strcmp(AllTrialsTable.Condition, 'Vonly') == 1 & strcmp(AllTrialsTable.TargetPos, 'none') == 0, :);
     Targets_Auvi = AllTrialsTable(strcmp(AllTrialsTable.Condition, 'Auvi') == 1 & strcmp(AllTrialsTable.TargetPos, 'none') == 0, :);
    
     ACC_targets(sub,1) = ( sum(Targets_Aonly.ACC) / size(Targets_Aonly,1) ) * 100; 
     ACC_targets(sub,2) = ( sum(Targets_Vonly.ACC) / size(Targets_Vonly,1) ) * 100; 
     ACC_targets(sub,3) = ( sum(Targets_Auvi.ACC) / size(Targets_Auvi,1) ) * 100; 
     
     %%% target-absent trials
     NoTargets_Aonly = AllTrialsTable(strcmp(AllTrialsTable.Condition, 'Aonly') == 1 & strcmp(AllTrialsTable.TargetPos, 'none') == 1, :);
     NoTargets_Vonly = AllTrialsTable(strcmp(AllTrialsTable.Condition, 'Vonly') == 1 & strcmp(AllTrialsTable.TargetPos, 'none') == 1, :);
     NoTargets_Auvi = AllTrialsTable(strcmp(AllTrialsTable.Condition, 'Auvi') == 1 & strcmp(AllTrialsTable.TargetPos, 'none') == 1, :);
    
     ACC_notargets(sub,1) = ( sum(NoTargets_Aonly.ACC) / size(NoTargets_Aonly,1) ) * 100; 
     ACC_notargets(sub,2) = ( sum(NoTargets_Vonly.ACC) / size(NoTargets_Vonly,1) ) * 100; 
     ACC_notargets(sub,3) = ( sum(NoTargets_Auvi.ACC) / size(NoTargets_Auvi,1) ) * 100; 
     
   
     
 end
     
%%%%% SAVE
save([PATHOUT filesep 'Behav_all_subj.mat'], 'RT', 'ACC_all', 'ACC_targets', 'ACC_notargets')
csvwrite([PATHOUT filesep 'RT_all_subj.csv'], RT)
csvwrite([PATHOUT filesep 'ACC_alltrials_subj.csv'], ACC_all)
csvwrite([PATHOUT filesep 'ACC_targets_subj.csv'], ACC_targets)     
csvwrite([PATHOUT filesep 'ACC_notargets_subj.csv'], ACC_notargets)  

% save as table (for JASP)
RTtable = array2table(RT,'VariableNames', {'Aonly', 'Vonly', 'Auvi'});
ACCtable = array2table(ACC_all,'VariableNames', {'Aonly', 'Vonly', 'Auvi'});
ACC_targets_table = array2table(ACC_targets,'VariableNames', {'Aonly', 'Vonly', 'Auvi'});
ACC_notargets_table = array2table(ACC_notargets,'VariableNames', {'Aonly', 'Vonly', 'Auvi'});

writetable(RTtable, [PATHOUT filesep 'RT_all_table_subj.csv'])
writetable(ACCtable, [PATHOUT filesep 'ACC_alltrials_table_subj.csv'])
writetable(ACC_targets_table, [PATHOUT filesep 'ACC_targets_alltrials_table_subj.csv'])
writetable(ACC_notargets_table, [PATHOUT filesep 'ACC_notargets_alltrials_table_subj.csv'])