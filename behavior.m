clear all;

% Paths
PATH_BEHAVIOR = '/mnt/data_fast/rub_seminar_2022/eprime_export/';


% Get file list
fl = dir([PATH_BEHAVIOR, '*.xlsx']);

counter = 0;
behav_res = [];

% Iterate datasets
for s = 1 : numel(fl)
    
    % Get id
    tmp = regexp(fl(s).name,'\d*','Match');
    id = tmp{1};

    % Open trial table
    trial_table = readtable([PATH_BEHAVIOR, 'S', id, '_eprime.xlsx'], 'ReadVariableNames', true);

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
        behav_res(counter, :) = [str2num(id), t - 24, cnd, cresp, resp, acc, rt];

    end

end