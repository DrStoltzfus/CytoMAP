function [tdat, Dreduc_opts] = func_tscan(Dat_Pre, opts)
addpath([fileparts(mfilename('fullpath')) filesep 'TSCAN']);

    if isempty(Dat_Pre)
        tdat = [];


    Dreduc_opts = cell(6,2);
    % Options for Mclust
    Dreduc_opts(:,1) = {'Column',...
        'CovarianceType', ...
        'SharedCovariance',...
        'StartCondition', ...
        'Scale', ...
        'Channel Name'};
    % Defaults
    Dreduc_opts(:,2) = {'2',...
        'diagonal', ...
        'false',...
        '1', ...
        '100', ...
        opts};

    dlg_title = 'Options';
    num_lines = 1;
    vAnswer = inputdlg(Dreduc_opts(:,1),dlg_title,num_lines,Dreduc_opts(:,2));

    if isempty(vAnswer)
        return
    end

    Dreduc_opts(:,2) = vAnswer;
    
    else
        Dreduc_opts = opts;
        scale = str2double(Dreduc_opts{5, 2});
        % Scale the data to be in the 0-100 range
        IND = range(Dat_Pre)<scale;
        Dat_Pre(:, IND) = scale * Dat_Pre(:, IND);
        % Shift the data so that it is positive
        IND = min(Dat_Pre) < scale;
        Dat_Pre(:, IND) = Dat_Pre(:, IND) - min(Dat_Pre(:, IND));

        procdata = preprocess(Dat_Pre);
        lpsmclust = exprmclust(procdata', Dreduc_opts(1:(end-2), :));
        tdat(:,1) = lpsmclust.pcareduceres(:,1);
        tdat(:,2) = lpsmclust.pcareduceres(:,2);
        plotmclust(lpsmclust);
    end
end