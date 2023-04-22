function [tdat, Dreduc_opts] = func_lle_hlle(Dat_Pre, opts)
addpath([fileparts(mfilename('fullpath')) filesep 'LLE']);

    if isempty(Dat_Pre)
        tdat = [];

        Dreduc_opts = cell(4,2);
        % Options for U-MAP
        Dreduc_opts(:,1) = {'function',...
            'number of neighbors', ...
            'max embedding dimensionality', ...
            'Channel Name'};
        % Defaults
        Dreduc_opts(:,2) = {'LLE',...
            '3', ...
            '2', ...
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

        if strcmp(Dreduc_opts{1, 2}, 'LLE')
            Dat_Pre= 10*Dat_Pre + 0.02*rand(size(Dat_Pre));
            [Y] = lle(Dat_Pre',str2double(Dreduc_opts{2, 2}), str2double(Dreduc_opts{3, 2}));
            tdat = Y';
    % %         size(tdat)
        elseif strcmp(Dreduc_opts{1, 2},'HLLE')
            size(Dat_Pre)
            Dat_Pre= 10*Dat_Pre + 0.05*rand(size(Dat_Pre));
            [Y, ~] = HLLE(Dat_Pre',str2double(Dreduc_opts{2, 2})*130, str2double(Dreduc_opts{3, 2}));
            tdat = Y';
    % %         size(tdat)
        end
    end