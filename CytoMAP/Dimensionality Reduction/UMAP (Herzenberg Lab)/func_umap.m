function [tdat, Dreduc_opts] = func_umap(Dat_Pre, opts)


    if isempty(Dat_Pre)
        tdat = [];
        Dreduc_opts = cell(5,2);
        % Options for U-MAP
        Dreduc_opts(:,1) = {'n_neighbors', ...
            'min_dist', ...
            'n_components', ...
            'metric', ...
            'Channel Name'};
        % Defaults
        Dreduc_opts(:,2) = {'30', ...
            '0.3', ...
            '2', ...
            'euclidean', ...
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

        addpath([fileparts(mfilename('fullpath')) filesep 'umap']);
        addpath([fileparts(mfilename('fullpath')) filesep 'util']);
        tdat = run_umap(Dat_Pre,...
                        Dreduc_opts{1, 1}, str2double(Dreduc_opts{1, 2}), ...
                        Dreduc_opts{2, 1}, str2double(Dreduc_opts{2, 2}), ...
                        Dreduc_opts{3, 1}, str2double(Dreduc_opts{3, 2}), ...
                        Dreduc_opts{4, 1}, Dreduc_opts{4, 2});
    end
end