function [tdat, Dreduc_opts] = func_pca(Dat_Pre, opts)
    if isempty(Dat_Pre)
        tdat = [];
        Dreduc_opts = cell(3,2);
        % Options for U-MAP
        Dreduc_opts(:,1) = {'Name', ...
            'Value', ...
            'Channel Name'};
        % Defaults
        Dreduc_opts(:,2) = {'Rows', ...
            'complete', ...
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

        [coeff,score,~] = pca(Dat_Pre,...
                        Dreduc_opts{1, 2}, ...
                        Dreduc_opts{2, 2});

        tdat = score(:, 1:2);
        figure;
        biplot(coeff(:,1:2),'scores',score(1:300,1:2));
    end
end

