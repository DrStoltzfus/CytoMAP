function [tdat, Dreduc_opts] = func_phate(Dat_Pre, opts)
addpath([fileparts(mfilename('fullpath')) filesep 'PHATE']);

    if isempty(Dat_Pre)
        tdat = [];

        Dreduc_opts = cell(2,2);
        % Options for U-MAP
        Dreduc_opts(:,1) = {'Function',...
            'Channel Name'};
        % Defaults
        Dreduc_opts(:,2) = {'PHATE_EB',...
            opts};
        dlg_title = 'Options';
        num_lines = 1;
        vAnswer = inputdlg(Dreduc_opts(:,1),dlg_title,num_lines,Dreduc_opts(:,2));
        if isempty(vAnswer)
            return
        end

        Dreduc_opts(:,2) = vAnswer;

% % %         ListOption = {'PHATE_EB','PHATE_DLA_TREE','PAHTE_TREE','PHATE_mESC'};
% % %         [Selection, ~] = listdlg('ListString',ListOption, 'SelectionMode','single');

    else
        Dreduc_opts = opts;

        %%
        switch Dreduc_opts{1,2}
            case 'PHATE_EB' 
                % PHATE_EB
                Dat_Pre= 10*Dat_Pre + 0.02*rand(size(Dat_Pre));
                tdat = phate(Dat_Pre, 't', 20);
                figure;
                scatter(tdat(:,1), tdat(:,2), 3, 'filled');
                set(gca,'xticklabel',[]);
                set(gca,'yticklabel',[]);
                axis tight
                xlabel 'PHATE1'
                ylabel 'PHATE2'
                title 'PHATE'
                drawnow
            case 'PHATE_DLA_TREE'
                % PHATE_DLA_TREE
                Dat_Pre= 10*Dat_Pre + 0.03*rand(size(Dat_Pre));
                tdat = phate(Dat_Pre, 't', 20, 'gamma', 0);
                figure;
                scatter(tdat(:,1), tdat(:,2), 10, 'filled');
                set(gca,'xticklabel',[]);
                set(gca,'yticklabel',[]);
                axis tight
                xlabel 'PHATE1'
                ylabel 'PHATE2'
                drawnow
            case 'PAHTE_TREE'
                % PAHTE_TREE
                Dat_Pre= 10*Dat_Pre + 0.03*rand(size(Dat_Pre));
                tdat = phate(Dat_Pre, 't', 32);
                figure;
                scatter(tdat(:,1), tdat(:,2), 10, 'filled');
            case 'PHATE_mESC'
                % PHATE_mESC
                Dat_Pre= 20*Dat_Pre + 0.04*rand(size(Dat_Pre));
                tdat = phate(Dat_Pre);
                figure;
                scatter(tdat(:,1),tdat(:,2),20,'filled')
                axis tight
                xlabel('PHATE 1')
                ylabel('PHATE 2')
                set(gca,'xticklabel',[])
                set(gca,'yticklabel',[])
                title 'PHATE 2D'
                drawnow
        end
    end
end