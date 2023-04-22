function Region_Prevalence_Ratio_Fig(app)


    %% Apply the following figure settings to the currently selected figure
    ax = gca;
    if ~iscell(ax.Title.String)
        ax.Title.String = {ax.Title.String};
    end
    PltDat = get(gca,'Children');
    INDRMV = [];
    for i=1:numel(PltDat)
        if strcmp(PltDat(i).Tag, 'boxplot')
            if isempty(INDRMV)
                INDRMV = i;
            else
                INDRMV = [INDRMV ,i];
            end
%                         PltDat(i) =[];
        end
    end
    PltDat(INDRMV) =[];
    inclX = 0;
    for i=1:numel(PltDat)
        X = PltDat(i).XData;
        Y = PltDat(i).YData;
        YLBL = PltDat(i).DisplayName;
        YLBL = strrep(YLBL, ',', '');
        XLBL = ax.XLabel.String;
        XLBL = strrep(XLBL, ',', '');
        if i==1
            HDR = {XLBL, YLBL};
            datMAT = [X; Y];
        else             
            % deal with non-unique X-Y pairs
            if size(datMAT, 2) < numel(Y)
                % pad datMAT with NaNs
                datMAT = [datMAT, NaN(size(datMAT, 1), numel(Y) - size(datMAT, 2))];
                datMAT = [datMAT; X; Y];
                HDR = [HDR, {XLBL}, {YLBL}];
                inclX = 1;
            elseif size(datMAT, 2) > numel(Y)
                % pad Y with NaNs
                Y = [Y, NaN(1, size(datMAT, 2) - numel(Y))];
                X = [X, NaN(1, size(datMAT, 2) - numel(X))];
                datMAT = [datMAT; X; Y];
                HDR = [HDR, {XLBL}, {YLBL}];
                inclX = 1;
            elseif size(datMAT, 2) == numel(Y) && inclX == 0
                datMAT = [datMAT; Y];
                HDR = [HDR, {YLBL}];
            elseif  size(datMAT, 2) == numel(Y) && inclX == 1
                datMAT = [datMAT; X; Y];
                HDR = [HDR, {XLBL}, {YLBL}];
            end
        end
    end
    
    RegNames = ax.XTickLabel;
    
    regsA = listdlg('PromptString',{'Select A Regions for A/B:'}, 'ListString', RegNames);
    
    regsB = listdlg('PromptString',{'Select B Regions for A/B:'}, 'ListString', RegNames);
    
% % %     regs = listdlg('PromptString',{'Select B Regions for A/B:'}, 'ListString', RegNames);

    
% %     grpNMS = {'TP2 A'; ...
% %               'TP2 B'; ...
% %               'TP2 C'; ...
% %               'TP2 D'; ...
% %               'TP3 A'; ...
% %               'TP3 B'; ...
% %               'TP3 C'; ...
% %               'TP3 D'};
% % %           
% % %     grpNMS = {'TP2 A'; ...
% % %       'TP2 B'; ...
% % %       'TP2 C'; ...
% % %       'TP2 D'};
          
    grpNMS = {'A'; ...
              'B'; ...
              'C'; ...
              'E'};


    Ngroups = numel(grpNMS);
    
% % %       
% % %     grpNMS = {'Sample A'; ...
% % %               'Sample B'; ...
% % %               'Sample C'};

    smplnmsHDR = HDR(2:end)';
    smplnms = fieldnames(app.data);
    IND = contains(strrep(smplnms, '_', ' '), smplnmsHDR);
    smplnms = smplnms(IND)
%     smplnms = fieldnames(app.data);
    
    
    smplVOL = zeros(numel(smplnms), 1);
    imgVOL = zeros(numel(smplnms), 1);
    for smp_i = 1:numel(smplnms)
        % put the volume with the right sample
        IND = find(strcmp(smplnmsHDR, strrep(smplnms{smp_i}, '_', ' ')));
% % %         smplVOL(IND) = app.data.(smplnms{smpl_i}).MetaData.SampleVolume;

        vol_i = app.data.(smplnms{smp_i}).MetaData.SampleVolumeHistoric(5);
% %         vol_i = app.data.(smplnms{smp_i}).MetaData.SampleVolumeHistoric(1);
% %         if vol_i==0
% %             vol_i = app.data.(smplnms{smp_i}).MetaData.SampleVolumeHistoric(2);
% %         end
        vol_f = app.data.(smplnms{smp_i}).MetaData.SampleVolumeHistoric(end);
        if isnan(vol_f)
            vol_f = app.data.(smplnms{smp_i}).MetaData.SampleVolumeHistoric(end-1);
        end
        smplVOL(IND) = (vol_f-vol_i)/vol_i;


        imgVOL(IND) = app.data.(smplnms{smp_i}).MetaData.ImagedVolume;
    end
    
    smplVOL
    imgVOL
    
    % Pull each group
    groupIND = zeros(Ngroups, size(datMAT,1));
    for group_i = 1:Ngroups
        groupIND(group_i,:) = startsWith(HDR, grpNMS{group_i});
    end
    
    
    %% Combine regions
% %     regs = [3,7,6];
    func_reg_plt(regsA, smplVOL, groupIND, datMAT, grpNMS, app, HDR, Ngroups)
    func_reg_plt(regsB, smplVOL, groupIND, datMAT, grpNMS, app, HDR, Ngroups)
    func_reg_Ratio(regsA, regsB, smplVOL, groupIND, datMAT, grpNMS, Ngroups)


    function func_reg_plt(regs, smplVOL, groupIND, datMAT, grpNMS, app, HDR, Ngroups)
        %% Plot regiond prevalence by sample
        fig = figure;
        clf
        ax = gca;
        fig.Color = app.GUIOPTS.bgclr;
        fig.InvertHardcopy = 'off';
        ax.Title.String = ['Region: ' num2str(regs)];
        
        % add an export data option
        % Create Menu bar
        ExportMenu = uimenu(fig);
        ExportMenu.Text = 'Export';
        % Create an export data option
        ExportPltDat = uimenu(ExportMenu);
        ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
        ExportPltDat.Text = 'Export Plot Data to .csv';

        hold on
        clrs = {'r', 'b', 'g', 'm','r', 'b', 'g', 'm'};
        datTMP = cell(Ngroups,2);
        for rows_i = 1:Ngroups
% % %             datTMP{rows_i,1} = datMAT(groupIND(rows_i, :)==1, reg_i)
            datTMP{rows_i,1} = sum(datMAT(groupIND(rows_i, :)==1, regs),2);
            
            datTMP{rows_i,2} = ones(numel(datTMP{rows_i,1}),1).*rows_i;
            for smpl_i = 1:sum(groupIND(rows_i, :))
                group = HDR(groupIND(rows_i, :)==1);
                plot(rows_i, datTMP{rows_i,1}(smpl_i), ...
                    '.', 'DisplayName', group{smpl_i}, 'Color', clrs{rows_i}, 'MarkerSize', 20)
            end
        end        
        
        boxplot(cell2mat(datTMP(:,1)), cell2mat(datTMP(:,2)) , 'Color', 'k')
        
        xlim([0 Ngroups+1])
        ax.XTick = 1:Ngroups;
        ax.XTickLabel = grpNMS;
        ax.YLabel.String = 'Region Prevalence';
        box off
        ylim([0 100])
        
        %% Now plot region prevalence versus sample volume
        fig = figure;
        clf
        ax = gca;
        fig.Color = app.GUIOPTS.bgclr;
        fig.InvertHardcopy = 'off';
        
        % add an export data option
        % Create Menu bar
        ExportMenu = uimenu(fig);
        ExportMenu.Text = 'Export';
        % Create an export data option
        ExportPltDat = uimenu(ExportMenu);
        ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
        ExportPltDat.Text = 'Export Plot Data to .csv';
        
        hold on
        
        IND_Order = zeros(numel(smplVOL),1);
        n = 0;
        for rows_i = 1:Ngroups
            for smpl_i = 1:sum(groupIND(rows_i, :))
                n = n+1;
                smpl_j = find(groupIND(rows_i, :))-1;
                smpl_j = smpl_j(smpl_i);
                IND_Order(n) = smpl_j;
                
                group = HDR(groupIND(rows_i, :)==1);
                plot(datTMP{rows_i,1}(smpl_i), smplVOL(smpl_j), ...
                    '.', 'DisplayName', group{smpl_i}, 'Color', clrs{rows_i}, 'MarkerSize', 20)
            end
        end
        x_dat = cell2mat(datTMP(:,1));
        y_dat = smplVOL(IND_Order);
        plot(x_dat, y_dat, 'ow')
        
        [R, P] = corr(x_dat, y_dat);
        Linfit = polyfit(x_dat, y_dat,1);
        x_fit = min(x_dat):((max(x_dat)-min(x_dat))/10):max(x_dat);
        y_fit = polyval(Linfit, x_fit);
        plot(x_fit, y_fit, ':k', 'DisplayName', 'Linear Fit')
               
        ax.Title.String = ['Region: ' num2str(regs) ...
            ' R = ' num2str(R) ...
            ' R^2 = ' num2str(R^2) ...
            ' p value = ' num2str(P)];

        ax.XLabel.String = 'Region Prevalence';
        ax.YLabel.String = 'Fold Change in Tumor Volume (vol_f-vol_i)/vol_i';
        box off

    end

    function func_reg_Ratio(regsA, regsB, smplVOL, groupIND, datMAT, grpNMS, Ngroups)
        %% Plot regiond prevalence by sample
        fig = figure;
        clf
        ax = gca;
        fig.Color = app.GUIOPTS.bgclr;
        fig.InvertHardcopy = 'off';
        ax.Title.String = ['Region: ' num2str(regsA) '\' num2str(regsB)];

        % add an export data option
        % Create Menu bar
        ExportMenu = uimenu(fig);
        ExportMenu.Text = 'Export';
        % Create an export data option
        ExportPltDat = uimenu(ExportMenu);
        ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
        ExportPltDat.Text = 'Export Plot Data to .csv';
        
        hold on
        clrs = {'r', 'b', 'g', 'm','r', 'b', 'g', 'm'};
        datTMP_A = cell(Ngroups,2);
        datTMP_B = cell(Ngroups,2);
        for rows_i = 1:Ngroups
                        
            datTMP_A{rows_i,1} = sum(datMAT(groupIND(rows_i, :)==1, regsA),2);
            datTMP_A{rows_i,2} = ones(numel(datTMP_A{rows_i,1}),1).*rows_i;
            
            datTMP_B{rows_i,1} = sum(datMAT(groupIND(rows_i, :)==1, regsB),2);
            datTMP_B{rows_i,2} = ones(numel(datTMP_B{rows_i,1}),1).*rows_i;
            
            for smpl_i = 1:sum(groupIND(rows_i, :))
                group = HDR(groupIND(rows_i, :)==1);
                plot(rows_i, datTMP_A{rows_i,1}(smpl_i)/datTMP_B{rows_i,1}(smpl_i), ...
                    '.', 'DisplayName', group{smpl_i}, 'Color', clrs{rows_i}, 'MarkerSize', 20)
            end
        end        
        
        boxplot(cell2mat(datTMP_A(:,1))./cell2mat(datTMP_B(:,1)), cell2mat(datTMP_A(:,2)) , 'Color', 'k')
        
% % %         xlim([0 Ngroups+1])
        ax.XTick = 1:Ngroups;
        ax.XTickLabel = grpNMS;
        ax.YLabel.String = 'Ratio a.u.';
        ax.YScale = 'log';

        box off
% % %         ylim([0 100])
        
        %% Now plot region prevalence versus sample volume
        fig = figure;
        clf
        ax = gca;
        fig.Color = app.GUIOPTS.bgclr;
        fig.InvertHardcopy = 'off';
        
        % add an export data option
        % Create Menu bar
        ExportMenu = uimenu(fig);
        ExportMenu.Text = 'Export';
        % Create an export data option
        ExportPltDat = uimenu(ExportMenu);
        ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
        ExportPltDat.Text = 'Export Plot Data to .csv';
        
        hold on
        
        IND_Order = zeros(numel(smplVOL),1);
        n = 0;
        for rows_i = 1:Ngroups
            for smpl_i = 1:sum(groupIND(rows_i, :))
                n = n+1;
                smpl_j = find(groupIND(rows_i, :))-1;
                smpl_j = smpl_j(smpl_i);
                IND_Order(n) = smpl_j;
                
                group = HDR(groupIND(rows_i, :)==1);
                plot(datTMP_A{rows_i,1}(smpl_i)/datTMP_B{rows_i,1}(smpl_i), smplVOL(smpl_j), ...
                    '.', 'DisplayName', group{smpl_i}, 'Color', clrs{rows_i}, 'MarkerSize', 20)
            end
        end
        
        rat_AB =  cell2mat(datTMP_A(:,1))./cell2mat(datTMP_B(:,1));
        y_dat = smplVOL(IND_Order);
        plot(rat_AB, y_dat, 'ow')
        
        %% Linear Fit
        [R, P] = corr(rat_AB, y_dat);
        Linfit = polyfit(rat_AB,y_dat, 1);
        x_fit = min(rat_AB):((max(rat_AB)-min(rat_AB))/10):max(rat_AB);
        y_fit = polyval(Linfit,x_fit);
        plot(x_fit, y_fit, ':k', 'DisplayName', 'Linear Fit')
        
        % Exponential Fit 
% %         X = smplVOL(IND_Order);
% %         Y = cell2mat(datTMP(:,1));
% %             mdl = fittype('a+b*exp(-c*(x-xo))', 'independent', 'x');
% %         mdl = @(a, b, c, d, x) a+b*exp(-c*(x-d));
% % 
% %         startpts = [0; 9; 1; -0.5];
% %         [f0, gof, output] = fit(rat_AB,y_dat,mdl,'StartPoint',startpts);
% % f0;
% % gof;
% % output;
% % % %             R = sqrt(gof.rsquare);
% % % %             P = gof.rmse;
% %         xx = linspace(min(rat_AB)+0.1*min(rat_AB),max(rat_AB)+0.1*max(rat_AB),50);
% %         plot(xx,f0(xx),':r', 'DisplayName', 'Exponential Fit');

               
        ax.Title.String = ['Region: (' num2str(regsA) ')\(' num2str(regsB) ')' ...
            ' R = ' num2str(R) ...
            ' R^2 = ' num2str(R^2) ...
            ' p value = ' num2str(P)];

        ax.XLabel.String = 'Ratio a.u.';
        ax.YLabel.String = 'Fold Change in Tumor Volume (vol_f-vol_i)/vol_i';
        ax.XScale = 'log';
        box off

    end


end