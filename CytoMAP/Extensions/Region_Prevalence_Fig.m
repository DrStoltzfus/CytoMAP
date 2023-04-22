function Region_Prevalence_Fig(app)


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
    
    regs = listdlg('PromptString',{'Select Regions to Plot:'}, 'ListString', RegNames);
    
    
%     grpNMS = {'TP2 A'; ...
%               'TP2 B'; ...
%               'TP2 C'; ...
%               'TP2 D'; ...
%               'TP3 A'; ...
%               'TP3 B'; ...
%               'TP3 C'; ...
%               'TP3 D'};
%     grpNMS = {'TP2 A'; ...
%             'TP2 B'; ...
%             'TP2 C'; ...
%             'TP2 D'};
    grpNMS = {'A'; ...
              'B'; ...
              'C'; ...
              'E'};


% % %       
% % %     grpNMS = {'Sample A'; ...
% % %               'Sample B'; ...
% % %               'Sample C'};
    Ngroups = numel(grpNMS);

    smplnmsHDR = HDR(2:end)';
    smplnms = fieldnames(app.data);
    IND = contains(strrep(smplnms, '_', ' '), smplnmsHDR);
    smplnms = smplnms(IND)
    
    smplVOL = zeros(numel(smplnms), 1);
    imgVOL = zeros(numel(smplnms), 1);
    for smpl_i = 1:numel(smplnms)
        % put the volume with the right sample
        IND = find(strcmp(smplnmsHDR, strrep(smplnms{smpl_i}, '_', ' ')));
% % %         smplVOL(IND) = app.data.(smplnms{smpl_i}).MetaData.SampleVolume;

% % %         vol_i = app.data.(smplnms{smpl_i}).MetaData.SampleVolumeHistoric(1);
        vol_i = app.data.(smplnms{smpl_i}).MetaData.SampleVolumeHistoric(5);

        vol_f = app.data.(smplnms{smpl_i}).MetaData.SampleVolumeHistoric(end);
        if isnan(vol_f)
            vol_f = app.data.(smplnms{smpl_i}).MetaData.SampleVolumeHistoric(end-1);
        end
        
% % %         if vol_i==0
% % %             vol_i = app.data.(smplnms{smpl_i}).MetaData.SampleVolumeHistoric(2);
% % %         end
        
        smplVOL(IND) = (vol_f-vol_i)/vol_i;


        imgVOL(IND) = app.data.(smplnms{smpl_i}).MetaData.ImagedVolume;
    end
    
    % Pull each group
    groupIND = zeros(Ngroups, size(datMAT,1));
    for group_i = 1:numel(grpNMS)
        groupIND(group_i,:) = startsWith(HDR, grpNMS{group_i});
    end
    
% %     groupIND(1,:) = startsWith(HDR, 'TP2 A');
% %     groupIND(2,:) = startsWith(HDR, 'TP2 B');
% %     groupIND(3,:) = startsWith(HDR, 'TP2 C');
% %     groupIND(4,:) = startsWith(HDR, 'TP2 D');
% %     
% %     groupIND(5,:) = startsWith(HDR, 'TP3 A');
% %     groupIND(6,:) = startsWith(HDR, 'TP3 B');
% %     groupIND(7,:) = startsWith(HDR, 'TP3 C');
% %     groupIND(8,:) = startsWith(HDR, 'TP3 D');
% %     [3, 6, 7, 10, 13, 14] [1,2,4,5,8,9,11,12]
    
    %% Plot one region at a time
    if true
        for reg_i = regs
            %% Plot regiond prevalence by sample
            fig = figure;
            clf
            ax = gca;
            fig.Color = app.GUIOPTS.bgclr;
            fig.InvertHardcopy = 'off';
            ax.Title.String = ['Region: ' RegNames{reg_i}];

            hold on
            clrs = {'r', 'b', 'g', 'm','r', 'b', 'g', 'm'};
            datTMP = cell(Ngroups,2);
            for rows_i = 1:Ngroups
                datTMP{rows_i,1} = datMAT(groupIND(rows_i, :)==1, reg_i);
                datTMP{rows_i,2} = ones(numel(datTMP{rows_i,1}),1).*rows_i;
                for smpl_i = 1:sum(groupIND(rows_i, :))
                    group = HDR(groupIND(rows_i, :)==1);
                    plot(rows_i, datTMP{rows_i,1}(smpl_i), ...
                        '.', 'DisplayName', group{smpl_i}, 'Color', clrs{rows_i}, 'MarkerSize', 20)
                end
            end        
            boxplot(cell2mat(datTMP(:,1)), cell2mat(datTMP(:,2)) , 'Color', 'k')

            xlim([0 Ngroups+1])
            ax.XTick = 1:numel(grpNMS);
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
                    plot(smplVOL(smpl_j), datTMP{rows_i,1}(smpl_i), ...
                        '.', 'DisplayName', group{smpl_i}, 'Color', clrs{rows_i}, 'MarkerSize', 20)
                end
            end

            plot(smplVOL(IND_Order), cell2mat(datTMP(:,1)), 'ow')
            plot([0, 0], [-200, 200], '--k')
            
            % Linear Fit
            [R, P] = corr(smplVOL(IND_Order), cell2mat(datTMP(:,1)));
            Linfit = polyfit(smplVOL(IND_Order),cell2mat(datTMP(:,1)),1);
            X = min(smplVOL):((max(smplVOL)-min(smplVOL))/10):max(smplVOL);
            Y = polyval(Linfit,X);
            plot(X, Y, ':k', 'DisplayName', 'Linear Fit')
            
%%%%%
% Swap X and Y
% only use time point 2


%%%%%
% %             % Exponential Fit 
% %             X = smplVOL(IND_Order);
% %             Y = cell2mat(datTMP(:,1));
% % % %             mdl = fittype('a+b*exp(-c*(x-xo))', 'independent', 'x');
% %             mdl = @(a, b, c, d, x) a+b*exp(-c*(x-d));
% %             
% %             startpts = [0; 9; 1; -0.5];
% %             [f0, gof, output] = fit(X,Y,mdl,'StartPoint',startpts);
% % f0;
% % gof;
% % output;
% %             R = sqrt(gof.rsquare);
% %             P = gof.rmse;
% %             xx = linspace(min(X)+0.1*min(X),max(X)+0.1*max(X),50);
% %             plot(xx,f0(xx),':r', 'DisplayName', 'Exponential Fit');
            % Nonlinear regression
            % Find some sort of p value for the exponential fit
%             mdl = @(p, x) p(1) + p(2)*exp(-p(3)*(x-p(4)));
%             opts = statset('nlinfit');
% opts;
%             opts.MaxFunEvals = 600;
%             opts.MaxIter = 400;
%             opts.TolFun = 1e-6;
%             opts.TolX = 1e-6;
%             opts.Robust = 'off';
%             [ahat,r,J,cov,mse] = nlinfit(X,Y,mdl,startpts, opts);
% ahat;
% r;
% J;
% cov;
% mse;
%             ci = nlparci(ahat,r,'Jacobian',J);
% mdl(ahat, xx);
%             plot(xx,mdl(ahat, xx),':b', 'DisplayName', 'Nonlinear regression');

            
            ax.Title.String = ['Region: ' RegNames{reg_i} ...
                ' R = ' num2str(R) ...
                ' R^2 = ' num2str(R^2) ...
                ' p value = ' num2str(P)];
            
% %             ax.Title.String = ['Region: ' num2str(reg_i) ...
% %                 ' R = ' num2str(R) ...
% %                 ' R^2 = ' num2str(R^2) ...
% %                 ' p value = ' num2str(P)];

            ax.YLim = [-2, 100];
            ax.YLabel.String = 'Region Prevalence';
            ax.XLabel.String = 'Fold Change in Tumor Volume (vol_f-vol_i)/vol_i';
            box off

        end
    end
    
    %% Combine regions
% %     regs = [3,7,6];
    
        %% Plot regiond prevalence by sample
        fig = figure;
        clf
        ax = gca;
        fig.Color = app.GUIOPTS.bgclr;
        fig.InvertHardcopy = 'off';
        ax.Title.String = ['Region: ' num2str(regs)];

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
        ax.XTick = 1:numel(grpNMS);
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
                plot(smplVOL(smpl_j), datTMP{rows_i,1}(smpl_i), ...
                    '.', 'DisplayName', group{smpl_i}, 'Color', clrs{rows_i}, 'MarkerSize', 20)
            end
        end

        plot(smplVOL(IND_Order), cell2mat(datTMP(:,1)), 'ow')
        
        [R, P] = corr(smplVOL(IND_Order), cell2mat(datTMP(:,1)));
        Linfit = polyfit(smplVOL(IND_Order),cell2mat(datTMP(:,1)),1);
        X = min(smplVOL):((max(smplVOL)-min(smplVOL))/10):max(smplVOL);
        Y = polyval(Linfit,X);
        plot(X, Y, ':k', 'DisplayName', 'Linear Fit')
               
% %         ax.Title.String = ['Region: ' RegNames ...
% %             ' R = ' num2str(R) ...
% %             ' R^2 = ' num2str(R^2) ...
% %             ' p value = ' num2str(P)];
        
        ax.Title.String = ['Region: ' num2str(regs) ...
            ' R = ' num2str(R) ...
            ' R^2 = ' num2str(R^2) ...
            ' p value = ' num2str(P)];

        ax.YLabel.String = 'Region Prevalence';
        ax.XLabel.String = 'Fold Change in Tumor Volume (vol_f-vol_i)/vol_i';
        box off



end