function Ripley_Distribution_Analysis(app)
    %% Build a user interface
    % Ripleys_K_Function calculates Ripleys K function
    %
    % Input:
    %   - app - Instance of CytoMAP
    %
    
    if ~Helper.any_sample(app)
        return;
    end

    %  Pull the phenotype names from the data table
    tmp = Helper.populate_table(app, 'smpls', {app.DataN.Value});
    n_pheno = sum(~cellfun('isempty', tmp(:, 3)));
    if n_pheno < size(tmp, 1)
        % More samples than phenotypes
        tbldat = cell(max(size(tmp, 1) + 1, n_pheno+1 ), 5);
        tbldat(1:n_pheno, 1)= {false};
        tbldat(1:n_pheno, 2) = tmp(1:n_pheno, 3);
        tbldat(end, 2) = {'Select All'};
        tbldat(1:numel(app.DataN.Items), 3) = {false};
        tbldat(strcmp(app.DataN.Value, app.DataN.Items), 3) = {true};
        tbldat(end, 3) = {false};
        tbldat(1:numel(app.DataN.Items), 4) = {1};
        tbldat(1:size(tmp, 1), 5) = tmp(:, 5);
        tbldat(end, 5) = {'Select All'};
    else
        % Number of phenotypes equal to size of table
        tbldat = cell(size(tmp, 1), 5);
        tbldat(:, 1)= {false};
        tbldat(1:size(tmp, 1), 2) = tmp(:, 3);
        tbldat(end, 2) = {'Select All'};
        tbldat(1:numel(app.DataN.Items), 3) = {false};
        tbldat(strcmp(app.DataN.Value, app.DataN.Items), 3) = {true};
        tbldat(end, 3) = {false};
        tbldat(1:numel(app.DataN.Items), 4) = {1};
        tbldat(1:size(tmp, 1), 5) = tmp(:, 5);
        tbldat(end, 5) = {'Select All'};
    end

    
    % Find the names of the user defined surfaces
    if max(contains(fieldnames(app.data.(app.DataN.Value).Surfaces), 'UDS'))
        UDSNames = fieldnames(app.data.(app.DataN.Value).Surfaces.UDS);
        UDSNames = [{'Surfaces', 'none'}, UDSNames(:)'];
        if isempty(UDSNames)
            UDSNames = {'Surfaces', 'none'};
        end
    else
        UDSNames = {'Surfaces', 'none'};
    end
    
    UIfig = uifigure('Name', 'Calculate Ripley"s K Function', 'Scrollable', 'on');
    Helper.func_SetCLR(app, UIfig, 'UIfigure');
    UIfig.Position = [10 10 800 800];
    
    t = uitable(UIfig);
    t.Data = tbldat;
    t.Position = [0 50 800 750];
    t.ColumnName = {'Include:', 'Phenotype', 'Include:', 'Group', 'Sample'};
    t.ColumnEditable = [true, false, true, true, false];
    t.ColumnWidth = {75 250 75 75 250};
    
    % Cell Types
    pnmsd = logical([t.Data{:, 1}]);
    pnmsd = t.Data(pnmsd, 2);
    if isempty(pnmsd)
        pnmsd = {'None'};
    else
        pnmsd = [{'None'}, pnmsd'];
    end
    if ~iscell(pnmsd)
        pnmsd = {pnmsd};
    end 
   
    % Select plot difference from
    popdiff = uidropdown(UIfig);
    popdiff.Items = pnmsd;
    popdiff.Value = popdiff.Items{1};
    popdiff.Position = [235+10+100+150 1 100 30];
    Helper.func_SetCLR(app, popdiff, 'UIField')
    % Select Surface label
    lbl = uilabel(UIfig); lbl.Text = sprintf('Plot Difference from:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = [235+10+100+150 32 200 15];
    
    t.CellEditCallback = @(dd, p) switched_smpl(app, dd, p, popdiff);
    
    % Select User defined surfaces drop down
    smplarea = uidropdown(UIfig);
    smplarea.Items = UDSNames(2:end);
    smplarea.Value = UDSNames{2};
    smplarea.Position = [10 1 100 30];
    Helper.func_SetCLR(app, smplarea, 'UIField')
    % Select Surface label
    lbl = uilabel(UIfig); lbl.Text = sprintf('Select Bounding Surface:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = [10 32 200 15];
    
    % Make a radius of neighborhood numeric option
    lbl = uilabel(UIfig); lbl.Text = 'Number of Distance Steps';
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = [235+10+100 32 300 15];
    Iterations = uieditfield(UIfig,'numeric');
    Iterations.Value = 25;
    Iterations.Position = [235+10+100 1 80 30];
    
    % Make a radius of neighborhood numeric option
    lbl = uilabel(UIfig); lbl.Text = 'Maximum Distance';
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = [235+10 32 300 15];
    MaximumDistance = uieditfield(UIfig,'numeric');
    MaximumDistance.Value = 2000;
    MaximumDistance.Position = [235+10 1 80 30];

    % Make a radius of neighborhood numeric option
    lbl = uilabel(UIfig); lbl.Text = 'Beginning Distance';
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = [135+10 32 300 15];
    BeginningDistance = uieditfield(UIfig,'numeric');
    BeginningDistance.Value = 5;
    BeginningDistance.Position = [135+10 1 80 30];
    
    % Create a Start button
    btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) RK_func_backend( ...
        app, t, UIfig, smplarea.Value, Iterations.Value, MaximumDistance.Value, BeginningDistance.Value, popdiff.Value));
    btn.Position = [800-60, 10, 50, 30];
    btn.Text = 'Ok';
    Helper.func_SetCLR(app, btn, 'button')

    function switched_smpl(app, dd, p, popdiff)
        % If the sample was changed
        if p.Indices(2) == 3
            % Make sure that at least one thing is selected
            if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 3)), 3)))
                dd.Data{p.Indices(1), 3} = true;
                return;
            end
            % Process 'Select All' button
            fill_select_all = false;
            if p.Indices(1) == find(strcmp(dd.Data(:, 5), 'Select All'))
                ind = ~cellfun('isempty', dd.Data(:, 3));
                dd.Data(ind, 3) = {logical(p.NewData)};
                % Make sure that at least one thing is selected
                if ~logical(p.NewData)
                    dd.Data{strcmp(dd.Data(:, 5), app.DataN.Value), 3} = true;
                end
                fill_select_all = logical(p.NewData);
            end
            ind = ~cellfun('isempty', dd.Data(:, 3));
            ind(end) = false;
            ind(ind == 1) = logical(cell2mat(dd.Data(ind, 3)));

            % Extract special ones (broadcast logic)
            logic = ~cellfun('isempty', dd.Data(:, 2));
            logic(logic) = ismember(dd.Data(logic, 2), 'Select All');

            keep_special = dd.Data(logic, 1:2);
            NSamp = sum(~cellfun('isempty', dd.Data(:, 3)))-1;
            keep_Group = dd.Data(1:NSamp, 4);
            
            tmp_s = cell(size(dd.Data, 1) - 1, 5);

            tmp_s(:, 1) = {1};
            tmp_s(:, 2) = dd.Data(1:end - 1, 1);
            tmp_s(:, 3) = dd.Data(1:end - 1, 2);

            tmp_s(:, 4) = dd.Data(1:end - 1, 3);
            tmp_s(:, 5) = dd.Data(1:end - 1, 5);
            tmp_s = Helper.populate_table(app, ...
                'smpls', dd.Data(ind, 5), ...
                'prev_table', tmp_s ...
            );

            non_empty_pheno = sum(~cellfun('isempty', tmp_s(:, 2)));
            if  non_empty_pheno < numel(app.DataN.Items)
                % Less phenotypes than samples.
                dd.Data = cell(max(non_empty_pheno + size(keep_special, 1), numel(app.DataN.Items) + 1), 5);

                dd.Data(1:non_empty_pheno, 1) = tmp_s(1:non_empty_pheno, 2);
                dd.Data(1:non_empty_pheno, 2) = tmp_s(1:non_empty_pheno, 3);
                dd.Data(end - size(keep_special, 1) + 1:end, 1:2) = keep_special;

                dd.Data(1:size(tmp_s, 1), 3) = tmp_s(:, 4);
                dd.Data(1:numel(keep_Group), 4) = keep_Group;
                dd.Data(1:size(tmp_s, 1), 5) = tmp_s(:, 5);
                dd.Data(end, [3 5]) = {fill_select_all, 'Select All'};
            else
                % Regular case.
                %Only re-draw the table if the size changes
                %TODo put something like this this in all of the table generation
                %sections to speed up table remaking
                if size(dd.Data, 1) ~= (size(tmp_s, 1) + 1)
                    dd.Data = cell(size(tmp_s, 1) + 1, 5);
                end
                dd.Data(1:end - 1, 1) = tmp_s(:, 2);
                dd.Data(1:end - 1, 2) = tmp_s(:, 3);
                dd.Data(end - size(keep_special, 1) + 1:end, 1:2) = keep_special;

                dd.Data(1:end - 1, 3) = tmp_s(:, 4);
                dd.Data(1:numel(keep_Group), 4) = keep_Group;
                dd.Data(1:end - 1, 5) = tmp_s(:, 5);
                dd.Data(end, [3 5]) = {fill_select_all, 'Select All'};
            end
        elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
            ind = ~cellfun('isempty', dd.Data(:, 2));
            dd.Data(ind, 1) = {logical(p.NewData)};
        end
        % Cell Types
        pnmstmp = logical([dd.Data{:, 1}]);
        pnmstmp = dd.Data(pnmstmp, 2);
        if isempty(pnmstmp)
            pnmstmp = {'None'};
        else
            pnmstmp = [{'None'}, pnmstmp'];
        end
        if ~iscell(pnmstmp)
            pnmstmp = {pnmstmp};
        end
        popdiff.Items = pnmstmp;
    end

    %% Run the Ripley K function
    function RK_func_backend(app, t, UIfig, smplarea, Iterations, MaximumDistance, BeginningDistance, popdiff)
        %% Pull Selected Variables
        % Cell Types
        IND = logical([t.Data{1:(end-1), 1}]);
        pnms = t.Data(1:(end-1), 2);
        pnms = pnms(IND, :);
        
        IND = logical([t.Data{1:(end-1), 3}]);
        Groups = t.Data(1:(end-1), 4); % Groups
        Groups = Groups(IND, :);
        
        IND = logical([t.Data{1:(end-1), 3}]);
        smpls = t.Data(1:(end-1), 5); % Samples
        smpls = smpls(IND, :);
                
        if ~iscell(smpls)
            smpls = {smpls};
        end
        if ~iscell(pnms)
            pnms = {pnms};
        end     
        if ~iscell(Groups)
            Groups = {Groups};
        end
        % Make sure all of the elements of samples are actually sample names
        INDrmv = [];
        for k = 1:numel(smpls)
            if ~ischar(smpls{k})
                INDrmv = [INDrmv, k];
            end
        end
        smpls(INDrmv) = [];
        
        INDrmv = [];
        for k = 1:numel(Groups)
            if isempty(Groups{k})
                INDrmv = [INDrmv, k];
            end
        end
        Groups(INDrmv) = [];
        Groups = cell2mat(Groups);

        %% Calculate Ripleys k-function
        % Details at https://pro.arcgis.com/en/pro-app/tool-reference/spatial-statistics/h-how-multi-distance-spatial-cluster-analysis-ripl.htm
        % define a bounding box around the actual tissue
        
%         smplarea = 'Tissue';
%         Iterations = 25;
%         MaximumDistance = 2000;
%         BeginningDistance = 5;
    
        DistanceIncrement = (MaximumDistance-BeginningDistance)/Iterations;

        % results = cell(numel(smpls)+1, 1)'
        results = {'Sample', 'CellType', 'Distance Threshold', 'sum(kij)', 'L(d)', 'Group'};
        for i = 1:numel(smpls)
            if ~isempty(pnms)
                % Find the names of the to data
                tags = Helper.gate_full2tag(app, pnms, smpls{i});
            end

            % Define the sample area
            if ~strcmp(smplarea, 'none')
                shp = app.data.(smpls{i}).Surfaces.UDS.(smplarea).Surf;
                A = shp.area;
            end

            for j = 1:numel(pnms)
                tag = tags{j};
                tagLogic = app.data.(smpls{i}).AllCells.(tag);

                dat = app.data.(smpls{i}).AllCells(tagLogic==1, {'X', 'Y', 'Z'});
% % %                 dat = dat(1:10, :); % For testing

                if ~isempty(dat)
                    if strcmp(smplarea, 'none')
                        if max(dat.Z) == min(dat.Z)
                            A = (max(dat.X)-min(dat.X))*(max(dat.Y)-min(dat.Y));
                        else
                            A = (max(dat.X)-min(dat.X))*(max(dat.Y)-min(dat.Y))*(max(dat.Z)-min(dat.Z));
                        end
                    else
                        % If the user defined a bonding surface
                        % exclude points outside the surface
                        if size(shp.Points, 2)==2
                            tf = inShape(shp,dat.X,dat.Y);
                        else
                            tf = inShape(shp,dat.X,dat.Y,dat.Z);
                        end
                        dat = dat(tf, :);
                    end
                end
                
                % Find the distanes between all cells 
                DIST = squareform(pdist(table2array(dat),'euclidean'));
                
                % use sort to get rid of the 0 diagonal elements 
                % (i.e. sum over j=1:n with j~=i)
                DIST = sort(DIST);
                n = size(dat, 1);
                const = A/(pi*(n*(n-1)));

                for d = BeginningDistance:DistanceIncrement:MaximumDistance
                    % weight factor equal to 1 if the distance between i and j is less
                    % than d and 0 otherwise
                    ksum = sum(sum(DIST(2:end, :) < d));
                    L_d = sqrt(const*ksum);
                    results = [results; {smpls{i}, pnms{j}, d, ksum, L_d, Groups(i)}];
                end
            end
        end

        %% Plot stuff
        fig1 = figure;
        fig1.Color = 'w';
        fig1.Visible = 'off';
        clf
        hold on
        title('Ripley"s K-Function')

        % add an export data option
        % Create Menu bar
        ExportMenu = uimenu(fig1);
        ExportMenu.Text = 'Export';
        % Create an export data option
        ExportPltDat = uimenu(ExportMenu);
        ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
        ExportPltDat.Text = 'Export Plot Data to .csv';
        
        fig2 = figure;
        fig2.Color = 'w';
        fig2.Visible = 'off';
        clf
        title('sum(kij)')
        hold on
        
        % add an export data option
        % Create Menu bar
        ExportMenu = uimenu(fig2);
        ExportMenu.Text = 'Export';
        % Create an export data option
        ExportPltDat = uimenu(ExportMenu);
        ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
        ExportPltDat.Text = 'Export Plot Data to .csv';
        
        if ~strcmp(popdiff, 'None')
            fig3 = figure;
            fig3.Color = 'w';
            fig3.Visible = 'off';
            clf
            title(['Difference from ' popdiff])
            hold on 
            
            % add an export data option
            % Create Menu bar
            ExportMenu = uimenu(fig3);
            ExportMenu.Text = 'Export';
            % Create an export data option
            ExportPltDat = uimenu(ExportMenu);
            ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
            ExportPltDat.Text = 'Export Plot Data to .csv';
            
        end

        clrs = jet(numel(pnms));
        % {'Sample', 'CellType', 'Distance Threshold', 'sum(kij)', 'L(d)', 'Group'};
        % strip the header
        results = results(2:end, :);
        GroupsU = unique(Groups, 'stable');
        for j = 1:numel(pnms)
            for k = 1:numel(GroupsU)
                logic = strcmp(results(:,2),pnms{j});
                logic2 = ([results{:,6}]==GroupsU(k))';
                logic = logical(logic.*logic2);
                
                figure(fig1.Number)
                X = [results{logic,3}];
                Y = [results{logic,5}];
                plot(X, Y, '.','Color', clrs(j, :), 'DisplayName', ['Group' num2str(GroupsU(k)) ' ' pnms{j}])
                % find the unique x values
                Xu = unique([results{logic,3}]);
                Ymean = zeros(numel(Xu), 1);
                % averaerage at each of those points
                for i_m = 1:numel(Xu) 
                    Ymean(i_m) = mean(Y(X==Xu(i_m)));
                end
                plot(Xu, Ymean, '-','Color', clrs(j, :), 'DisplayName', ['Group' num2str(GroupsU(k)) ' Avg ' pnms{j}])

                ax = gca;
%                 ax.YLim = [0 1000];
                ax.XLim = [0, MaximumDistance];
                ax.XLabel.String = 'Interaction Radius, um';
                ax.YLabel.String = 'Ripley"s K value';

                figure(fig2.Number)
                plot(zeros(sum(logic),1)+j, [results{logic,4}], 'o','Color', clrs(j, :), 'DisplayName', ['Group' num2str(GroupsU(k)) ' ' pnms{j}])
                ax = gca;
                ax.XTick = 1:1:numel(pnms);
                ax.XTickLabel = pnms;
            %     ax.YLim([0 2])
                ax.XLim = [0, numel(pnms)+1];
                ax.XTickLabelRotation = 30;
                ax.TickLabelInterpreter  = 'none';
                ax.YLabel.String = 'sum(kij)';

                if ~strcmp(popdiff, 'None')
                    logicd = strcmp(results(:,2),popdiff);
                    logic2 = ([results{:,6}]==GroupsU(k))';
                    logicd = logical(logicd.*logic2);
                
                    figure(fig3.Number)
                    Yd = [results{logicd,5}];
                    plot(X, Y-Yd, '.','Color', clrs(j, :), 'DisplayName', ['Group' num2str(GroupsU(k)) ' ' pnms{j}])
                    % find the unique x values
                    Ymean = zeros(numel(Xu), 1);
                    % averaerage at each of those points
                    for i_m = 1:numel(Xu) 
                        Ymeand(i_m) = mean(Y(X==Xu(i_m)) - Yd(X==Xu(i_m)));
                    end
                    plot(Xu, Ymeand, '-','Color', clrs(j, :), 'DisplayName', ['Group' num2str(GroupsU(k)) ' Avg ' pnms{j}])

                    ax = gca;
%                     ax.YLim = [-100 1000];
                    ax.XLim = [0, MaximumDistance];
                    ax.XLabel.String = 'Interaction Radius, um';
                    ax.YLabel.String = ['Difference from ' popdiff ' Ripley Values'];
                end
            end
        end
        fig1.Visible = 'on';
        fig2.Visible = 'on';
        fig3.Visible = 'on';
        figure(fig1.Number)

    end
end