function Cell_Density_Model_V1(app)
    %% Build a user interface
    % Makes bar graphs and gives a convienyent way to export data from
    % them
    %
    % Input:
    %   - app - Instance of CytoMAP
    %
    
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            if ~Helper.any_sample(app)
                return;
            end

            %% Build options for sorting
            [MFI_TYPE, sample] = Helper.find_MFI(app);
            if ~iscell(sample)
                sample = {sample};
            end
            tdata = Helper.populate_table(app, 'smpls', sample, 'MFI', 'AllCells', 'fill_checkbox', false);
            tdata((end-2):end, 5) = {false};
            tdata([tdata{:,1}]==1, 1) =   {1}; 
            tdata(end + 1, :) = {[], [], [], 1, false, 'Number of Cells', [], []}; 
            tdata(end + 1, :) = {[], false, 'Select All', [], false, 'Select All', false, 'Select All'}; 

            % Build the clustering options menus
            UIfig = uifigure('Name', 'Model Y Axis Density Options', 'Scrollable', 'on');
            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 1000 800];

            % Select Data type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Data Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10 45 400 15];
            DataType = uidropdown(UIfig);
            DataType.Items = {'Raster Scanned Neighborhoods', 'Cell Centered Neighborhoods',  'Individual Cells'};
            DataType.Value = DataType.Items{3};
            DataType.Position = alpha*[10 10 175 30];
            Helper.func_SetCLR(app, DataType, 'button')
            
            % Select X axis
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select X-Axis:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+10 45 400 15];
            XAxisDat = uidropdown(UIfig);
            XAxisDat.Items = tdata(:, 6);
            XAxisDat.Value = XAxisDat.Items{1};
            XAxisDat.Position = alpha*[10+175+10 10 175 30];
            Helper.func_SetCLR(app, XAxisDat, 'button')

            % Create the table of options
            t = uitable(UIfig);
            t.Data = tdata;
            t.Position = alpha*[0 75 1000 725];
            t.ColumnName = {'Group', 'Include','Phenotype (Must be in all samples)', ...
                            'Weight', 'Include', 'Y-Axis','Include', 'Sample'};
            t.ColumnEditable = [true true false true true false true false];
            t.ColumnWidth = {alpha*50, alpha*50, alpha*350, alpha*50, alpha*50, alpha*125, alpha*50, alpha*125};
            t.CellEditCallback = @(dd, p) switched_sample(app, DataType, XAxisDat, dd, p);

            DataPrep = uidropdown(UIfig);
            DataPrep.Items = Constants.neigh_cell_norm;
            DataPrep.Value = DataPrep.Items(3);
            DataPrep.Position = alpha*[10+175+175+10+75 40 175 30];
            Helper.func_SetCLR(app, DataPrep, 'button')
            DataPrep.Visible = 'off';
                       
            % Select Data Preperation for MFI
            DataPrepMFI = uidropdown(UIfig);
            DataPrepMFI.Items = Constants.cell_mfi_norm;
            DataPrepMFI.Value = DataPrepMFI.Items(1);
            DataPrepMFI.Position = alpha*[10+175+175+10+75 10 175 30];
            Helper.func_SetCLR(app, DataPrepMFI, 'button')

            % Create a Norm per sample or per dataset label
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Normalize per:');
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+10 50+5 75 20];
            % Create a Norm per sample or per dataset button
            NormPer = uibuttongroup(UIfig, 'Visible','off');
            NormPer.Position = alpha*[10+175+175+10, 10, 75, 45];
            Helper.func_SetCLR(app, NormPer, 'UICpopup')
            % Create two radio buttons in the button group.
            r1 = uiradiobutton(NormPer, 'Text','Sample',...
                'Position',alpha*[5, 25, 75, 20]);
            Helper.func_SetCLR(app, r1, 'table')
            r2 = uiradiobutton(NormPer, 'Text','Dataset',...
                'Position',alpha*[5, 5, 75, 20]);
            Helper.func_SetCLR(app, r2, 'table')
            NormPer.Visible = 'on';
            
            % Create Bin Options
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Number of Bins:');
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+10+75+175 50+5 100 20];
            NBins = uieditfield(UIfig,'numeric');
            NBins.Value = 20;
            NBins.Position = alpha*[10+175+175+10+75+175 40 75 20];
            
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Bin Length:');
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+10+75+175 25 75 20];
            LBin = uieditfield(UIfig,'numeric');
            LBin.Value = 5;
            LBin.Position = alpha*[10+175+175+10+75+175 10 75 20];
            
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Bin Start:');
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+10+75+175+100 50+5 100 20];
            RangeMin = uieditfield(UIfig,'numeric');
            RangeMin.Value = 0;
            RangeMin.Position = alpha*[10+175+175+10+75+175+100 40 75 20];
            
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Bin End:');
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+10+75+175+100 25 100 20];
            RangeMax = uieditfield(UIfig,'numeric');
            RangeMax.Value = 100;
            RangeMax.Position = alpha*[10+175+175+10+75+175+100 10 75 20];
            
            NBins.ValueChangedFcn = @(~,~) func_bins(LBin, NBins, RangeMin, RangeMax, 'NBins');
            LBin.ValueChangedFcn = @(~,~) func_bins(LBin, NBins, RangeMin, RangeMax, 'LBin');
            RangeMin.ValueChangedFcn = @(~,~) func_bins(LBin, NBins, RangeMin, RangeMax, 'RangeMin');
            RangeMax.ValueChangedFcn = @(~,~) func_bins(LBin, NBins, RangeMin, RangeMax, 'RangeMax');

            
            
            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap(t, app, DataType, XAxisDat, DataPrep, DataPrepMFI, NormPer, LBin, NBins, RangeMin, RangeMax));
            btn.Position = alpha*[1000-100, 5, 100, 50];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')

            DataType.ValueChangedFcn = @(~, ~) ChangeDataType(app, t, DataType, XAxisDat, DataPrep, DataPrepMFI);

%             ChangeDataType(app, t, DataType, DataPrep, DataPrepMFI)
            function func_bins(LBin, NBins, RangeMin, RangeMax, tag)
             switch tag
                 case 'NBins'
                     NBins.Value = round(NBins.Value);
                     if NBins.Value <= 0
                         NBins.Value = 1;
                     end
                     LBin.Value = (RangeMax.Value-RangeMin.Value)/NBins.Value;

                 case 'LBin'
                      if LBin.Value <= 0
                         LBin.Value = 1;
                     end
                     NBins.Value = round((RangeMax.Value-RangeMin.Value)/LBin.Value);
                     func_bins(LBin, NBins, RangeMin, RangeMax, 'NBins')
                     
                 case 'RangeMin'
                      if RangeMin.Value >= RangeMax.Value
                         RangeMin.Value = RangeMax.Value-1;
                      end
                     func_bins(LBin, NBins, RangeMin, RangeMax, 'NBins')
                 case 'RangeMax'
                     if RangeMax.Value <= RangeMin.Value
                         RangeMax.Value = RangeMin.Value+1;
                     end
                     func_bins(LBin, NBins, RangeMin, RangeMax, 'NBins')
             end
            end
        
            function ChangeDataType(app, t, DataType, XAxisDat, DataPrep, DataPrepMFI)
                % changes visiblity of manual choice numerical input.
                p.Indices = [7 7];
                p.EditData = 0;
                switched_sample(app, DataType, XAxisDat, t, p);
                if strcmp(DataType.Value, 'Individual Cells')
                    DataPrep.Visible = 'off';
                    DataPrepMFI.Items = Constants.cell_mfi_norm;
                    DataPrepMFI.Value = DataPrepMFI.Items(1);
                else
                    DataPrep.Visible = 'on';
                    DataPrepMFI.Items = Constants.neigh_mfi_norm;
                    DataPrepMFI.Value = DataPrepMFI.Items(3);
                end
            end

            function switched_sample(app, DataType, XAxisDat, dd, p)
                if p.Indices(2) == 7
                    % Make sure that at least one thing is selected
                    if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 7)), 7)))
                        dd.Data{p.Indices(1), 7} = true;
                        return;
                    end
                    % Process Select All button
                    fill_select_all = false;
                    if p.Indices(1) == find(strcmp(dd.Data(:, 8), 'Select All')) && p.Indices(2) == 7
                        ind = ~cellfun('isempty', dd.Data(:, 7));
                        dd.Data(ind, 7) = {logical(p.NewData)};
                        if ~logical(p.NewData)
                            dd.Data{strcmp(dd.Data(:, 8), app.DataN.Value), 7} = true;
                        end
                        fill_select_all = logical(p.NewData);
                    end
                    if strcmp(DataType.Value, 'Individual Cells')
                        scan_type = 'AllCells';
                    elseif strcmp(DataType.Value, 'Cell Centered Neighborhoods')
                        scan_type = 'MFICCN';
                    elseif strcmp(DataType.Value, 'Raster Scanned Neighborhoods')
                        scan_type = 'MFIRSN';
                    end

                    ind = ~cellfun('isempty', dd.Data(:, 7));
                    ind(end) = false;
                    ind(ind) = logical(cell2mat(dd.Data(ind, 7)));
                    smpls = dd.Data(ind, 8);
                    if isempty(smpls)
                        smpls = {app.DataN.Value};
                    end

                    tdataTMP = Helper.populate_table(app, ...
                        'smpls', smpls, ...
                        'MFI', scan_type, ...
                        'prev_table', dd.Data ...
                    );

                    tdataTMP((end-2):end,5) = {false};
                    
                    switch scan_type
                        case 'AllCells'
                            tdataTMP(end + 1, :) =  {[], [], [], 1, false, 'Number of Objects', [], []};  
                        case 'MFICCN'
                            tdataTMP(end + 1, :) =  {[], [], [], 1, false, 'Number of Neighborhoods', [], []};  
                        case 'MFIRSN'
                            tdataTMP(end + 1, :) =  {[], [], [], 1, false, 'Number of Neighborhoods', [], []};  
                    end
                    
                    tdataTMP(end + 1, :) = {[], false, 'Select All', [], false, 'Select All', fill_select_all, 'Select All'}; 
                    dd.Data = tdataTMP;
                    XAxisDat.Items = tdataTMP(:, 6);
                    
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 3), 'Select All')) && p.Indices(2) == 2
                    ind = ~cellfun('isempty', dd.Data(:, 3));
                    dd.Data(ind, 2) = {logical(p.NewData)};
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 6), 'Select All')) && p.Indices(2) == 5
                    ind = ~cellfun('isempty', dd.Data(:, 5));
                    dd.Data(ind, 5) = {logical(p.NewData)};
                end
            end

            function backend_wrap(optionstable, app, DataType, XAxisDat, DataPrep, DataPrepMFI, NormPer, LBin, NBins, RangeMin, RangeMax)
                tmp_data = optionstable.Data(1:end - 1, :);
                INDSmpls = [tmp_data{1:numel(app.DataN.Items), 7}]==1;
                SortNames = tmp_data(1:numel([tmp_data{:,2}]),3);
                SortNames = SortNames([tmp_data{:,2}]==1,:);
                MFIList = tmp_data(1:numel([tmp_data{:,5}]),6);
                MFIList = MFIList([tmp_data{:,5}]==1,:);

                groups = [tmp_data{[tmp_data{:,2}]==1,1}];
                if numel(unique(groups))==1
                    groups = 1:numel(groups);
                end
                
                func_backend( ...
                    app, ...
                    tmp_data(INDSmpls, 8), ... Samples
                    SortNames, ... Phenotypes
                    groups, ... Phenotype Groups
                    MFIList, ... MFIs
                    [tmp_data{[tmp_data{:,5}]==1,4}], ... MFI weights
                    DataType.Value, ...
                    XAxisDat.Value, ...
                    DataPrep.Value, ...
                    DataPrepMFI.Value, ...
                    NormPer.SelectedObject.Text, ...
                    LBin.Value,  ...
                    NBins.Value, ...
                    RangeMin.Value, ...
                    RangeMax.Value ...
                );
            end
        
    %% Run the Ripley K function
    function func_backend(app, smplnms, SortNames, pheno_groups, MFIList, MFI_weights, DataType, XAxisDat, DataPrep, DataPrepMFI, NormOpt, LBin, NBins, RangeMin, RangeMax)
        %% Load the data
        NPlots = numel(MFIList);
        
        if ~isempty(MFIList)
            MFI_weights = [MFI_weights, 1];
            MFILoad = Helper.valid_channel([MFIList; XAxisDat]);
            INDnum = contains(MFIList, 'Number of Objects');
            if any(INDnum)
                INDnum = find(INDnum);
                MFILoad(INDnum) = [];
                MFI_weights(INDnum) = [];
            elseif any(contains(MFIList, 'Number of Neighborhoods'))
                INDnum = find(contains(MFIList, 'Number of Neighborhoods'));
                MFILoad(INDnum) = [];
                MFI_weights(INDnum) = [];
            else
                INDnum = [];
            end
            
        end

        
        loopn = 0;
        fignum = cell(NPlots, 1);
        groups = unique(pheno_groups, 'stable');
        CellNames = cell(numel(groups), 1);

        Dat_PlotXBox = cell((max(groups)-min(groups)+1), NPlots);
        Dat_PlotYBox = cell((max(groups)-min(groups)+1), NPlots);  
        
        for groupi = 1:numel(groups)
            
            % Pull the cell names belonging to this group
            GroupNames = SortNames(pheno_groups==groups(groupi));
            if numel(GroupNames)==1
                CellNames(groups(groupi)) = GroupNames;
            else
                CellNames{groups(groupi)} = strjoin(GroupNames, ' and ');
            end
            % Not Individual Cells
            if ~strcmp(DataType, 'Individual Cells')
                [Dat_Pre, ~, INDDatPre, INDzrs, INDons] = Helper.func_loaddata(app, ...
                    smplnms, ...
                    GroupNames, ...
                    MFILoad, ...
                    [], ...
                    MFI_weights, ...
                    DataType, ...
                    DataPrep, ...
                    NormOpt, ...
                    1, ...
                    DataPrepMFI); %<Change back to 0 maybe
                if strcmp(DataType, 'Raster Scanned Neighborhoods')
                    type = 'RSN';
                elseif strcmp(DataType, 'Cell Centered Neighborhoods')
                    type = 'CCN';
                end
            elseif strcmp(DataType, 'Individual Cells')
                [Dat_Pre, ~, INDDatPre, INDzrs, INDons] = Helper.func_loaddata(app, ...
                    smplnms, ...
                    GroupNames, ...
                    MFILoad, ...
                    [], ...
                    MFI_weights, ...
                    DataType, ...
                    'Cellularity: Number of Cells / Neighborhood', ...
                    NormOpt, ...
                    1, ...
                    DataPrepMFI); 
                type = 'Cells';
                Dat_Pre = Dat_Pre(:, (numel(GroupNames)+1):end);
            end
            
            %% Calculate Stuff
            
            
           
            %% Plot stuff
            for MFI_i = 1:NPlots
                
                MFI_X = size(Dat_Pre,2);
                
                if loopn==0
                    
                    fig = figure;
                    fig.Name = MFIList{MFI_i};
                    fignum{MFI_i} =  fig.Number;
                    
                    % add an export data option
                    % Create Menu bar
                    ExportMenu = uimenu(fig);
                    ExportMenu.Text = 'Export';
                    % Create an export data option
                    ExportPltDat = uimenu(ExportMenu);
                    ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                    ExportPltDat.Text = 'Export Plot Data to .csv';

                    fig.Color = app.GUIOPTS.bgclr;
                    fig.InvertHardcopy = 'off';
                    
                else
                   fig = figure(fignum{MFI_i});
                end
                hold on
                        for smlp_i = 1:numel(smplnms)
                                    
                            Dat_PlotX = Dat_Pre(INDDatPre{smlp_i}{1}(1):INDDatPre{smlp_i}{1}(2), MFI_X);
                            [Dat_PlotX, INDsort] = sort(Dat_PlotX);

                            if isempty(INDnum) || MFI_i ~= INDnum
                                %% Case for MFI plot
                                Dat_PlotY = Dat_Pre(INDDatPre{smlp_i}{1}(1):INDDatPre{smlp_i}{1}(2), MFI_i);
                                % sort the data in order of x values
                                Dat_PlotY = Dat_PlotY(INDsort);
                                
                                plot(Dat_PlotX, Dat_PlotY, '.', 'DisplayName', [smplnms{smlp_i} '_' CellNames{groups(groupi)}])
                                
                                % Bin the data
                                BinEdge = RangeMin;
                                Dat_PlotX_Bined = zeros(NBins, 1);
                                Dat_PlotY_Bined = zeros(NBins, 1);
                                Dat_PlotY_Bined_Num = zeros(NBins, 1);
                                for Bin_i = 1:NBins
                                    % center position of this bin
                                    Dat_PlotX_Bined(Bin_i) = BinEdge+(LBin/2);
                                    Dat_PlotY_BinedIND = BinEdge <= Dat_PlotX & Dat_PlotX <= (BinEdge + LBin);
                                    Dat_PlotY_Bined(Bin_i)  = mean(Dat_PlotY(Dat_PlotY_BinedIND), 'omitnan');
                                    Dat_PlotY_Bined_Num(Bin_i)  = numel(Dat_PlotY(Dat_PlotY_BinedIND));
                                    BinEdge = BinEdge+LBin;
    % % %                                 plot(Dat_PlotX(Dat_PlotY_BinedIND), Dat_PlotY(Dat_PlotY_BinedIND), '.r', 'DisplayName', [smplnms{smlp_i} '_' CellNames{groups(groupi)}])
                                end
                                plot(Dat_PlotX_Bined, Dat_PlotY_Bined, '-o', 'DisplayName', [smplnms{smlp_i} '_' CellNames{groups(groupi)} '_Mean_MFI_Of_Binned'])
                                plot(Dat_PlotX_Bined, Dat_PlotY_Bined_Num, '--', 'Color', [0.8,0.8,0.8], 'DisplayName', [smplnms{smlp_i} '_' CellNames{groups(groupi)} '_Number_in_Bin'])
                                
                            elseif ~isempty(INDnum) || MFI_i == INDnum
                                %% Case for number of objects plot
                                % Bin the data
                                BinEdge = RangeMin;
                                Dat_PlotX_Bined = zeros(NBins, 1);
                                Dat_PlotY_Bined_Num = zeros(NBins, 1);
                                for Bin_i = 1:NBins
                                    % center position of this bin
                                    Dat_PlotX_Bined(Bin_i) = BinEdge+(LBin/2);
                                    Dat_PlotY_BinedIND = BinEdge <= Dat_PlotX & Dat_PlotX <= (BinEdge + LBin);
                                    Dat_PlotY_Bined_Num(Bin_i)  = sum(Dat_PlotY_BinedIND);
                                    BinEdge = BinEdge+LBin;
                                end
                                % Normalize the data
%                                 Dat_PlotY_Bined_Num = Dat_PlotY_Bined_Num./max(Dat_PlotY_Bined_Num);
                                
                                plot(Dat_PlotX_Bined, Dat_PlotY_Bined_Num, '.', 'DisplayName', [smplnms{smlp_i} '_' CellNames{groups(groupi)} '_Number_in_Bin'])
% %                                 plot(Dat_PlotX_Bined, Dat_PlotY_Bined_Num./max(Dat_PlotY_Bined_Num), '.', 'DisplayName', [smplnms{smlp_i} '_' CellNames{groups(groupi)} '_Number_in_Bin'])

                                %% fit with two poisson distributions
                                
% P(x) = ( (lam_a^x)*exp(-lam_a) ) / factorial(x) 
                                mdl = @(p0, p1, p2, p3, p4, p5, x) ...
                                    p1.*x.*(exp(-(x-p0).^2/p2^2)) + ...
                                    p3.*(exp(-(x-p4).^2/p5^2));
                                
                                mdSl = @(p0, p1, p2, x) ...
                                    p1.*x.*(exp(-(x-p0).^2/p2^2));
                                
                                mdS2 = @(p3, p4, p5, x) ...
                                    p3.*(exp(-(x-p4).^2/p5^2));
                                
                                [maxnum, maxIND] = max(Dat_PlotY_Bined_Num);
                                maxpos = Dat_PlotX_Bined(maxIND);
                                
                                startpts = [(maxpos-0.2*maxpos), maxnum, 5, maxnum-0.1*maxnum, (maxpos+0.5*maxpos), 50];
                                
% %                                 startpts = [2.0, 75, 4, 300, 53, 45];
                                
                                [f0, gof, output] = fit(Dat_PlotX_Bined, Dat_PlotY_Bined_Num , mdl,'StartPoint',startpts);
                                
% %                                 figure(100)
% %                                 hold on
% % %                                 gof
% % %                                 output
% %                                 plot(groups(groupi)-0.2, f0.p0, '.')
% %                                 plot(groups(groupi), f0.p4, '.')
% %                                 plot(groups(groupi)+0.2, mean(Dat_PlotX), '.')
                                
                                figure(fignum{MFI_i})
                                smoothedx = min(Dat_PlotX_Bined):((max(Dat_PlotX_Bined)-min(Dat_PlotX_Bined))/10000):max(Dat_PlotX_Bined);
                                
                                plot(smoothedx, f0(smoothedx), 'k-')
                                plot(smoothedx, mdSl(f0.p0,f0.p1,f0.p2,smoothedx), 'r:')
                                plot(smoothedx, mdS2(f0.p3,f0.p4,f0.p5,smoothedx), 'm:')
%                                 plot(Dat_PlotX_Bined, f0(Dat_PlotX_Bined), 'k-')
                                
                                
                            end
                            
                        end % end loop through samples
                
                % Label Stuff
                if loopn==0
                    ax = gca;
                    ax.YLabel.String = MFIList{MFI_i};
                    ax.XLabel.String = XAxisDat;
                    if isempty(INDnum) || MFI_i ~= INDnum
                        ax.YLim = [0 max(Dat_PlotY)];
                    else
                        ax.YLim = [0 max(Dat_PlotY_Bined_Num)];
                    end
                    
                end
                
                if isempty(INDnum) || MFI_i ~= INDnum
                    ax.YLim = [0 max(max(Dat_PlotY), ax.YLim(2))];
                else
                    ax.YLim = [0 max(max(Dat_PlotY_Bined_Num), ax.YLim(2))];
                end
                
            end % end plot loop
            loopn = loopn+1;  
        end % end loop through groups            
        
        for MFI_i = 1:numel(MFIList)
            fig = figure(fignum{MFI_i});
            
            box off
        end
        
    end
end