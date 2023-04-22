function Nearest_Neighbor_Analysis(app)
    % Nearest neighbor function finds the probablility of cells neighboring
    % each other It is useful to understand how a
    % population/region relates to global average, or other
    % populations.
    %
    % Input:
    %   - app - Instance of CytoMAP.

    alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

    if ~Helper.any_sample(app)
        return;
    end
    
    if ~Helper.any_net(app)
        return;        
    end

    % Build the clustering options menus
    UIfig = uifigure('Name', 'Nearest Neighbor Calculation', 'Scrollable', 'on');
    Helper.func_SetCLR(app, UIfig, 'UIfigure');
    UIfig.Position = alpha*[10 10 1365 800];

    % Initialize the table of options
    t = uitable(UIfig);
    t.Position = alpha*[0 70 1360 730];

    % Select Neighborhood type
    lbl = uilabel(UIfig); lbl.Text = sprintf('Select Model:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = alpha*[10 45 175 15];
    DataType = uidropdown(UIfig);
    DataType.Items = fieldnames(app.net);
    if ~isempty(DataType.Items)
        DataType.Value = DataType.Items{1};
    end
    DataType.Position = alpha*[10 10 175 30];
    Helper.func_SetCLR(app, DataType, 'button')
    p = struct('Indices', [], 'EditData', 1, 'NewData', []);
    DataType.ValueChangedFcn = @(~, ~) switched_sample(app, DataType, t, p);

    %Build the table
    [~, sample] = Helper.find_MFI(app);
    tmpData = Helper.populate_table(app, 'smpls', {sample}, 'fill_checkbox', false);
    
    % Initialize table
    tData = cell(max(app.net.(DataType.Value).NR+1, size(tmpData, 1)) + 1, 10);
    % Cells
    tData(1:numel(tmpData(:, 2)), 1)= tmpData(:, 2);
    tData(1:sum(~cellfun(@isempty,tmpData(:, 3))), 2) = {1};
    tData(1:numel(tmpData(:, 3)), 3) = tmpData(:, 3);
    tData(end, [1,3]) = {false, 'Select All'};
    % Sample
    tData(1:numel(tmpData(:, 4)), 4) = tmpData(:, 4);
    tData(1:sum(~cellfun(@isempty,tmpData(:, 5))), 5) = {1};
    tData(1:numel(tmpData(:, 5)), 6) = tmpData(:, 5);
    tData(end, [4, 6]) = {false, 'Select All'};
    % Regions
    tData(1, 7) = {false};
    tData(2:app.net.(DataType.Value).NR+1, 7) = {true};
    tData(2:app.net.(DataType.Value).NR+1, 8) = {1};
    tData(1:app.net.(DataType.Value).NR+1, 9) = num2cell(0:app.net.(DataType.Value).NR);

    % Populate the table of options
    t.Data = tData;
    t.ColumnName = { ...
        'Include', 'Group', 'Phenotype', ...
        'Include', 'Group', 'Sample', ...
        'Include', 'Group', 'Region' ...
    };
    t.ColumnEditable = [true false false true true false true true false];
    t.ColumnWidth = {alpha*50,alpha*50, alpha*350, alpha*50, alpha*50, alpha*225, alpha*50, alpha*50, alpha*100};
    t.CellEditCallback = @(dd, p) switched_sample(app, DataType, dd, p);

    % Select Heatmap type
    lbl = uilabel(UIfig); lbl.Text = sprintf('Select Heatmap Type:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = alpha*[10+175+175+75+10 45 400 15];
    MAPType = uidropdown(UIfig);
    MAPType.Items = {'Individual Heatmap for each Sample', ...
                      'Combined Heatmap of all Samples'};
    MAPType.Value = {'Individual Heatmap for each Sample'};
    MAPType.Position = alpha*[10+175+175+75+10 10 175 30];
    Helper.func_SetCLR(app, MAPType, 'button')

    % Select Color Scale
    lbl = uilabel(UIfig); lbl.Text = sprintf('Color Scale:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = alpha*[10+175+175+175+75+10 45 100 15];
    MAPScale = uidropdown(UIfig);
    MAPScale.Items = {'linear', ...
                      'log'};
    MAPScale.Value = {'linear'};
    MAPScale.Position = alpha*[10+175+175+175+75+10 10 100 30];
    Helper.func_SetCLR(app, MAPScale, 'button')

    % Select Calculation Type
    lbl = uilabel(UIfig); lbl.Text = sprintf('Calculation:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = alpha*[10+175+175+175+100+75+10 45 100 15];
    CalcType = uidropdown(UIfig);
    CalcType.Items = {'Nearest-Neighbor Prevalence', ...
                      ' '};
    CalcType.Value = CalcType.Items{1};
    CalcType.Position = alpha*[10+175+175+175+100+75+10 10 175 30];
    Helper.func_SetCLR(app, CalcType, 'button')

    % Select Data transform 
    lbl = uilabel(UIfig); lbl.Text = sprintf('Neighbors:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = alpha*[10+175+175+175+100+175+75+10 45 100 15];
    NeighType = uidropdown(UIfig);
    NeighType.Items = {'Distance', ...
                      'Number'};
    NeighType.Value = {'Number'};
    NeighType.Position = alpha*[10+175+175+175+100+175+75+10 10 100 30];
    Helper.func_SetCLR(app, lbl, 'button')

    % Select Confidence Interval
    lbl = uilabel(UIfig); lbl.Text = sprintf('Neighbors within:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = alpha*[10+175+175+175+100+175+75+10+100 45 100 15];
    NNeigh = uieditfield(UIfig, 'numeric' ...
    );
    NNeigh.Value = 10;
    NNeigh.Position = alpha*[10+175+175+175+100+175+75+10+100 10 100 30];
    NNeigh.Editable = 'on';
    Helper.func_SetCLR(app, NNeigh, 'UIField')

    % Create a plot types button group
    plottypes = uibuttongroup(UIfig, 'Visible','off');
    plottypes.Position = [10+175, 5, 175, 65];
    Helper.func_SetCLR(app, plottypes, 'UICpopup')
    % Create two radio buttons in the button group.
    r1 = uicheckbox(plottypes, 'Text','Percentage Heatmap',...
        'Position',[5, 2.5, 150, 20], 'Value', true);
    Helper.func_SetCLR(app, r1, 'table')
    r2 = uicheckbox(plottypes, 'Text','Deviation Heatmap',...
        'Position',[5, 2.5+20, 150, 20], 'Value', true);
    Helper.func_SetCLR(app, r2, 'table')
    r3 = uicheckbox(plottypes, 'Text','Stacked Bar Graph',...
        'Position',[5, 2.5+2*20, 150, 20], 'Value', true);
    Helper.func_SetCLR(app, r3, 'table')
%     r4 = uicheckbox(plottypes, 'Text','Fold Change',...
%         'Position',[5+150, 22.5, 150, 20], 'Value', true);
%     Helper.func_SetCLR(app, r4, 'table')
%     r5 = uicheckbox(plottypes, 'Text','Normalized Mean',...
%         'Position',[5+150+150, 2.5, 150, 20], 'Value', true);
%     Helper.func_SetCLR(app, r5, 'table')
    plottypes.Visible = 'on';
    
    % Create a push button
    btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap( ...
        t, app, DataType, MAPType, MAPScale, CalcType, NeighType, NNeigh, plottypes ...
    ));
    btn.Position = alpha*[1365-100, 5, 100, 50];
    btn.Text = 'Ok';
    Helper.func_SetCLR(app, btn, 'button')

    function switched_sample(app, DataType, dd, p)
        % if the samples index was changed
        if isempty(p.Indices) || p.Indices(2) == 4
            % Make sure that at least one thing is selected
            if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, p.Indices(2))), p.Indices(2))))
                dd.Data{p.Indices(1), p.Indices(2)} = true;
                return;
            end
            % Process Select All button
            fill_select_all = false;
            if ~isempty(p.Indices) && p.Indices(1) == find(strcmp(dd.Data(:, 6), 'Select All'))
                ind = ~cellfun('isempty', dd.Data(:, p.Indices(2)));
                dd.Data(ind, p.Indices(2)) = {logical(p.NewData)};
                if ~logical(p.NewData)
                    dd.Data{strcmp(dd.Data(:, 6), app.DataN.Value), p.Indices(2)} = true;
                end
                fill_select_all = logical(p.NewData);
            end
            ind = ~cellfun('isempty', dd.Data(:, 4));
            ind(end) = false;
            ind(ind) = logical(cell2mat(dd.Data(ind, 4)));

            %Find which samples were selected
            fcsz = ~cellfun('isempty', dd.Data(:, 3));
            fcsz(end) = false;
            tmpDat = cell(size(dd.Data, 1) - 1, 5);
            tmpDat(fcsz, 1) = {1};
            tmpDat(fcsz, 2) = dd.Data(fcsz, 1);
            tmpDat(:, 3) = dd.Data(1:end - 1, 3);
            
            tmpDat(:, 4) = dd.Data(1:end - 1, 4);
            tmpDat(:, 5) = dd.Data(1:end - 1, 6);

            smpls = dd.Data(ind, 7);
            if isempty(smpls)
                smpls = {app.DataN.Value};
            end

            tmpDat = Helper.populate_table(app, ...
                'smpls', smpls, ...
                'prev_table', tmpDat ...
            );

            % Keep column which is not conserved in populate table
            keep_6_col = dd.Data(:, 5);
            prev_smpls = dd.Data(:, 6);
            keep_regions = dd.Data(:, 7:end);

            dd_ph = ~cellfun('isempty', dd.Data(:, 1));
            dd_ph(end) = false;
            dd_ph = dd.Data(dd_ph, 3);
            tmp_ph = ~cellfun('isempty', tmpDat(:, 2));
            tmp_ph = tmpDat(tmp_ph, 3);
            if isempty(p.Indices) && app.net.(DataType.Value).NR+1 < size(dd.Data, 1)
                dd.Data(:, 7:end) = [];
                dd.Data(1, 7) = {false};
                dd.Data(2:app.net.(DataType.Value).NR+1, 7) = {true};
                dd.Data(2:app.net.(DataType.Value).NR+1, 8) = {1};
                dd.Data(1:app.net.(DataType.Value).NR+1, 9) = num2cell(0:app.net.(DataType.Value).NR);
            elseif ~Helper.setequal(dd_ph, tmp_ph) || isempty(p.Indices)
                dd.Data = cell(max(app.net.(DataType.Value).NR + 1, size(tmpDat, 1)) + 1, 9);
                % Cells
                fcsz = ~cellfun('isempty', tmpDat(:, 3));
                logic_scale = false(size(dd.Data, 1), 1);
                logic_scale(1:size(fcsz, 1)) = fcsz;
                dd.Data(logic_scale, 1) = tmpDat(fcsz, 2);
                dd.Data(1:numel(tmpDat(:, 3)), 2) = tmpDat(:, 3);
                dd.Data(end, 1:2) = {false, 'Select All'};
                % Sample
                dd.Data(1:numel(tmpDat(:, 4)), 5) = tmpDat(:, 4);
                if sum(~cellfun('isempty', keep_6_col)) == sum(~cellfun('isempty', tmpDat(:, 5)))
                    if numel(keep_6_col) > size(dd.Data, 1)
                        dd.Data(:, 6) = keep_6_col(1:size(dd.Data, 1));
                    else
                        dd.Data(1:numel(keep_6_col), 6) = keep_6_col;
                    end
                else
                    prev_idx = ~cellfun('isempty', keep_6_col);
                    new_idx = ~cellfun('isempty', dd.Data(:, 7));

                    prev_smpls = prev_smpls(prev_idx);
                    keep_6_col = keep_6_col(prev_idx);
                    new_smpls  = dd.Data(new_idx, 7);
                    for new_i=1:numel(new_smpls)
                        find_idx = find(strcmp(prev_smpls, new_smpls(new_i)));
                        if isempty(find_idx)
                            dd.Data(new_i, 6) = keep_6_col(find_idx);
                        else
                            dd.Data(new_i, 6) = {1};
                        end
                    end
                end
                dd.Data(1:numel(tmpDat(:, 5)), 6) = tmpDat(:, 5);
                dd.Data(end, [4, 6]) = {fill_select_all, 'Select All'};
                % Regions
                if ~isempty(p.Indices)
                    if size(keep_regions, 1) > size(dd.Data, 1)
                        dd.Data(:, 7:end) = keep_regions(1:size(dd.Data,1), :);
                    else
                        dd.Data(1:size(keep_regions, 1), 7:end) = keep_regions;
                    end
                else
                    dd.Data(1,7) = {false};
                    dd.Data(2:app.net.(DataType.Value).NR+1, 7) = {true};
                    dd.Data(2:app.net.(DataType.Value).NR+1, 8) = {1};
                    dd.Data(1:app.net.(DataType.Value).NR+1, 9) = num2cell(0:app.net.(DataType.Value).NR);
                end
            end
        elseif p.Indices(1) == find(strcmp(dd.Data(:, 3), 'Select All')) && p.Indices(2) == 1
            ind = ~cellfun('isempty', dd.Data(:, p.Indices(2)));
            dd.Data(ind, p.Indices(2)) = {logical(p.NewData)};
        end
    end

    function backend_wrap(t, app, DataType, MAPType, MAPScale, CalcType, NeighType, NNeigh, plottypes)
        % remove the select all button
        tmp_data = t.Data(1:end - 1, :);
        backend_Calculations_func( ...
            app, ...
            cell2mat(tmp_data([tmp_data{:, 7}]==1, 8)'), ... Region Groups
            cell2mat(tmp_data([tmp_data{:, 7}]==1, 9)'), ... NRegions
            tmp_data([tmp_data{:, 4}]'==1, 6), ...           Samples
            cell2mat(tmp_data([tmp_data{:, 4}]==1, 5)'), ... Sample Groups
            tmp_data([tmp_data{:, 1}]'==1, 3), ...           Phenotypes
            tmp_data([tmp_data{:, 1}]'==1, 2), ...           Phenotype Groups
            DataType.Value, ...
            MAPType.Value, ...
            MAPScale.Value, ...
            CalcType.Value, ...
            NeighType.Value, ...
            NNeigh.Value, ...
            plottypes ...
        );
    end % end of backend wrapper

    function backend_Calculations_func(app, reg_groups, NRegions, smplnms, smp_groups, phenos, phe_groups, DataType, MAPType, MAPScale, CalcType, NeighType, NNeigh, plottypes)
        % CELLINTERACTION_BACKEND Back-End of Cell Interaction function. Allows user to plot
        % different statistics as different variations of a heatmap in
        % between different cell populations/cell regions in different
        % samples/groups. It is useful to understand how a
        % population/region relates to global average, or other
        % populations.
        %
        % Input:
        %   - app - Instance of CytoMAP.
        %   - groups - array - Groupings of regions. All of the regions
        %       in the same group will be plotted on the same figure.
        %   - NRegions - array - All of the possible regions.
        %   - smplnms - cell - Samples from which data to make heatmap
        %       will be pulled.
        %   - smp_groups - array - Groupings of samples. Samples in the
        %       same group will be included on the same plot, if the
        %       Combine Heatmap of all Samples is chosen.
        %   - phenos - cell - Phenotypes on which heatmap will be
        %       defined (alongside with MFIs).
        %   - MFIList - cell - MFIs on which heatmap will be defined
        %       (alongside with Phenotypes).
        %   - DataType - string - Name of the Model to use in creating
        %       the heatmap. Has to exist as field under app.net.
        %   - DataPrep - string - How to normalize and prepare the
        %       data. Same as in Helper.Func_DataPrep.
        %   - MAPType - string - Whether to make one `Combined Heatmap
        %       of all Samples` or multiple `Individual Heatmap for
        %       each Sample`.
        %   - MAPScale - string - Whether to apply apply `log` or
        %       `linear` scale to heatmap in the end.
        %   - CalcType - string - What sort of heatmap to create.
        %       Current Options are:
        %           - Pearson Correlation Coefficient
        %           - Spearman Correlation Coefficient
        %           - Kendall Correlation Coefficient
        %           - Covariance
        %           - Cross-correlation
        %           - Conditional Phenotype/Channel value given region
        %           - Conditional Phenotype/Channel value given 
        %               Phenotype/Channel
        %   - DatScale - string - Whether to apply apply `log` or
        %       `linear` scale to data before processing in the end.
        %   - ConfInt - numeric in range (0, 1] - Whether to build
        %       Confidence interval or not, and if so, then how big
        %       they should be.

        % Make a waitbar
        vPD = waitbar(0, 'Initializing', ...
            'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));

        % Pull which plot types were selected
        values = [plottypes.Children.Value];
% %             % Plot Stacked Bargraphs?
% %             values(1);
% %             % Deviation Heatmap?
% %             values(2);
% %             % Percentage Heatmap?
% %             values(3);

        
        %% Load the data

        MFIList = [{'X'}, {'Y'}, {'Z'}];

        [Dat_PreMAIN, ~, INDDatPre, INDzrs, INDons] = Helper.func_loaddata(app, ...
            smplnms, ...
            phenos, ...
            MFIList, ...
            [], ...
            [], ...
            'Individual Cells', ...
            'Cellularity: Number of Cells / Neighborhood', ...
            'Sample', ...
            1, ...
            'Cellularity: Number of Cells / Neighborhood'); 
                
        %% Calculate stuff
        smp_Ngroups = numel(unique(smp_groups));
        for Sgroup_i = 1:smp_Ngroups
            smplnms_i = smplnms(ismember(smp_groups,Sgroup_i));
            smplsIND = find(ismember(smp_groups,Sgroup_i));
            for smpl_i=1:numel(smplnms_i)
                %% Prepare the data the data
                if ~isvalid(vPD)
                    return;
                end
                waitbar(0.1*(smpl_i/numel(smplnms_i)),vPD, 'Loading Data!');

                CellNms = phenos;

                % Load the Region data
                ROW2 = app.data.(smplnms_i{smpl_i}).AllCells.([Constants.other_tag, DataType]);
                ROW2 = ROW2(INDons{smpl_i});

                %% Combine all data Samples from specific groups if selected
                switch MAPType
                    case 'Individual Heatmap for each Sample'
                        Dat_Pre = Dat_PreMAIN(INDDatPre{smplsIND(smpl_i)}{1}(1):INDDatPre{smplsIND(smpl_i)}{1}(2),:);
                    case 'Combined Heatmap of all Samples'
                        if smpl_i==1
                            CombinedDat = Dat_PreMAIN(INDDatPre{smplsIND(smpl_i)}{1}(1):INDDatPre{smplsIND(smpl_i)}{1}(2),:);
                            CombinedROW2 = ROW2;
                        else
                            CombinedDat = [CombinedDat; Dat_PreMAIN(INDDatPre{smplsIND(smpl_i)}{1}(1):INDDatPre{smplsIND(smpl_i)}{1}(2),:)];
                            CombinedROW2 = [CombinedROW2; ROW2];
                        end
                        if smpl_i==numel(smplnms_i)
                            Dat_Pre = CombinedDat;
                            ROW2 = CombinedROW2;
                        end
                end

                % Only include selected regions
                Ngroups = numel(unique(reg_groups));

                %% Do the Calculations if it is the last sample or you
                % are making individual plots
                if strcmp(MAPType, 'Individual Heatmap for each Sample') || smpl_i==numel(smplnms_i)
                    
                    for group_i = 1:Ngroups
                        
                        waitbar(0.1+0.5*(smpl_i/numel(smplnms_i)),vPD, 'PROCESSING!');
                        
                        NRegions_i = NRegions(ismember(reg_groups,group_i));
                        pcst_i = Dat_Pre(ismember(ROW2, NRegions_i), :);
                        ROW2_i = ROW2(ismember(ROW2, NRegions_i), :);
                        YLabel = CellNms;
                        switch CalcType
                            case 'Nearest-Neighbor Prevalence'
                                switch NeighType
                                    case 'Distance'
                                        % Pull the positions of the selected cells
                                        Position = pcst_i(:, (numel(phenos)+1):end);
                                        
                                        % Find expected percentage of neighbors of each type
                                        Dat_Phe = pcst_i(:, 1:(numel(phenos)));
                                        ExpectedPhe = sum(Dat_Phe)/size(Dat_Phe, 1);

                                        NeighINX = cell(size(Dat_Phe, 1), 1);
                                        SzLim = 100;
                                        Nloop = ceil(size(Dat_Phe, 1)/SzLim);
                                        for loopi = 1:Nloop
                                            
                                            waitbar((loopi*smpl_i)/(Nloop*Ngroups), vPD, sprintf([ ...
                                                'Calculating distances; ' num2str(size(Dat_Phe, 1)*size(Dat_Phe, 1)) ' Calculations']));
                                            
                                            st =((loopi-1)*SzLim+1);
                                            en = min((st+SzLim-1), size(Dat_Phe, 1));
                                        
                                            % Find the neighbors within the radius 
                                            NeighINX(st:en) = rangesearch(Position,Position(st:en, :),NNeigh,'Distance','euclidean', 'NSMethod', 'exhaustive');
                                        end
                                        
                                        % find the percentage of neighbors of each type for each cell
                                        NeighTypes = zeros(size(Dat_Phe));
                                        for cell_i = 1:size(Dat_Phe, 1)
                                            % Number of cells within NNeigh distance of cell_i
                                            TotalCells = numel(NeighINX{cell_i});
                                            if TotalCells>1
                                                NNeighbors = sum(Dat_Phe(NeighINX{cell_i},:));
                                            else
                                                NNeighbors = Dat_Phe(NeighINX{cell_i},:);
                                            end
                                            NeighTypes(cell_i, :) = NNeighbors/TotalCells;
                                        end
                                        % loop through the cell types and find their average neighbors
                                        MeanTypes = zeros(size(Dat_Phe, 2), size(Dat_Phe, 2));
                                        for  celltype_i = 1:size(Dat_Phe, 2)
                                            if isempty(NeighTypes(Dat_Phe(:,celltype_i)==1, :))
                                                %No selected cells of this type are in any selcted regions
                                            else
                                                MeanTypes(celltype_i, :) = mean(NeighTypes(Dat_Phe(:,celltype_i)==1, :));
                                            end
                                            
                                        end
                                        
                                        ttla = ['Interaction among neighbors within ' num2str(NNeigh) ' um'];
                                        ttlb = ['Raw percentage of neighbors within ' num2str(NNeigh) ' um'];
                                                                                
                                        RatioMatrix = (MeanTypes-ExpectedPhe)./ExpectedPhe;
                                        RatioMatrix(:, ExpectedPhe==0) = NaN;
                                        RatioMatrix(ExpectedPhe==0, :) = NaN;
                                        P_Matrix = MeanTypes;
                                        
                                        % Check the percentages
                                        '----'
                                        P_Matrix
                                        sum(P_Matrix)
                                        sum(P_Matrix, 2)

                                    case 'Number'
                                        % Pull the positions of the cells
                                        Position = pcst_i(:, (numel(phenos)+1):end);
                                        
                                        % Find expected percentage of neighbors of each type
                                        Dat_Phe = pcst_i(:, 1:(numel(phenos)));
                                        ExpectedPhe = sum(Dat_Phe)/size(Dat_Phe, 1);
                                        
                                        % Find the index of the K nearest neighbors
                                        NeighINX = knnsearch(Position,Position,'K',NNeigh,'Distance','euclidean');
                                        
                                        % Define your figure titles
                                        ttla = ['Interaction among ' num2str(NNeigh) ' nearest neighbors'];
                                        ttlb = 'Raw percentage of neighbors';

                                        % find the percentage of neighbors of each type for each cell
                                        NeighTypes = zeros(size(Dat_Phe));
                                        for cell_i = 1:size(Dat_Phe, 1)
                                            NeighTypes(cell_i, :) = sum(Dat_Phe(NeighINX(cell_i,:),:))/NNeigh;                                            
                                        end
                                        
                                        % loop through the cell types and find their average neighbors
                                        MeanTypes = zeros(size(Dat_Phe, 2), size(Dat_Phe, 2));
                                        for  celltype_i = 1:size(Dat_Phe, 2)
                                            MeanTypes(celltype_i, :) = mean(NeighTypes(Dat_Phe(:,celltype_i)==1, :));
                                        end
                                        
                                        RatioMatrix = (MeanTypes-ExpectedPhe)./ExpectedPhe;
                                        RatioMatrix(:, ExpectedPhe==0) = NaN;
                                        RatioMatrix(ExpectedPhe==0, :) = NaN;
                                        P_Matrix = MeanTypes;
                                        
                                        % Check the percentages
                                        '----'
                                        P_Matrix
                                        sum(P_Matrix)
                                        sum(P_Matrix, 2)
                                        

                                end
                            case 'Pair plots of raw data'
                                region = cell(numel(ROW2_i), 1);
                                for k = 1:numel(ROW2_i)
                                    region{k} = ['R' num2str(ROW2_i(k))];
                                end
                                fig = figure;
                                fig.Color = app.GUIOPTS.btnbgclr;
                                fig.InvertHardcopy = 'off';
                                clf
                                hold on
                                Helper.pairplot(pcst_i, YLabel, region, app.net.(DataType).cmap((unique(ROW2_i, 'stable')+1), :))
                        end
%%%%%%%%%%% add the neighbor data back to AllCells
                        if true
                            switch MAPType
                                case 'Individual Heatmap for each Sample'
                                    
%                                     [size(app.data.(smplnms_i{smpl_i}).AllCells, 1); ...
%                                         size(INDons{smpl_i}, 1); ...
%                                         sum(INDons{smpl_i}); ...
%                                         size(Dat_Pre, 1); ...
%                                         size(pcst_i, 1); ...
%                                         size(Dat_Phe, 1), ...
%                                         ]
%                                     size(NeighTypes)

                                    for  celltype_i = 1:size(Dat_Phe, 2)
                                        NearestNeighborName = [Constants.other_tag '_NN_' Helper.valid_var(CellNms{celltype_i})];
                                        tmpvect = zeros(size(app.data.(smplnms_i{smpl_i}).AllCells, 1), 1);
                                        
                                        IND = find(INDons{smpl_i});
                                        IND =IND(ismember(ROW2, NRegions_i));
                                        
                                        tmpvect(IND) = NeighTypes(:, celltype_i);
                                        app.data.(smplnms_i{smpl_i}).AllCells.(NearestNeighborName) = tmpvect;
                                    end
                                    
                                case 'Combined Heatmap of all Samples'
%                                     if smpl_i==1
%                                         CombinedDat = Dat_PreMAIN(INDDatPre{smplsIND(smpl_i)}{1}(1):INDDatPre{smplsIND(smpl_i)}{1}(2),:);
%                                         CombinedROW2 = ROW2;
%                                     else
%                                         CombinedDat = [CombinedDat; Dat_PreMAIN(INDDatPre{smplsIND(smpl_i)}{1}(1):INDDatPre{smplsIND(smpl_i)}{1}(2),:)];
%                                         CombinedROW2 = [CombinedROW2; ROW2];
%                                     end
%                                     if smpl_i==numel(smplnms_i)
%                                         Dat_Pre = CombinedDat;
%                                         ROW2 = CombinedROW2;
%                                     end
                            end
                            
                        end
%%%%%%%%%%%%%%

                        % Display the log-values if selected by the user
                        if strcmp(MAPScale, 'log') && ~strcmp(CalcType, 'Pair plots of raw data')
                            RatioMatrix = real(log(RatioMatrix));
                            P_Matrix = real(log(P_Matrix));
                        end
                        if ~strcmp(CalcType, 'Pair plots of raw data')
                            %% plots on plots on plots
                            % Plot the correlation heatmap
                            fig = figure;
                            % add an export data option
                            % Create Menu bar
                            ExportMenu = uimenu(fig);
                            ExportMenu.Text = 'Export';
                            % Create an export data option
                            ExportPltDat = uimenu(ExportMenu);
                            ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                            ExportPltDat.Text = 'Export Plot Data to .csv';

                            fig.Color = app.GUIOPTS.btnbgclr;
                            fig.InvertHardcopy = 'off';
                            clf
                            % add an export data option
                            % Create Menu bar
                            ExportMenu = uimenu(fig);
                            ExportMenu.Text = 'Export';

                            % Create an export data option
                            ExportPltDat = uimenu(ExportMenu);
                            ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Corr_Heatmap');
                            ExportPltDat.Text = 'Export Plot Data to .csv';

                            % Plot the heatmap
                            hm = heatmap(RatioMatrix);
                            hm.Colormap = app.GUIOPTS.redbluecmap;
                            hm.XDisplayLabels = CellNms;
                            hm.YDisplayLabels = YLabel;

                            limits = [-1, 1]*max(abs(min(min(RatioMatrix))), abs(max(max(RatioMatrix))));
                            hm.ColorLimits = limits;

                            if strcmp(MAPType, 'Individual Heatmap for each Sample')
                                title(sprintf([ttla '\n' strrep(smplnms_i{smpl_i}, '_', ' ') '\n Regions: ' num2str(NRegions_i)]))
                            elseif strcmp(MAPType, 'Combined Heatmap of all Samples')
                                if smp_Ngroups==1
                                    title(sprintf([ttla ' \n All Samples' '\n Regions: ' num2str(NRegions_i)]))
                                else
                                    title(sprintf([ttla ' \n Group ' num2str(Sgroup_i) ' Samples' '\n Regions: ' num2str(NRegions_i)]))
                                end
                            end

                            % Plot the p-Values, if there are any
                            if contains(CalcType, 'Nearest-Neighbor Prevalence')
                                %% Plot the p-values heatmap
                                fig = figure;
                                % add an export data option
                                % Create Menu bar
                                ExportMenu = uimenu(fig);
                                ExportMenu.Text = 'Export';
                                % Create an export data option
                                ExportPltDat = uimenu(ExportMenu);
                                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                                ExportPltDat.Text = 'Export Plot Data to .csv';

                                fig.Color = app.GUIOPTS.btnbgclr;
                                fig.InvertHardcopy = 'off';
                                clf
                                % add an export data option
                                % Create Menu bar
                                ExportMenu = uimenu(fig);
                                ExportMenu.Text = 'Export';

                                % Create an export data option
                                ExportPltDat = uimenu(ExportMenu);
                                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Corr_Heatmap');
                                ExportPltDat.Text = 'Export Plot Data to .csv';

                                hm = heatmap(P_Matrix);
                                hm.Colormap = app.GUIOPTS.redbluecmap;
                                hm.XDisplayLabels = CellNms;
                                hm.YDisplayLabels = CellNms;

                                limits = [0, 1];
                                hm.ColorLimits = limits;

                                if strcmp(MAPType, 'Individual Heatmap for each Sample')
                                    title(sprintf([ttlb '\n' strrep(smplnms_i{smpl_i}, '_', ' ') '\n Regions: ' num2str(NRegions_i)]))
                                elseif strcmp(MAPType, 'Combined Heatmap of all Samples')

                                    if smp_Ngroups==1
                                        title(sprintf([ttlb ' \n All Samples' '\n Regions: ' num2str(NRegions_i)]))
                                    else
                                        title(sprintf([ttlb ' \n Group ' num2str(Sgroup_i) ' Samples' '\n Regions: ' num2str(NRegions_i)]))
                                    end
                                end
                            end

                            if group_i == 1 && Ngroups > 1
                                RatioMatrix_G1 = RatioMatrix;
                            end
                            % if the regions have groups
                            if group_i == Ngroups && Ngroups > 1
                                %% Plot the difference between the correlation heatmaps
                                fig = figure;
                                % add an export data option
                                % Create Menu bar
                                ExportMenu = uimenu(fig);
                                ExportMenu.Text = 'Export';
                                % Create an export data option
                                ExportPltDat = uimenu(ExportMenu);
                                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                                ExportPltDat.Text = 'Export Plot Data to .csv';

                                fig.Color = app.GUIOPTS.btnbgclr;
                                fig.InvertHardcopy = 'off';
                                clf
                                hm = heatmap(RatioMatrix_G1-RatioMatrix);
                                hm.Colormap = app.GUIOPTS.redbluecmap;
                                hm.XDisplayLabels = CellNms;
                                hm.YDisplayLabels = CellNms;

                                limits = [min(min(RatioMatrix_G1-RatioMatrix)), max(max(RatioMatrix_G1-RatioMatrix))];
                                [~,MAX] = max(abs(limits));
                                [~,MIN] = min(abs(limits));
                                limits(MIN) = -1*limits(MAX);
                                hm.ColorLimits = limits;

                                if strcmp(MAPType, 'Individual Heatmap for each Sample')
                                    title(sprintf(['Difference of the ' ttla 'between groups\n' strrep(smplnms_i{smpl_i}, '_', ' ')]))
                                elseif strcmp(MAPType, 'Combined Heatmap of all Samples')
                                    title(sprintf(['Difference of the ' ttla ' between groups\n All Samples']))
                                end
                            end

                            % if the samples have groups
                            if Sgroup_i == 1 && smp_Ngroups > 1
                                RatioMatrix_G1 = RatioMatrix;
                            end
                            
                            if Sgroup_i == smp_Ngroups && smp_Ngroups > 1
                                %% Plot the difference between the correlation heatmaps
                                fig = figure;
                                % add an export data option
                                % Create Menu bar
                                ExportMenu = uimenu(fig);
                                ExportMenu.Text = 'Export';
                                % Create an export data option
                                ExportPltDat = uimenu(ExportMenu);
                                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                                ExportPltDat.Text = 'Export Plot Data to .csv';

                                fig.Color = app.GUIOPTS.btnbgclr;
                                fig.InvertHardcopy = 'off';
                                clf
                                hm = heatmap(RatioMatrix_G1-RatioMatrix);
                                hm.Colormap = app.GUIOPTS.redbluecmap;
                                hm.XDisplayLabels = CellNms;
                                hm.YDisplayLabels = CellNms;

                                limits = [min(min(RatioMatrix_G1-RatioMatrix)), max(max(RatioMatrix_G1-RatioMatrix))];
                                [~,MAX] = max(abs(limits));
                                [~,MIN] = min(abs(limits));
                                limits(MIN) = -1*limits(MAX);
                                hm.ColorLimits = limits;

                                if strcmp(MAPType, 'Individual Heatmap for each Sample')
                                    title(sprintf(['Difference of the ' ttla 'between groups\n' strrep(smplnms_i{smpl_i}, '_', ' ')]))
                                elseif strcmp(MAPType, 'Combined Heatmap of all Samples')
                                    title(sprintf(['Difference of the ' ttla ' between groups\n All Samples']))
                                end
                            end
                        end % end if not pair plots
                    end % End if plot
                end % end region groups loop

            end
        end
        close(vPD)

        function cancel_waitbar_callback(hObject)
            delete(ancestor(hObject, 'figure'));
        end
    end % end of neares neighbor backend

end  