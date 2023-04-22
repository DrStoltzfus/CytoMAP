classdef Statistics
% Statistics defines functions which help visualize the processed data and
% regions
%
% Last edited by CS 2019/01/24

    methods(Static)

        function Import_Definitions_Func
% % %             import Plotting.*;
% % %             import Analysis.*;
% % %             import Helper.*;
        end

        function heatmap_func(app, web)
            % HEATMAP_FUNC Front-End of making heatmap function. Allows user to create
            % heatmap of Phenotypes/Channels within specific samples
            % (either seperate or combined throughout samples).
            %
            % Input:
            %   - app - Instance of CytoMAP
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

            if ~Helper.any_sample(app) || ~Helper.any_net(app)
                return;
            end

            [MFI_TYPE, sample] = Helper.find_MFI(app);
            if isempty(MFI_TYPE)
                errordlg("In order to run this operation at least one of your samples had to be scanned beforehand.");
                return;
            end
            tmpData = Helper.populate_table(app, 'MFI', MFI_TYPE, 'smpls', {sample}, 'fill_checkbox', false);
            % Cells
            tData = cell(size(tmpData, 1) + 1, 6);
            tData(1:size(tmpData, 1), 1)= tmpData(:, 2);
            tData(1:size(tmpData, 1), 2) = tmpData(:, 3);

% %             % Manually defined regions
% %             if any(startsWith(fieldnames(app.data.(tmpData{[tmpData{:, 7}]==1, 8}).(MFI_TYPE)), Constants.neigh_tag))
% %                 IND_nei = startsWith(tmpData(:, 6), Constants.neigh_tag);
% %                 tmpData(IND_nei, 6) = strrep(tmpData(IND_nei, 6), Constants.neigh_tag, 'Manually Defined Regions: ');
% %             end

            % Channels
            tData(1:size(tmpData, 1), 3) = tmpData(:, 5);
            tData(1:size(tmpData, 1), 4) = tmpData(:, 6);

            % Samples
            tData(1:size(tmpData, 1), 5) = tmpData(:, 7);
            tData(1:size(tmpData, 1), 6) = tmpData(:, 8);
            tData(end, :) = {false, 'Select All', false, 'Select All', false, 'Select All'}; 

            % Build the clustering options menus
            UIfig = uifigure('Name', 'Plot Heatmap Options', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end
            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 1000 800];

            % Select Neighborhood type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Clustering Model:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+150 45 400 15];
            DataType = uidropdown(UIfig);
            % if there are any manually defined neighborhoods add the option to use them        
% % %             if any(startsWith(fieldnames(app.data.(tData{[tData{:,5}]==1, 6}).(MFI_TYPE)), Constants.neigh_tag))
            dattypes = [fieldnames(app.net); {'Manually Defined Regions (RSN)'; 'Manually Defined Regions (CCN)'}];
% % %             else
% % %                 dattypes = fieldnames(app.net);
% % %             end
            DataType.Items = dattypes;
            DataType.Value = DataType.Items{1};
            DataType.Position = alpha*[10+175+150 15 175 30];
            DataType.BackgroundColor = app.GUIOPTS.bgclr;

            % Create the table of options
            t = uitable(UIfig);
            t.Data = tData;
            t.Position = alpha*[0 70 1000 730];
            t.ColumnName = {'Include in Plot','Phenotype (Must be in all samples)', ...
                            'Include in Plot', 'Channel MFI','Plot HM for', 'Sample'};
            t.ColumnEditable = [true false true false true false];
            t.ColumnWidth = {alpha*100, alpha*350, alpha*100, alpha*125, alpha*100, alpha*200};
            t.CellEditCallback = @(dd, p) switched_sample(app, DataType, dd, p);

            % Select Data Preperation
        % % %             lbl = uilabel(UIfig); lbl.Text = sprintf('Select Data Preperation:');
        % % %             Helper.func_SetCLR(app, lbl, 'label')
        % % %             lbl.Position = alpha*[10+75 45 400 15];
            DataPrep = uidropdown(UIfig);
            DataPrep.Items = Constants.neigh_cell_norm;
            DataPrep.Value = {'Composition: Number of Cells / Total cells in Neighborhood'};
            DataPrep.Position = alpha*[10+75 37 175 30];
            DataPrep.BackgroundColor = app.GUIOPTS.bgclr;

            % Select Data Preperation for MFI
            DataPrepMFI = uidropdown(UIfig);
            DataPrepMFI.Items = Constants.neigh_mfi_norm;
            DataPrepMFI.Value = {'MFI normalized to max MFI per neighborhood'};
            DataPrepMFI.Position = alpha*[10+75 7 175 30];
            DataPrepMFI.BackgroundColor = app.GUIOPTS.bgclr;
            Helper.func_SetCLR(app, DataPrepMFI, 'button')

            % Select Heatmap type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Heatmap Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+150 45 400 15];
            MAPType = uidropdown(UIfig);
            MAPType.Items = {'Individual Heatmap for each Sample', ...
                              'Combined Heatmap of all Samples'};
            MAPType.Value = {'Individual Heatmap for each Sample'};
            MAPType.Position = alpha*[10+175+175+150 15 175 30];
            MAPType.BackgroundColor = app.GUIOPTS.bgclr;

            % Select Color Scale
            lbl = uilabel(UIfig); lbl.Text = sprintf('Scale:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+150+175 45 100 15];
            MAPScale = uidropdown(UIfig);
            MAPScale.Items = {'linear', ...
                              'log'};
            MAPScale.Value = {'log'};
            MAPScale.Position = alpha*[10+175+175+150+175 15 100 30];
            MAPScale.BackgroundColor = app.GUIOPTS.bgclr;

            % Create a Norm per sample or per dataset label
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Normalize per:');
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10 50 150 20];
            % Create a Norm per sample or per dataset button
            NormPer = uibuttongroup(UIfig, 'Visible','off');
            NormPer.Position = alpha*[10, 7, 75, 45];
            Helper.func_SetCLR(app, NormPer, 'UICpopup')
            % Create two radio buttons in the button group.
            r1 = uiradiobutton(NormPer, 'Text','Sample',...
                'Position',alpha*[5, 25, 75, 20]);
            Helper.func_SetCLR(app, r1, 'table')
            r2 = uiradiobutton(NormPer, 'Text','Dataset',...
                'Position',alpha*[5, 5, 75, 20]);
            Helper.func_SetCLR(app, r2, 'table')
            NormPer.Visible = 'on';

            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap(t, app, DataType, DataPrep, MAPType, MAPScale, NormPer, DataPrepMFI));
            btn.Position = alpha*[1000-150, 10, 100, 50];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')
            
            DataType.ValueChangedFcn = @(btn, event) ChangeDataType(app, t, DataType);

            function ChangeDataType(app, t, DataType)
                % changes visiblity of manual choice numerical input.
                p.Indices = [1 5];
                p.EditData = 0;
                switched_sample(app, DataType, t, p);
            end
            
            function switched_sample(app, DataType, dd, p)
                if p.Indices(2) == 5
                    % Make sure that at least one thing is selected
                    if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 5)), 5)))
                        dd.Data{p.Indices(1), 5} = true;
                        return;
                    end
                    % Process Select All button
                    fill_select_all = false;
                    if p.Indices(1) == find(strcmp(dd.Data(:, 6), 'Select All'))
                        ind = ~cellfun('isempty', dd.Data(:, 5));
                        dd.Data(ind, 5) = {logical(p.NewData)};
                        if ~logical(p.NewData)
                            dd.Data{strcmp(dd.Data(:, 6), app.DataN.Value), 5} = true;
                        end
                        fill_select_all = logical(p.NewData);
                    end
                    if strcmp(DataType.Value, 'Manually Defined Regions (RSN)')
                        scan_type = 'MFIRSN';
                    elseif strcmp(DataType.Value, 'Manually Defined Regions (CCN)')
                        scan_type = 'MFICCN';
                    else
                        if strcmp(app.net.(DataType.Value).userdata.DataType, 'Cell Centered Neighborhoods')
                            scan_type = 'MFICCN';
                        elseif strcmp(app.net.(DataType.Value).userdata.DataType, 'Raster Scanned Neighborhoods')
                            scan_type = 'MFIRSN';
                        elseif strcmp(app.net.(DataType.Value).userdata.DataType, 'Individual Cells')
                            scan_type = 'AllCells';
                        end
                    end
                    ind = dd.Data(:, 5);
                    ind = ind(~cellfun('isempty', ind));
                    ind = cell2mat(ind);

                    if isempty(dd.Data(ind, 6))
                        dd.Data(p.Indices(1), p.Indices(2)) = {true};
                        return;
                    end

                    tmpDat = cell(size(dd.Data, 1) - 1, 8);

                    tmpDat(1:end, 1) = {1}; % Stub. Needed for populate table
                    tmpDat(1:end, 2) = dd.Data(1:end-1, 1);
                    tmpDat(1:end, 3) = dd.Data(1:end-1, 2);
                    tmpDat(1:end, 4) = {1}; % Stub. Needed for populate table
                    tmpDat(1:end, 5) = dd.Data(1:end-1, 3);
                    tmpDat(1:end, 6) = dd.Data(1:end-1, 4);
                    tmpDat(1:end, 7) = dd.Data(1:end-1, 5);
                    tmpDat(1:end, 8) = dd.Data(1:end-1, 6);

                    tmpDat = Helper.populate_table(app, ...
                        'smpls', dd.Data(ind, 6), ...
                        'MFI', scan_type, ...
                        'prev_table', tmpDat ...
                    );

                    % Only re-make the data table if you have to
                    INDddCL = ~cellfun('isempty', dd.Data(:, 2));
                    INDddCL(end) = false;
                    INDtmpCL = ~cellfun('isempty', tmpDat(:, 3));
                    INDddMFI = ~cellfun('isempty', dd.Data(:, 4));
                    INDddMFI(end) = false;
                    INDtmpMFI = ~cellfun('isempty', tmpDat(:, 6));

                    if ~Helper.setequal(dd.Data(INDddCL, 2), tmpDat(INDtmpCL, 3)) || ~Helper.setequal(dd.Data(INDddMFI, 4), tmpDat(INDtmpMFI, 6))
                        dd.Data = cell(size(tmpDat, 1) + 1, 6);
                        dd.Data(1:numel(tmpDat(:, 2)), 1) = tmpDat(:, 2);
                        dd.Data(1:numel(tmpDat(:, 3)), 2) = tmpDat(:, 3);
                        dd.Data(1:numel(tmpDat(:, 5)), 3) = tmpDat(:, 5);
                        dd.Data(1:numel(tmpDat(:, 6)), 4) = tmpDat(:, 6);
                        dd.Data(1:numel(tmpDat(:, 7)), 5) = tmpDat(:, 7);
                        dd.Data(1:numel(tmpDat(:, 8)), 6) = tmpDat(:, 8);
                        dd.Data(end, :) = {false, 'Select All', false, 'Select All', fill_select_all, 'Select All'}; 
                    end
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
                    ind = ~cellfun('isempty', dd.Data(:, 2));
                    dd.Data(ind, 1) = {logical(p.NewData)};
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 4), 'Select All')) && p.Indices(2) == 3
                    ind = ~cellfun('isempty', dd.Data(:, 3));
                    dd.Data(ind, 3) = {logical(p.NewData)};
                end
            end

            function backend_wrap(t, app, DataType, DataPrep, MAPType, MAPScale, NormPer, DataPrepMFI)
                tmp_data = t.Data(1:end-1, :);
                Statistics.heatmap_backend( ...
                    app, ...
                    tmp_data([tmp_data{:, 5}]'==1, 6), ... Samples
                    tmp_data([tmp_data{:, 1}]'==1, 2), ... Phenotypes
                    tmp_data([tmp_data{:, 3}]'==1, 4), ... MFIs
                    DataType.Value, ...
                    DataPrep.Value, ...
                    MAPType.Value, ...
                    MAPScale.Value, ...
                    NormPer.SelectedObject.Text, ...
                    DataPrepMFI.Value ...
                );
            end
        end

        function heatmap_backend(app, smplnms, phenos, MFIList, DataType, DataPrep, MAPType, MAPScale, NormOpt, DataPrepMFI, web)
            % HEATMAP_BACKEND Back-End of making heatmap function. Allows user to create
            % heatmap of Phenotypes/Channels within specific samples
            % (either seperate or combined throughout samples).
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smplnms - cell - Samples from which data to make heatmap
            %       will be pulled.
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
            %   - NormOpt - string, weather to normalize by sample or by
            %       dataset
            %
            % Check Samples-Models compatibility

            % Check to make sure something was selected
            if isempty(phenos) && isempty(MFIList)
                return
            end

            if nargin <= 10
                web = 0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

            % make sure the selected model or regions are in all sampes
            [smplnms, type] = Helper.check_model(app, smplnms, DataType, MFIList);
            
            switch type
                case 'Raster Scanned Neighborhoods'
                    type_sub = 'MFIRSN';
                case 'Cell Centered Neighborhoods'
                    type_sub = 'MFICCN';
                case 'Individual Cells'
                    type_sub = 'AllCells';
            end
            
            % Pull the manual region names
            if startsWith(DataType, 'Manually Defined Regions')
                % Pull all available neighborhood names
                Man_names = fieldnames(app.data.(smplnms{1}).(type_sub));
                indnms = startsWith(Man_names, Constants.neigh_tag);
                Man_names = Man_names(indnms);
                % Make sure those names are in all samples
                for smpl_i=2:numel(smplnms)
                    % Pull all neighborhood names
                    Man_names_tmp = fieldnames(app.data.(smplnms{smpl_i}).(type_sub));
                    indnms = startsWith(Man_names_tmp, Constants.neigh_tag);
                    Man_names_tmp = Man_names_tmp(indnms);
                    % if the sets are not equal take the more restrictive one
                    if any(~contains(Man_names, Man_names_tmp))
                        Man_names = Man_names_tmp;
                    end
                end
                % have user select manual regions
                [indx,tf] = listdlg('ListString',Man_names);
                if ~tf
                    return
                end
                ManRegs = Man_names(indx);
            end
            
            %% Load the data                
            [Dat_PreMAIN, ~, INDDatPre, ~, ~] = Helper.func_loaddata(app, ...
                smplnms, ...
                phenos, ...
                Helper.valid_channel(MFIList), ...
                [], ...
                [], ...
                type, ...
                DataPrep, ...
                NormOpt, ...
                0, ...
                DataPrepMFI);
            switch MAPType
                case 'Combined Heatmap of all Samples'
                    Dat_Pre = Dat_PreMAIN;
            end

            %% Plot the heatmap
            for i=1:numel(smplnms)

                % Make a waitbar
                vPD = waitbar(0, 'Plotting heatmap', ...
                    'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));
                vPD.Position(4) = 85;
                
                CellNms = phenos;
                switch MAPType
                    case 'Individual Heatmap for each Sample'
                        Dat_Pre = Dat_PreMAIN(INDDatPre{i}{1}(1):INDDatPre{i}{1}(2),:);
                end

                % Load the Cell data
                if startsWith(DataType, 'Manually Defined Regions')
                    for reg_i = 1:numel(ManRegs)
                        if reg_i==1
                            ROW2 = app.data.(smplnms{i}).(type_sub).(ManRegs{reg_i});
                        else
                            reglogic = app.data.(smplnms{i}).(type_sub).(ManRegs{reg_i});
                            ROW2 = ROW2 + (reg_i+max(max(ROW2))).*reglogic;
                        end
                    end
                    RowMAP = jet(max(max(ROW2)) + 1);
                else
                    ROW2 = app.data.(smplnms{i}).(type_sub).([Constants.other_tag, DataType]);
                    RowMAP = app.net.(DataType).cmap;

                end

                if strcmp(DataPrep, 'Composition+Density: Number of Cells / Total cells in Neighborhood; Density of Neighborhood')
                    CellNms = [CellNms; 'Total Cellular Density'];
                end

                % Load the MFI data
                CellNms = [CellNms', MFIList'];

                switch MAPType
                    case 'Combined Heatmap of all Samples'
                        if i==1
                            CombinedROW2 = ROW2;
                        else
                            CombinedROW2 = [CombinedROW2; ROW2];
                        end
                        if i==numel(smplnms)
                            [CombinedROW2, IND2] = sort(CombinedROW2);

                            if strcmp(MAPScale, 'log')
                                Dat_Pre = real(log(Dat_Pre));
                            end
                                Dat_Pre = [CombinedROW2, CombinedROW2, Dat_Pre(IND2, :)];
                                Dat_Pre(CombinedROW2==0, :) = []; 
                            ROW2 = CombinedROW2;
                        end
                    case 'Individual Heatmap for each Sample'
                        % This will only import the data that was used by the NN to
                        % cluster the neighborhoods, it should also sort them into the
                        % same order as was used in the NN
                        if strcmp(MAPScale, 'log')
                            Dat_Pre = real(log(Dat_Pre));
                        end

                        %Add colorbars on top of graph
                        [ROW2, IND] = sort(ROW2);
                        Dat_Pre = [ROW2, ROW2, Dat_Pre(IND, :)]; %add information of color
                        Dat_Pre(ROW2==0, :) = []; % remove neighborhoods with no cells
                end

                % Do the plot if its the last sample or you are making
                % plots for every sample
                if strcmp(MAPType, 'Individual Heatmap for each Sample') || i == numel(smplnms)
                    % Plot the data heatmap
                    fig = figure;
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
                    if web==1
                        fig.Visible = 'off';
                    end

                    if strcmp(DataType, 'Manually Defined Regions')
                        NMS = [strrep(strrep(strrep(strrep(CellNms, 'Gate_', ''), 'POS', '+'), 'NEG', '-'), 'Ch', '')];
                        NMS = strrep(NMS,  Constants.neigh_tag, 'Region: ');
                        NMS = strrep(NMS,  '_', '');
                        Dat_Pre = Dat_Pre(:, 3:end);
                        NMS = [NMS(:, (end-numel(ManRegs)+1):end), NMS(:, 1:(end-numel(ManRegs)))];
                        dat_regions = Dat_Pre(:, (end-numel(ManRegs)+1):end);
% %                         unique(dat_regions)
% %                         dat_regions = max(max(Dat_Pre(:, 1:(end-numel(ManRegs))))).*(dat_regions./(max(max(dat_regions))));
                        Dat_Pre = [dat_regions, Dat_Pre(:, 1:(end-numel(ManRegs)))];
                    else
                        NMS = [{'', 'Region:'}, strrep(strrep(strrep(strrep(strrep(CellNms, 'Gate_', ''), '_', ' '), 'POS', '+'), 'NEG', '-'), 'Ch', '')];
                    end

                    clf
                    % Create the image
                    imagesc(Dat_Pre);

                    cmap = app.GUIOPTS.redbluecmap;
                    cmap(1:round(size(cmap,1)/2), :) = [];
                    colormap(cmap)

                    c=colorbar;
                    if strcmp(MAPScale, 'log')
                        c.Label.String = ['log( ' DataPrep ' )'];
                    elseif strcmp(MAPScale, 'linear')
                        c.Label.String =  DataPrep ;
                    end

                    axis tight
                    ax = gca;
                    ax.XColor = app.GUIOPTS.txtclr;
                    ax.YColor = app.GUIOPTS.txtclr;
                    ax.XTick = 1:size(Dat_Pre, 2);
                    ax.XTickLabel = NMS;
                    ax.YTick = [];
                    ax.YTickLabel = [];
                    ax.YColor = app.GUIOPTS.bgclr;
                    ax.YDir = 'Normal';
                    ax.XDir = 'Normal';

                    if strcmp(DataPrep, 'Standardize: subtract Mean, divide by standard deviation')
                        limits = [min(min(Dat_Pre(:, 3:end))), max(max(Dat_Pre(:, 3:end)))]-1;
                        [~,MAX] = max(abs(limits));
                        [~,MIN] = min(abs(limits));
                        limits(MIN) = -1*limits(MAX); % Make the color limits symetrical
                        colormap(app.GUIOPTS.redbluecmap)
                        caxis(limits);
                    else
                        limits = [min(min(Dat_Pre(:, 3:end))), max(max(Dat_Pre(:, 3:end)))];
                        if limits(1) == -Inf
                            limits(1) = -6;
                        end
                        if limits(2) > 0
                            limits(1) = 0;
                        end
                        caxis(limits);
                    end
                    % Function that changes the Y axis order
                    heat_order_btn = uicontrol(fig, ...
                        'Style','pushbutton', ...
                        'Position', alpha*[Constants.plt_menu_x0 Constants.plt_menu_y0 130 25], ...
                        'String', 'Change Y axis Order' ...
                    );
                    heat_order_btn.Callback = @(~, ~) change_heat_order(app, ax);

                    if ~strcmp(DataType, 'Manually Defined Regions')
                        % Function that changes the region order on the plot
                        heat_order_btn_2 = uicontrol(fig, ...
                            'Style','pushbutton', ...
                            'Position', alpha*[Constants.plt_menu_x0+130 Constants.plt_menu_y0 130 25], ...
                            'String', 'Change Region Order' ...
                        );
                        heat_order_btn_2.Callback = @(~, ~) change_color_order(app, ax, RowMAP);
                        % Function that downsamples the plotted data
                        heat_order_btn_3 = uicontrol(fig, ...
                            'Style','pushbutton', ...
                            'Position', alpha*[Constants.plt_menu_x0+130+130 Constants.plt_menu_y0 130 25], ...
                            'String', 'Downsample Plot' ...
                        );
                        heat_order_btn_3.Callback = @(~, ~) change_factor(app, ax, RowMAP);
                        % Make a button that resamples region columns such that
                        % they are the same size o the plot
                        heat_order_btn_4 = uicontrol(fig, ...
                            'Style','pushbutton', ...
                            'Position', alpha*[Constants.plt_menu_x0+130+130+130 Constants.plt_menu_y0 130 25], ...
                            'String', 'Evenly Space Regions' ...
                        );
                        heat_order_btn_4.Callback = @(~, ~) even_columns(ax, RowMAP);
                    end
                    if strcmp(MAPType, 'Individual Heatmap for each Sample')
                        title(sprintf(['Sorted Neighborhoods \n' strrep(smplnms{i}, '_', ' ')]))
                    elseif strcmp(MAPType, 'Combined Heatmap of all Samples')
                        title(sprintf('Sorted Neighborhoods \n All Samples'))
                    end

                    box off
                    hold on
                    if ~strcmp(DataType, 'Manually Defined Regions')
                        ax.XLim = [1.5 numel(CellNms)+2.5];

                        %Add colorbars on top of graph
                        icolor=Dat_Pre(:,1);
                        for ik=unique(icolor(icolor~=0))'
                            for ij = 1.75:0.05:2.25
                                plot([ij ij], [find(icolor(icolor~=0)==ik,1, 'first'), find(icolor(icolor~=0)==ik,1, 'last')], ...
                                    '-', 'Color', RowMAP(ik+1,:), 'LineWidth', 10)
                            end
                        end
                    end
                    camroll(-90)
                    clear pcst

                end % End if plot


                waitbar(1,vPD, 'Done!');
                Helper.func_closeVPD(app, vPD)
            end

            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end

            function change_heat_order(app, ax, web)
                if nargin<3
                    web=0;
                end
                alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

                names = ax.XTickLabel(2:end);
                % Make sure that order list is synchronized with phenos from table.
                uifig = uifigure('Position', alpha*[100 100 300 250]);
                if web==1
                    uifig.Visible='OFF';
                end
                heat_axis_btn = uilistbox(uifig, ...
                    'Items', names, ...
                    'Position', alpha*[0 50 uifig.Position(3) uifig.Position(4) - 50] ...
                );
                uibutton(uifig, 'Position', alpha*[0 25 uifig.Position(3) / 2.0 25], 'Text', 'Move Up', 'ButtonPushedFcn', @(btn, p) move_up(heat_axis_btn, ax));
                uibutton(uifig, 'Position', alpha*[0 0 uifig.Position(3) / 2.0 25], 'Text', 'Move Top', 'ButtonPushedFcn', @(btn, p) move_top(heat_axis_btn, ax));
                uibutton(uifig, 'Position', alpha*[uifig.Position(3) / 2.0, 25, uifig.Position(3) / 2.0, 25], 'Text', 'Move Down', 'ButtonPushedFcn', @(btn, p) move_down(heat_axis_btn, ax));
                uibutton(uifig, 'Position', alpha*[uifig.Position(3) / 2.0, 0, uifig.Position(3) / 2.0, 25], 'Text', 'Move Bottom', 'ButtonPushedFcn', @(btn, p) move_bottom(heat_axis_btn, ax));

                function move_up(btn, ax)
                    idx = find(strcmp(btn.Items, btn.Value));
                    order = 1:numel(btn.Items);
                    if idx == 1 || idx == 2
                        return;
                    end

                    order(idx) = idx - 1;
                    order(idx - 1) = idx;
                    tmp = btn.Items(idx - 1);
                    btn.Items(idx - 1) = {btn.Value};
                    btn.Items(idx) = tmp;

                    order = cat(2, 1, order + 1);
                    ax.Children(end).CData = ax.Children(end).CData(:, order);
                    ax.XTickLabel = ax.XTickLabel(order);
                end

                function move_down(btn, ax)
                    idx = find(strcmp(btn.Items, btn.Value));
                    order = 1:numel(btn.Items);
                    if idx == numel(btn.Items) || idx ==1
                        return;
                    end

                    order(idx) = idx + 1;
                    order(idx + 1) = idx;
                    tmp = btn.Items(idx + 1);
                    btn.Items(idx + 1) = {btn.Value};
                    btn.Items(idx) = tmp;

                    order = cat(2, 1, order + 1);
                    ax.Children(end).CData = ax.Children(end).CData(:, order);
                    ax.XTickLabel = ax.XTickLabel(order);
                end

                function move_top(btn, ax)
                    idx = find(strcmp(btn.Items, btn.Value));
                    order = 1:numel(btn.Items);
                    if idx == 1  || idx ==2
                        return;
                    end

                    order(idx) = 2;
                    order(2) = idx;
                    tmp = btn.Items(2);
                    btn.Items(2) = {btn.Value};
                    btn.Items(idx) = tmp;

                    order = cat(2, 1, order + 1);
                    ax.Children(end).CData = ax.Children(end).CData(:, order);
                    ax.XTickLabel = ax.XTickLabel(order);
                end

                function move_bottom(btn, ax)
                    idx = find(strcmp(btn.Items, btn.Value));
                    order = 1:numel(btn.Items);
                    if idx == numel(btn.Items)  || idx ==1
                        return;
                    end

                    order(idx) = numel(btn.Items);
                    order(end) = idx;
                    tmp = btn.Items(end);
                    btn.Items(end) = {btn.Value};
                    btn.Items(idx) = tmp;

                    order = cat(2, 1, order + 1);
                    ax.Children(end).CData = ax.Children(end).CData(:, order);
                    ax.XTickLabel = ax.XTickLabel(order);
                end
            end

            function change_color_order(app, ax, RowMAP, web)
                if nargin<4
                    web=0;
                end
                alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

                names=cellstr(int2str(unique(ax.Children(end).CData(:,1))));

                % Make sure that order list is synchronized with phenos from table.
                uifig = uifigure('Position', alpha*[100 100 300 250]);
                if web==1
                    uifig.Visible='OFF';
                end

                heat_axis_btn = uilistbox(uifig, ...
                    'Items', names, ...
                    'Position', alpha*[0 50 uifig.Position(3) uifig.Position(4) - 50] ...
                );
                uibutton(uifig, 'Position', alpha*[0 25 uifig.Position(3) / 2.0 25], ...
                    'Text', 'Move Left', 'ButtonPushedFcn', @(btn, p) move_left(heat_axis_btn,ax, RowMAP));
                uibutton(uifig, 'Position', alpha*[0 0 uifig.Position(3) / 2.0 25], ...
                    'Text', 'Move Right', 'ButtonPushedFcn', @(btn, p) move_right(heat_axis_btn, ax, RowMAP));


                function move_left(btn, ax, RowMAP)
                    CDataTEMP = ax.Children(end).CData;
                    CDataTEMP = CDataTEMP';

                    % Move left one
                    RowMove = str2double(btn.Value);

                    if RowMove~=1
                        btn.Value=btn.Items{str2double(btn.Value)-1};
                        % find indeces of current group

                        INDMove = [find(CDataTEMP(1,:)==RowMove, 1 ), find(CDataTEMP(1,:)==RowMove, 1, 'last' )];
                        % Find index pf previous group
                        INDTo = find(CDataTEMP(1,:)==(RowMove-1), 1 );

                        CDataTEMP(1, INDMove(1):INDMove(2)) = RowMove-1;
                        CDataTEMP(1, INDTo:(INDMove(1)-1)) = RowMove;

                        CDataTEMP = [CDataTEMP(:, 1:INDTo-1),...
                                    CDataTEMP(:, INDMove(1):INDMove(2)), ...
                                    CDataTEMP(:, INDTo:(INDMove(1)-1)),...
                                    CDataTEMP(:, (INDMove(2)+1):end)];
                    end

                    CDataTEMP = CDataTEMP';
                    % Put the data back
                    ax.Children(end).CData = CDataTEMP;
                    % erase all of the lines
                    delete(ax.Children(1:end-1));
                    %Add colorbars on top of graph
                    color=CDataTEMP(:,2);
                    for k=unique(color(color~=0))'
                        for j = 1.75:0.05:2.25
                            plot([j j], [find(color(color~=0)==k,1, 'first'), find(color(color~=0)==k,1, 'last')], ...
                                '-', 'Color', RowMAP(k+1,:), 'LineWidth', 10);
                        end
                    end
                end

                function move_right(btn, ax, RowMAP)
                    CDataTEMP = ax.Children(end).CData;
                    CDataTEMP = CDataTEMP';
                    % Move left one
                    RowMove = str2double(btn.Value);

                    if RowMove~=max(unique(CDataTEMP(1,:)))
                        btn.Value=btn.Items{str2double(btn.Value)+1};
                        % find indeces of current group

                       INDMove = [find(CDataTEMP(1,:)==RowMove, 1 ), find(CDataTEMP(1,:)==RowMove, 1, 'last' )];
                       % Find index pf previous group
                       INDTo = find(CDataTEMP(1,:)==(RowMove+1), 1 );
                       INDToend = find(CDataTEMP(1,:)==(RowMove+1), 1, 'last' );

                       CDataTEMP(1, INDMove(1):INDMove(2)) = RowMove+1;
                       CDataTEMP(1, INDTo:INDToend) = RowMove;

                       CDataTEMP = [CDataTEMP(:, 1:INDMove(1)-1),...
                                    CDataTEMP(:, INDTo:INDToend),...
                                    CDataTEMP(:, INDMove(1):INDMove(2)), ...
                                    CDataTEMP(:, (INDToend+1):end)];

                    end

                    CDataTEMP = CDataTEMP';
                    % Put the data back
                    ax.Children(end).CData = CDataTEMP;
                    % erase all of the lines
                    delete(ax.Children(1:end-1));
                    %Add colorbars on top of graph
                    color=CDataTEMP(:,2);
                    for k=unique(color(color~=0))'
                        for j = 1.75:0.05:2.25
                            plot([j j], [find(color(color~=0)==k,1, 'first'), find(color(color~=0)==k,1, 'last')], ...
                                '-', 'Color', RowMAP(k+1,:), 'LineWidth', 10)
                        end
                    end
                end
            end

            function change_factor(app, ax, RowMAP, web)
                if nargin<4
                    web=0;
                end
                alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

                factors={'2','3','5','10','20'};

                % Make sure that order list is synchronized with phenos from table.
                uifig = uifigure('Position', alpha*[100 100 300 250]);
                if web==1
                    uifig.Visible='OFF';
                end

                heat_axis_btn = uilistbox(uifig, ...
                    'Items', factors, ...
                    'Position', alpha*[0 50 uifig.Position(3) uifig.Position(4) - 50] ...
                );
                uibutton(uifig, 'Position', alpha*[0 25 uifig.Position(3) / 2.0 25], ...
                    'Text', 'Down Sample', 'ButtonPushedFcn', @(btn, p) down_sample(heat_axis_btn,ax, RowMAP));

                function down_sample(btn, ax, RowMAP)
                    factor = str2double(btn.Value);

                    CDataTEMP = ax.Children(end).CData;
                    CDataTEMP = CDataTEMP';

                    %r = round(size(CDataTEMP,2)/factor); %Decimation factor
                    y = zeros( size(CDataTEMP,1), ceil(size(CDataTEMP,2)/factor) );

                    for y_i = 1:size(CDataTEMP,1)
                         y(y_i, :) = downsample(CDataTEMP(y_i, :),factor);
                    end
                    CDataTEMP = y;
                    CDataTEMP = CDataTEMP';
                    % Put the data back
                    ax.Children(end).CData = CDataTEMP;

                     % erase all of the lines
                    delete(ax.Children(1:end-1));
                    %Add colorbars on top of graph
                    color=CDataTEMP(:,2);
                    for k=unique(color(color~=0))'
                        for j = 1.75:0.05:2.25
                            plot([j j], [find(color(color~=0)==k,1, 'first'), find(color(color~=0)==k,1, 'last')], ...
                                '-', 'Color', RowMAP(k+1,:), 'LineWidth', 10)
                        end
                    end

                end
            end

            function even_columns(ax, RowMAP)
                %take the data out of ax    
                CDataTEMP = ax.Children(end).CData;
                color=CDataTEMP(:,2);               
                width = ones(1,size(unique(color(color~=0))',2));
                for k=unique(color(color~=0))'
                    width(k) = find(color(color~=0)==k,1,'last')-find(color(color~=0)==k,1,'first');
                end
                z=[];
                idx=1;
                for k=unique(color(color~=0))'
                    z=[z;downsample(CDataTEMP(idx:idx+width(k)-1,:),round(width(k)/min(width)))];
                    idx = idx+width(k);
                end
                CDataTEMP = z;
                % Put the data back
                ax.Children(end).CData = CDataTEMP;

                % erase all of the lines
                delete(ax.Children(1:end-1));
                color=CDataTEMP(:,2);
                %Add colorbars on top of graph
                for k=unique(color(color~=0))'
                    for j = 1.75:0.05:2.25
                        plot([j j], [find(color(color~=0)==k,1, 'first'), find(color(color~=0)==k,1, 'last')], ...
                            '-', 'Color', RowMAP(k+1,:), 'LineWidth', 10)
                    end
                end
            end
        end
    
        function interact_func(app, web)
            % INTERACT_FUNC Front-End of Interaction Plots. Allows user to see
            % dependencies between different regions of chosen model, and
            % how they are correlated with each other.
            %
            % Input:
            %   - app - Instance of CytoMAP
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            if ~Helper.any_sample(app) || ~Helper.any_net(app)
                return;
            end
            tData = cell(numel(app.DataN.Items),2);

            tData(:,1) = {true};
            tData(:,2) = app.DataN.Items;

            % Build the clustering options menus
            UIfig = uifigure('Name', 'Region Interaction Options', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end
            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 400 800];

            % Create the table of options
            t = uitable(UIfig);
            t.Data = tData;
            t.Position = alpha*[0 70 1000 730];
            t.ColumnName = {'','Sample'};
            t.ColumnEditable = [true false];
            t.ColumnWidth = {alpha*100, alpha*200};

            % Select Neighborhood type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Neighborhood Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10 40 400 15];

            DataType = uidropdown(UIfig);
            DataType.Items = fieldnames(app.net);
            DataType.Value = DataType.Items{1};
            DataType.Position = alpha*[10 5 235 30];
            DataType.BackgroundColor = app.GUIOPTS.bgclr;

            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Statistics.interact_backend( ...
                app, ...
                t.Data([t.Data{:, 1}]'==1, 2), ... Samples
                DataType.Value ...
            ));
            btn.Position = alpha*[400-100, 5, 100, 50];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')
        end

        function interact_backend(app, smplnms, DataType)
            % INTERACT_BACKEND Back-End of Interaction Plots. Allows user to see
            % dependencies between different regions of chosen model, and
            % how they are correlated with each other.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smplnms - cell - Names of samples to use in visualization
            %       of regions.
            %   - DataType - string - Name of the Model to use in creating
            %       the interaction graph. Has to exist as field under
            %       app.net.
            
            is_in = ones(numel(smplnms), 1);
            for smpl_idx=1:numel(smplnms)
                if strcmp(app.net.(DataType).userdata.DataType, 'Raster Scanned Neighborhoods')
                    if ~ismember(...
                            [Constants.other_tag, DataType], ...
                            app.data.(smplnms{smpl_idx}).MFIRSN.Properties.VariableNames...
                    )
                        is_in(smpl_idx) = 0;
                    end
                elseif strcmp(app.net.(DataType).userdata.DataType, 'Cell Centered Neighborhoods')
                    if ~ismember(...
                            [Constants.other_tag, DataType], ...
                            app.data.(smplnms{smpl_idx}).MFICCN.Properties.VariableNames...
                    )
                        is_in(smpl_idx) = 0;
                    end
                end
            end
            is_in = logical(is_in);
            if ~all(is_in)
                if any(is_in)
                    warndlg( ...
                        "Samples: " + join(smplnms(~is_in), ", ") + ...
                        " do not contain the chosen model." + newline + ...
                        "However samples: " + join(smplnms(is_in), ",") + ...
                        " contain this model, and process will continue with only these samples." ...
                    );
                    smplnms = smplnms(is_in);
                else
                    errordlg( ...
                        "All of the samples chosen do not contain the chosen model." + newline + ...
                        "The program will abort." + newline + ...
                        "In order to choose these samples, reuse the model on these datasets."...
                    );
                    return;
                end
            end

            for i=1:numel(smplnms)
                %% Load the samples data
                if strcmp(app.net.(DataType).userdata.DataType, 'Raster Scanned Neighborhoods')
                    DatSub = app.data.(smplnms{i}).MFIRSN;
                    ROWSUB = app.data.(smplnms{i}).MFIRSN.([Constants.other_tag, DataType]);
                    r = 2 * sqrt( ...
                        (DatSub.X(2) - DatSub.X(1))^2 + ...
                        (DatSub.Y(2) - DatSub.Y(1))^2 + ...
                        (DatSub.Z(2) - DatSub.Z(1))^2 ...
                    );
                    RowMAP = app.net.(DataType).cmap;
                elseif strcmp(app.net.(DataType).userdata.DataType, 'Cell Centered Neighborhoods')
                    DatSub = app.data.(smplnms{i}).MFICCN;
                    ROWSUB = app.data.(smplnms{i}).MFICCN.([Constants.other_tag, DataType]);
                    r = app.rwindowCCN;
                    RowMAP = app.net.(DataType).cmap;
                end
                %% Create the nearest neighborhood table
                vPD = waitbar(0, 'Calculculating the percentages of bordering neighborhoods', ...
                    'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));
                vPD.Position(4) = 85;
                
                clear NBorder;
                NBorder = table;
                for j = unique(ROWSUB)'
                    NBorder.(['R' num2str(j)]) = 0.*DatSub.X;
                end
                %%%%%%%%%%%%%%%%%%%%%%%% This is a bad way to do this Matrix algebra would be better here
                % Find all bordering neighborhoods
                for k=1:numel(ROWSUB)
                    Logical = (DatSub.X-DatSub.X(k)).^2 + ...
                              (DatSub.Z-DatSub.Z(k)).^2 + ...
                              (DatSub.Y-DatSub.Y(k)).^2 <= r^2;
                    Logical = Logical & ~((DatSub.X==DatSub.X(k)) & (DatSub.Y==DatSub.Y(k)));
                    %Calculate percentage of border
                    for j = unique(ROWSUB)'
                        if sum(Logical)~=0
                            NBorder.(['R' num2str(j)])(k) = sum(ROWSUB(Logical)==j)/sum(Logical);
                        else
                            NBorder.(['R' num2str(j)])(k) = 0;
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%
                %% Make an interaction table and heatmap
                if ~isvalid(vPD)
                    return;
                end
                waitbar(0.5, vPD, 'Calulating interactions');
                NTable = table;
                NTable.Region = strcat('R', num2str(unique(ROWSUB)));
                for j=unique(ROWSUB)'
                    Tsum = (sum(table2array(NBorder(ROWSUB==j,:)))./sum(ROWSUB==j))';
                    NTable.(['R' num2str(j)]) = Tsum;
                end
                % Heatmap representation plot
                fig = figure;
                fig.Color = 'w';

                % add an export data option
                % Create Menu bar
                ExportMenu = uimenu(fig);
                ExportMenu.Text = 'Export';
                % Create an export data option
                ExportPltDat = uimenu(ExportMenu);
                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Corr_Heatmap');
                ExportPltDat.Text = 'Export Plot Data to .csv';
                
%                 hm = heatmap(log2(100.*table2array(NTable(:,2:end))));
                hm = heatmap(100.*table2array(NTable(:,2:end)));
                hm.Colormap = app.GUIOPTS.redbluecmap;
                hm.XLabel = 'Region #';
                hm.YLabel = 'Region #';
                hm.Title = strrep(smplnms{i}, '_', '');
                hm.XDisplayLabels = strcat('R', num2str(unique(ROWSUB))); % don't use 'Stable' here
                hm.YDisplayLabels = strcat('R', num2str(unique(ROWSUB)));
                
                %%Store the calculated matrices
                if isfield(app.RegInt, smplnms{i})
                    if ~isfield(app.RegInt.(smplnms{i}), DataType)
                        app.RegInt.(smplnms{i}).(DataType) = struct;
                    end
                else
                    app.RegInt.(smplnms{i}) = struct;
                    app.RegInt.(smplnms{i}).(DataType) = struct;
                end
                
                app.RegInt.(smplnms{i}).(DataType).NTable = NTable;
                app.RegInt.(smplnms{i}).(DataType).RowMAP = RowMAP;
                app.RegInt.(smplnms{i}).(DataType).ROWSUB = ROWSUB;
                                
                %% Make a network plot of the region interaction
                if ~isvalid(vPD)
                    return;
                end
                % ToDo 
                % Make the percent overlap an option
                % Make upper vs. lower and option
                % Make new interaction metric that uses neares neighbor on the neighborhoods
             
                % define the percent overlap
                olp = 0.005;
                
                waitbar(0.75, vPD, 'Making network plot');
                
                % Threshold for connecting two regions
                conn = 0.005;
                

                fig = figure;
                fig.Name = smplnms{i};
                fig.Tag = DataType;
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
                
                % Add a Z-range selection
                RNG = uicontrol('Style','edit');
                RNG.String = '0.5';
                RNG.Position = [365 38 85 15];
                Helper.func_SetCLR(app, RNG, 'UICpopup')
                RNG.Callback = @testfunc;
                
                Statistics.ForceDirectedGraph(app, fig, str2double(RNG.String))

                Helper.func_closeVPD(app, vPD)
            end
            
            function testfunc(src, event)

                Statistics.ForceDirectedGraph(app, src.Parent, str2double(src.String))
            end
                
            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end
        end
        
        function ForceDirectedGraph(app, fig, conn)
                %% Make a network plot of the region interaction
                % Convert from percentage
                conn = conn/100;
                
                NTable = app.RegInt.(fig.Name).(fig.Tag).NTable;
                RowMAP = app.RegInt.(fig.Name).(fig.Tag).RowMAP;
                ROWSUB = app.RegInt.(fig.Name).(fig.Tag).ROWSUB;
                
                figure(fig.Number)
                cla
                
                nodes = 1:numel(unique(ROWSUB));
                NodeLabel = cellstr(strcat('   R', num2str(unique(ROWSUB))));
                edgs = table2array(NTable(:,2:end)); % draw lines between boardering regions
                edgs(edgs < conn) = 0;
                G = graph(edgs, 'lower', 'OmitSelfLoops');
                WDTH = (abs(G.Edges.Weight)-min(abs(G.Edges.Weight)))/max((abs(G.Edges.Weight)-min(abs(G.Edges.Weight))));
                
                plt = plot(G, ...
                    'EdgeColor',[0.9 0.9 0.9], ...% make the edges gray
                    'EdgeAlpha',0.9);
                plt.LineWidth=7.*WDTH+0.1;        % make the linewidth equal to the bordering neighborhoods
                %Label Nodes
                labelnode(plt,nodes,NodeLabel);
                layout(plt,'force');
                %Change the marker sizes
                Tsum = zeros(numel(unique(ROWSUB)),1);
                n=0;
                for j=unique(ROWSUB)'
                    n=n+1;
                    Tsum(n) = sum(ROWSUB==j);
                end
                plt.MarkerSize = 1+50*Tsum/max(Tsum);

                %Change the marker colors
                plt.NodeColor = RowMAP((unique(ROWSUB)+1), :);
                axis off
                %% Add annotation / description
                str = sprintf(['Sample: ' strrep(fig.Name, '_', ' ') ...
                             '\nNode Color = Region' ...
                             '\nNode Size = Number of elements/neighborhoods in region' ...
                             '\nLine Thickness = Percentage of shared border' ...
                             '\nLine = lower half of matrix above percent overlap threshold:']);
                dim = [0 0 .3 .3];
                annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none');
                
        end
        
        function reg_stats_func(app, web)
             % REG_STATUS Front-End of Regions Statistics function. Allows user to
            % explore different statistics on per sample basis, or on
            % cross-sample basis. Since it's an ensamble of methods,
            % multiple plots will be created, each showing a specific
            % property of regions.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
% % %             if ~Helper.any_sample(app) || ~Helper.any_net(app)
% % %                 return;
% % %             end

            %% Build the GUI for Region Statistics
            [mfi, smpl] = Helper.find_MFI(app);
            if isempty(mfi)
                errordlg("In order to run this operation at least one of your samples had to be scanned beforehand.");
                return;
            end
            tmp = Helper.populate_table(app, 'smpls', {smpl}, 'MFI', mfi);

            tData = cell(size(tmp, 1) + 1, 6);
            tData(1:size(tmp, 1), 1) = tmp(:, 2);
            tData(1:size(tmp, 1), 2) = tmp(:, 3);
            tData(1:size(tmp, 1), 3) = tmp(:, 5);
            tData(1:size(tmp, 1), 4) = tmp(:, 6);
            tData(1:size(tmp, 1), 5) = tmp(:, 7);
            tData(1:size(tmp, 1), 6) = tmp(:, 8);
            tData(end, :) = {false, 'Select All', false, 'Select All', false, 'Select All'};

            % Build the clustering options menus
            UIfig = uifigure('Name', 'Region Statistics Options', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end
            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 1000 800];

            % Create the table of options
            t = uitable(UIfig);
            t.Data = tData;
            t.Position = alpha*[0 70 1000 730];
            t.ColumnName = {'Include in Plot','Phenotype (Must be in all samples)', ...
                            'Include in Plot', 'Channel MFI','Plot', 'Sample'};
            t.ColumnEditable = [true false true false true false];
            t.ColumnWidth = {alpha*100, alpha*350, alpha*100, alpha*125, alpha*100, alpha*200};
            
            % Create a Norm per sample or per dataset label
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Generate plots:');
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = [10+175 50+5 300 20];
            
            % Create a plot types button group
            plottypes = uibuttongroup(UIfig, 'Visible','off');
            plottypes.Position = [10+175, 10, 460, 45];
            Helper.func_SetCLR(app, plottypes, 'UICpopup')
            % Create two radio buttons in the button group.
            r1 = uicheckbox(plottypes, 'Text','Cluster Composition',...
                'Position',[5, 2.5, 150, 20], 'Value', true);
            Helper.func_SetCLR(app, r1, 'table')
            r2 = uicheckbox(plottypes, 'Text','Region Prevalence',...
                'Position',[5+150+150, 22.5, 150, 20], 'Value', true);
            Helper.func_SetCLR(app, r2, 'table')
            r3 = uicheckbox(plottypes, 'Text','Fold Change',...
                'Position',[5+150, 2.5, 150, 20], 'Value', true);
            Helper.func_SetCLR(app, r3, 'table')
            r4 = uicheckbox(plottypes, 'Text','Diff. from Mean',...
                'Position',[5+150, 22.5, 150, 20], 'Value', true);
            Helper.func_SetCLR(app, r4, 'table')
            r5 = uicheckbox(plottypes, 'Text','Normalized Mean',...
                'Position',[5, 22.5, 150, 20], 'Value', true);
            Helper.func_SetCLR(app, r5, 'table')
            r6 = uicheckbox(plottypes, 'Text','Prevalence Stacked Bars',...
                'Position',[5+150+150, 2.5, 150, 20], 'Value', true);
            Helper.func_SetCLR(app, r6, 'table')
            plottypes.Visible = 'on';

            % Select Neighborhood type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Region Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10 45 175 15];
            DataType = uidropdown(UIfig);
            % Exclude individual cell models
            modelnames = fieldnames(app.net);
            indexc = [];
% %             for model_i = 1:numel(modelnames)
% %                 if strcmp(app.net.(modelnames{model_i}).userdata.DataType, 'Individual Cells')
% %                     indexc = [indexc, model_i];
% %                 end
% %             end
            modelnames(indexc) = [];
            % if there are any manually defined neighborhoods add the option to use them        
            dattypes = [modelnames; {'Manually Defined Regions (RSN)'; 'Manually Defined Regions (CCN)'}];

            
            DataType.Items = dattypes;
            DataType.Value = DataType.Items{1};
            DataType.Position = alpha*[10 10 175 30];
            DataType.BackgroundColor = app.GUIOPTS.bgclr;

            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap(t, app, DataType, plottypes));
            btn.Position = alpha*[1000-100, 5, 100, 50];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')

            % set the callback functions
            t.CellEditCallback = @(dd, p) switched_sample(app, DataType, dd, p);
            
            DataType.ValueChangedFcn = @(btn, event) ChangeDataType(app, t, DataType);

            function ChangeDataType(app, t, DataType)
                % changes visiblity of manual choice numerical input.
                p.Indices = [1 5];
                p.EditData = 0;
                switched_sample(app, DataType, t, p);
            end
            
            function backend_wrap(t, app, DataType, plottypes)
                tmp_data = t.Data(1:end-1, :);
                tmp_pheno = tmp_data(~cellfun('isempty', tmp_data(:, 2)), 2);
                tmp_pheno = tmp_pheno([tmp_data{:, 1}]'==1);
                MFIList = tmp_data(~cellfun('isempty', tmp_data(:, 4)), 4);
                MFIList = MFIList([tmp_data{:, 3}]'==1);
 
                Statistics.reg_stats_backend( ...
                    app, ...
                    tmp_data([tmp_data{:, 5}]'==1, 6), ... Sample
                    tmp_pheno, ... Phenotypes
                    DataType.Value, ...
                    plottypes, ...
                    MFIList ...
                );
            end

            function switched_sample(app, DataType, dd, p)
                if p.Indices(2) == 5
                    % Make sure that at least one thing is selected
                    if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 5)), 5)))
                        dd.Data{p.Indices(1), 5} = true;
                        return;
                    end
                    % Process Select All button
                    fill_select_all = false;
                    if p.Indices(1) == find(strcmp(dd.Data(:, 6), 'Select All'))
                        ind = ~cellfun('isempty', dd.Data(:, 5));
                        dd.Data(ind, 5) = {logical(p.NewData)};
                        if ~logical(p.NewData)
                            dd.Data{strcmp(dd.Data(:, 6), app.DataN.Value), 5} = true;
                        end
                        fill_select_all = logical(p.NewData);
                    end
                    ind = ~cellfun('isempty', dd.Data(:, 5));
                    ind(end) = false;
                    ind(ind) = logical(cell2mat(dd.Data(ind, 5)));

                    smpls = dd.Data(ind, 6);
                    if isempty(smpls)
                        [~, smpls] = Helper.find_MFI(app);
                        smpls = {smpls};
                    end
                    
                    if strcmp(DataType.Value, 'Manually Defined Regions (RSN)')
                        scan_type = 'MFIRSN';
                    elseif strcmp(DataType.Value, 'Manually Defined Regions (CCN)')
                        scan_type = 'MFICCN';
                    else
                        if strcmp(app.net.(DataType.Value).userdata.DataType, 'Cell Centered Neighborhoods')
                            scan_type = 'MFICCN';
                        elseif strcmp(app.net.(DataType.Value).userdata.DataType, 'Raster Scanned Neighborhoods')
                            scan_type = 'MFIRSN';
                        elseif strcmp(app.net.(DataType.Value).userdata.DataType, 'Individual Cells')  
                            scan_type = 'AllCells';
                        end
                    end

                    tmpDat = cell(size(dd.Data, 1) - 1, 8);

                    tmpDat(1:end, 1) = {1}; % Stub. Needed for populate table
                    tmpDat(1:end, 2) = dd.Data(1:end-1, 1);
                    tmpDat(1:end, 3) = dd.Data(1:end-1, 2);
                    tmpDat(1:end, 4) = {1}; % Stub. Needed for populate table
                    tmpDat(1:end, 5) = dd.Data(1:end-1, 3);
                    tmpDat(1:end, 6) = dd.Data(1:end-1, 4);

                    tmpDat = Helper.populate_table(app, ...
                        'smpls', smpls, ...
                        'MFI', scan_type, ...
                        'prev_table', tmpDat);

                    dd_ph = ~cellfun('isempty', dd.Data(:, 1));
                    dd_ph(end) = false;
                    dd_ph = dd.Data(dd_ph, 2);
                    dd_mfi = ~cellfun('isempty', dd.Data(:, 3));
                    dd_mfi(end) = false;
                    dd_mfi = dd.Data(dd_mfi, 4);

                    tmp_ph = ~cellfun('isempty', tmpDat(:, 2));
                    tmp_ph = tmpDat(tmp_ph, 3);
                    tmp_mfi = ~cellfun('isempty', tmpDat(:, 5));
                    tmp_mfi = tmpDat(tmp_mfi, 6);

                    if ~Helper.setequal(dd_ph, tmp_ph) || ~Helper.setequal(dd_mfi, tmp_mfi)
                        dd.Data = cell(size(tmpDat, 1) + 1, 6);
                        dd.Data(1:numel(tmpDat(:, 2)), 1) = tmpDat(:, 2);
                        dd.Data(1:numel(tmpDat(:, 3)), 2) = tmpDat(:, 3);
                        dd.Data(1:numel(tmpDat(:, 5)), 3) = tmpDat(:, 5);
                        dd.Data(1:numel(tmpDat(:, 6)), 4) = tmpDat(:, 6);
                        dd.Data(1:numel(tmpDat(:, 7)), 5) = tmpDat(:, 7);
                        dd.Data(1:numel(tmpDat(:, 8)), 6) = tmpDat(:, 8);
                        
                        dd.Data(end, :) = {false, 'Select All', false, 'Select All', fill_select_all, 'Select All'};
                    end
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
                    ind = ~cellfun('isempty', dd.Data(:, 1));
                    dd.Data(ind, 1) = {logical(p.NewData)};
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 4), 'Select All')) && p.Indices(2) == 3
                    ind = ~cellfun('isempty', dd.Data(:, 3));
                    dd.Data(ind, 3) = {logical(p.NewData)};
                end
            end
        end
        
        function reg_stats_backend(app, smplnms, Phenotypes, DataType, plottypes, MFIList)
            % REG_STATS_BACKEND Back-End of Regions Statistics function. Allows user to
            % explore different statistics on per sample basis, or on
            % cross-sample basis. Since it's an ensamble of methods,
            % multiple plots will be created, each showing a specific
            % property of regions.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - smplnms - cell - Samples from which data to make heatmap
            %       will be pulled.
            %   - Phenotypes - cell - Phenotypes on which region statistics
            %       will be calculated.
            %   - DataType - string - Name of the Model which will be used
            %       to pull regions from. Has to exist as field under
            %       app.net.
            % Pull which plot types were selected
            values = [plottypes.Children.Value];
            
            marker = {'o', 'o', 'o', 'o', 's', 's', 's', 'x', 'x', 'x', '*', '.', '+', 'p'};
            % Just in case
            marker = [marker,marker,marker,marker];
            marker = [marker,marker,marker,marker];

            % make sure the selected model or regions are in all sampes
            [smplnms, type] = Helper.check_model(app, smplnms, DataType, MFIList);
            
            switch type
                case 'Raster Scanned Neighborhoods'
                    type_sub = 'MFIRSN';
                case 'Cell Centered Neighborhoods'
                    type_sub = 'MFICCN';
                case 'Individual Cells'
                    type_sub = 'AllCells';
            end
            
            % Pull the manual region names
            if startsWith(DataType, 'Manually Defined Regions')
                % Pull all available neighborhood names
                Man_names = fieldnames(app.data.(smplnms{1}).(type_sub));
                indnms = startsWith(Man_names, Constants.neigh_tag);
                Man_names = Man_names(indnms);
                % Make sure those names are in all samples
                for smpl_i=2:numel(smplnms)
                    % Pull all neighborhood names
                    Man_names_tmp = fieldnames(app.data.(smplnms{smpl_i}).(type_sub));
                    indnms = startsWith(Man_names_tmp, Constants.neigh_tag);
                    Man_names_tmp = Man_names_tmp(indnms);
                    % if the sets are not equal take the more restrictive one
                    if any(~contains(Man_names, Man_names_tmp))
                        Man_names = Man_names_tmp;
                    end
                end
                % have user select manual regions
                [indx,tf] = listdlg('ListString',Man_names);
                if ~tf
                    return
                end
                ManRegs = Man_names(indx);
                ManDispRegs = strrep(strrep(ManRegs, Constants.neigh_tag, ''), '_', ' ');
            end
            
            
            % Pull the manual region names
            if startsWith(DataType, 'Manually Defined Regions')                
                for smpl_i=1:numel(smplnms)
                    regLogMAT = double(table2array(app.data.(smplnms{smpl_i}).(type_sub)(:, ManRegs)));
                    ROW2 = Helper.bi2de(regLogMAT);
                    if smpl_i==1
                        regs = unique(ROW2);
                        NGroups = numel(unique(ROW2));
                        RowMAP = jet(NGroups + 1);
                    else
                        regs = unique([regs; ROW2]);
                        NGroups = max(NGroups, numel(unique(ROW2)));
                        RowMAP = jet(NGroups + 1);
                    end
                end
                regs = regs(2:end);
                logic_RegCombos = Helper.de2bi(regs);
                NGroups = NGroups-1;

                ManDispRegs_nms = cell(numel(regs), 1);
                for regname_i = 1:numel(regs)
                    regnmTMP = ManDispRegs(logic_RegCombos(regname_i,:)==1);
                    if numel(regnmTMP)==1
                        ManDispRegs_nms{regname_i} = regnmTMP{1};
                    else
                        regnmTMP = strcat(regnmTMP, '+');
                        ManDispRegs_nms{regname_i} = [regnmTMP{:}];
                    end
                end
            else
                RowMAP = app.net.(DataType).cmap;
                NGroups = app.net.(DataType).NR;
                regs = [1:NGroups];
            end

            CellNms = [Phenotypes; MFIList];
            CellNms = strrep(CellNms, Constants.neigh_tag, 'Region: ');
            CellNms = strrep(CellNms, '_', '');
            % Convert to channel tags/valid ignore names
            MFIList = Helper.valid_channel(MFIList);
            MFIList(ismember(MFIList, Constants.ignore_names)) = Helper.valid_var(MFIList(ismember(MFIList, Constants.ignore_names)));
 
            % Load in the data
            rmvzeros = 0;
            NormOpt = 'Sample';
            cellprep = 'Cellularity: Number of Cells / Neighborhood';
            mfiprep = 'Sum MFI per neighborhood';
               
            [Dat_Pre, ~, INDDatPre, ~, ~] = Helper.func_loaddata( ...
                app, ...
                smplnms, ...
                Phenotypes, ...
                MFIList, ...
                [], ...
                [], ...
                type, ...
                cellprep, ...
                NormOpt, ...
                rmvzeros, ...
                mfiprep);
            
            if strcmp(type_sub, 'AllCells')
                % remove the binary cell type classifiers from the data table
%                 Dat_Pre = Dat_Pre(:, (numel(Phenotypes)+1):end);
            end
            
            % Initialize the matrices
            MeanCellsM = zeros(numel(smplnms),numel(CellNms),NGroups);
            Percentage = zeros(numel(smplnms),NGroups);
            
            % Prepare the data
            % For each region
            for j = 1:NGroups
                reg_j = regs(j);
                % For each sample
                for i=1:numel(smplnms)

                    % Pull in the sample data
                    DatSub = Dat_Pre(INDDatPre{i}{1}(1):INDDatPre{i}{1}(2),:);

                    % Load correct Row value (ROW is reduduntant and
                    % can be substituted later, for str value)
                    if startsWith(DataType, 'Manually Defined Regions')
                        regLogMAT = double(table2array(app.data.(smplnms{i}).(type_sub)(:, ManRegs)));
                        % convert from binary matrix to decimal
                        ROW = Helper.bi2de(regLogMAT);
                        
% % %                         for reg_i = 1:numel(ManRegs)
% % %                             if reg_i==1
% % %                                 ROW = app.data.(smplnms{i}).(type_sub).(ManRegs{reg_i});
% % %                             else
% % %                                 reglogic = app.data.(smplnms{i}).(type_sub).(ManRegs{reg_i});
% % %                                 ROW = ROW + (reg_i+max(max(ROW))).*reglogic;
% % %                             end
% % %                         end
                    else
                        ROW  = app.data.(smplnms{i}).(type_sub).([Constants.other_tag, DataType]);
                    end
                    % sort all of the scan data in the same order
                    if i==1
                        if isempty(DatSub(ROW==reg_j, :))
                            MeanCells = DatSub(1, :);
                        else
                            MeanCells = DatSub(ROW==reg_j, :);
                        end
                    else  % Order all data in the same way
                        if isempty(DatSub(ROW==reg_j, :)) % If there are no neighborhoods classified as region j
                            MeanCells = DatSub(1, :);
                        else
                            MeanCells = DatSub(ROW==reg_j, :);
                        end
                    end
                    if size(MeanCells, 1)~=1
                        MeanCells = mean(MeanCells);
                    end
                    
                    if isempty(DatSub(ROW==reg_j, :)) % If there are no neighborhoods classified as region j
                        MeanCellsM(i, :, j) = NaN;
                    else
                        MeanCellsM(i, :, j) = MeanCells;
                    end
                    %Calculate the percentage of sample making up that region
                    Percentage(i, j) =  numel(ROW(ROW==reg_j))/numel(ROW(ROW~=0));
                end
            end
            %% Mean makeup of each region across multiple tumors
            % Initialize figures and values
            
% %             values(5); % Plot Region prevalence wisker plot?
% %             values(6); % Plot cluster composition?
% %             values(3); % PlotDifference from mean?
% %             values(4); % Plot fold change?
% %             values(2); % Plot Normalized Max?
% %             values(1); % plot region prevalence box plot 

            % PlotDifference from mean
            if values(3)==1
                fig = figure;           
                clf
                fig.Color = app.GUIOPTS.bgclr;
                fig.InvertHardcopy = 'off';
                figN = fig.Number;
                ax = gca;
                ax.Title.String = 'Difference from mean';
                % add an export data option
                % Create Menu bar
                ExportMenu = uimenu(fig);
                ExportMenu.Text = 'Export';
                % Create an export data option
                ExportPltDat = uimenu(ExportMenu);
                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                ExportPltDat.Text = 'Export Plot Data to .csv';
                % Create a convert plot option
                convertplt = uimenu(ExportMenu);
                convertplt.MenuSelectedFcn = @(btn,event) Plt_Helper.plot_converter(app, 'Line_Point_Plot');
                convertplt.Text = 'Convert Plot to Heatmap';
            end
            % Plot fold change
            if values(4)==1
                fig = figure;
                clf
                fig.Color = app.GUIOPTS.bgclr;
                fig.InvertHardcopy = 'off';
                figN2 = fig.Number;
                ax = gca;
                ax.Title.String = 'Fold Change';
                % add an export data option
                % Create Menu bar
                ExportMenu = uimenu(fig);
                ExportMenu.Text = 'Export';
                % Create an export data option
                ExportPltDat = uimenu(ExportMenu);
                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                ExportPltDat.Text = 'Export Plot Data to .csv';
                % Create a convert plot option
                convertplt = uimenu(ExportMenu);
                convertplt.MenuSelectedFcn = @(btn,event) Plt_Helper.plot_converter(app, 'Line_Point_Plot');
                convertplt.Text = 'Convert Plot to Heatmap';
            end
            
            % Plot Normalized to Max
            if values(2)==1
                fig = figure;          
                clf
                fig.Color = app.GUIOPTS.bgclr;
                fig.InvertHardcopy = 'off';
                figN3 = fig.Number;
                ax = gca;
                ax.Title.String = 'Normalized to Max';
                % add an export data option
                % Create Menu bar
                ExportMenu = uimenu(fig);
                ExportMenu.Text = 'Export';
                % Create an export data option
                ExportPltDat = uimenu(ExportMenu);
                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                ExportPltDat.Text = 'Export Plot Data to .csv';
                % Create a convert plot option
                convertplt = uimenu(ExportMenu);
                convertplt.MenuSelectedFcn = @(btn,event) Plt_Helper.plot_converter(app, 'Line_Point_Plot');
                convertplt.Text = 'Convert Plot to Heatmap';
            end
                        
            MaxCellsM2 = max(max(MeanCellsM, [], 3, 'omitnan'), [], 1, 'omitnan'); %Find the max cellularity across the tissues, and regions
            MeanCellsM2 = mean(mean(MeanCellsM, 3, 'omitnan'), 1, 'omitnan'); %Find the average cellularity across the tissues, and regions
            MeanCellsM3 = reshape(mean(MeanCellsM, 1, 'omitnan'), numel(CellNms),NGroups); % Find the average cellularity for each region (across all tissues)
            for j = 1:NGroups
                reg_j = regs(j);
                
                % Plot cluster composition
                if values(6)==1
                    fig = figure;           
                    clf
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

                    if numel(smplnms)~=1
                        plot(mean(MeanCellsM(:, :, j), 'omitnan'), '-o' , 'Color', RowMAP(j+1,:), 'DisplayName', 'Mean of all Samples')
                        boxplot(MeanCellsM(:, :, j) , 'Color', 'k')
                    else
                        plot(MeanCellsM(:, :, j), '-' , 'Color', RowMAP(j+1,:), 'DisplayName', 'Mean of all Samples')
                    end
                    n=0;
                    for i=1:numel(smplnms)
                        n=n+1;
                        plot(MeanCellsM(i, :, j), marker{n}, 'MarkerSize', 15, 'DisplayName', strrep(smplnms{i}, '_', ' '))
                    end
                    ax = gca;
                    box off
                    ax.Color = 'w';
                    xlim([0 numel(CellNms)+1]);
                    ylabel('Mean / neighborhood')
                    ax.XTick = 1:numel(CellNms);
                    ax.XTickLabel = CellNms;
                    xtickangle(45)
                    
                    if startsWith(DataType, 'Manually Defined Regions')
                        title(['Region: ' ManDispRegs_nms{j}]);
                    else
                        if ismember('RegNames', fieldnames(app.net.(DataType)))                  
                            title(app.net.(DataType).RegNames{j});
                        else
                            title(['Region #' num2str(j)]);
                        end
                    end
                
                    legend show
                end
                
                % Plot the Deviation from the mean
                if values(3)==1
                    figure(figN)
                    hold on
                    plt = plot(MeanCellsM3(:,j)'-MeanCellsM2, '-o' , 'Color', RowMAP(j+1,:), 'LineWidth', 3);
                    
                    if startsWith(DataType, 'Manually Defined Regions')
                        plt.DisplayName = ['Region: ' ManDispRegs_nms{j}];
                    else
                        if ismember('RegNames', fieldnames(app.net.(DataType)))                  
                            plt.DisplayName = app.net.(DataType).RegNames{j};
                        else
                            plt.DisplayName = ['Region #' num2str(j)];
                        end
                    end
                    
                    xlim([1 numel(CellNms)])
                    ax = gca;
                    ax.Color = 'w';
                    ax.XTick = 1:numel(CellNms);
                    set(ax,'TickLabelInterpreter','none')
                    ax.XTickLabel = CellNms;
                    xtickangle(45)
                    ylabel('( Mean / neighborhood ) - ( All sample mean / neighborhood )');
                    grid on
                end
                
                % Plot the fold change            
                if values(4)==1
                    figure(figN2)
                    hold on
                    plt = plot(MeanCellsM3(:,j)'./MeanCellsM2, '-o' , 'Color', RowMAP(j+1,:), 'LineWidth', 3);
                    if startsWith(DataType, 'Manually Defined Regions')
                        plt.DisplayName = ['Region: ' ManDispRegs_nms{j}];
                    else
                        if ismember('RegNames', fieldnames(app.net.(DataType)))                  
                            plt.DisplayName = app.net.(DataType).RegNames{j};
                        else
                            plt.DisplayName = ['Region #' num2str(j)];
                        end
                    end
                    xlim([1 numel(CellNms)])
                    ax = gca;
                    ax.Color = 'w';
                    ax.XTick = 1:numel(CellNms);
                    set(ax,'TickLabelInterpreter','none')
                    ax.XTickLabel = CellNms;
                    xtickangle(45)
                    ylabel('$$Fold Change = { Mean / neighborhood \over  All samples mean / neighborhood }$$', ...
                            'Interpreter', 'Latex')
                    grid on
                end
                
                % Plot Normalized Max
                if values(2)==1
                    figure(figN3)
                    hold on
                    plt = plot(MeanCellsM3(:,j)'./MaxCellsM2, '-o' , 'Color', RowMAP(j+1,:), 'LineWidth', 3);
                    if startsWith(DataType, 'Manually Defined Regions')
                        plt.DisplayName = ['Region: ' ManDispRegs_nms{j}];
                    else
                        if ismember('RegNames', fieldnames(app.net.(DataType)))                  
                            plt.DisplayName = app.net.(DataType).RegNames{j};
                        else
                            plt.DisplayName = ['Region #' num2str(j)];
                        end
                    end
                    xlim([1 numel(CellNms)])
                    ax = gca;
                    ax.Color = 'w';
                    ax.XTick = 1:numel(CellNms);
                    set(ax,'TickLabelInterpreter','none')
                    ax.XTickLabel = CellNms;
                    xtickangle(45)
                    ylabel('$$Normalized Mean = { Mean / neighborhood \over  All samples max / neighborhood }$$', ...
                            'Interpreter', 'Latex')
                    grid on
                end
            end
            
            % Plot Region prevalence wisker plot
            if values(5)==1
                fig = figure;
                clf
                fig.Color = app.GUIOPTS.bgclr;
                fig.InvertHardcopy = 'off';
                if numel(smplnms)~=1
                    boxplot(100.*Percentage , 'Color', 'k')
                end
                
                ax = gca;
                ax.Title.String = 'Region Prevalence';
                hold on
                % add an export data option
                % Create Menu bar
                ExportMenu = uimenu(fig);
                ExportMenu.Text = 'Export';
                % Create an export data option
                ExportPltDat = uimenu(ExportMenu);
                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                ExportPltDat.Text = 'Export Plot Data to .csv';

                n=0;
                for i=1:numel(smplnms)
                    n=n+1;
                    plot(100.*Percentage(i, :), marker{n}, 'DisplayName', strrep(smplnms{i}, '_', ' '), 'MarkerSize', 15)
                end
                legend show
                ax = gca;
                box off
                ax.Color = 'w';
                ylabel('Percentage of sample neighborhoods belonging to region');
                ylim([0 100]);
                xlabel('Region number');
                ax.XTick = 1:1:(size(Percentage, 2));
                
                if startsWith(DataType, 'Manually Defined Regions')
                    ax.XTickLabels = ManDispRegs_nms;
                    ax.XTickLabelRotation = 45;
                else
                    if ismember('RegNames', fieldnames(app.net.(DataType)))                  
                        ax.XTickLabels = app.net.(DataType).RegNames;
                        ax.XTickLabelRotation = 45;
                    else
                        
                    end
                end
                
            end
            
            % plot region prevalence stacked bar plot
            if values(1)==1
                fig = figure;
                clf
                fig.Color = app.GUIOPTS.bgclr;
                fig.InvertHardcopy = 'off';
                ax = gca;
                ax.Title.String = 'Region Prevalence Stacked Bar Plot';
                hold on
                % add an export data option
                % Create Menu bar
                ExportMenu = uimenu(fig);
                ExportMenu.Text = 'Export';
                % Create an export data option
                ExportPltDat = uimenu(ExportMenu);
                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                ExportPltDat.Text = 'Export Plot Data to .csv';

                % Deal with if only one sample is selected
                if size(Percentage,1)==1
                    plt = bar([Percentage; zeros(1, size(Percentage,2))],1,...
                        'stacked', 'EdgeColor', 'white', 'FaceColor','flat');  
                else
                    plt = bar(Percentage,1,...
                        'stacked', 'EdgeColor', 'white', 'FaceColor','flat');
                end
                % Put in the channel colors
                for n=1:size(Percentage, 2)
                    plt(n).CData = RowMAP(n+1,:);
                    
                    if startsWith(DataType, 'Manually Defined Regions')
                        plt(n).DisplayName = ['Region: ' ManDispRegs_nms{n}];
                    else
                        if ismember('RegNames', fieldnames(app.net.(DataType)))                  
                            plt(n).DisplayName = app.net.(DataType).RegNames{n};
                        else
                            plt(n).DisplayName = ['R' num2str(n)];
                        end
                    end
                end

                box off
                ax.XLabel.String = 'Sample';
                ax.XLim = [0 numel(smplnms)+1];
                ax.YLim = [0 1];
                ax.XTick = 0:1:(numel(smplnms)+1);
                ax.XTickLabels = strrep([{''}; smplnms(:); {''}], '_', '');
                ax.XTickLabelRotation = 45;
                ax.YLabel.String = 'Region Prevalence';
                %%%%%
                legend off
                box off
                ax.Color = 'w';
                %%%%%%%%%
% % %                 clf
% % %                 hold off
% % %                 Percentage
% % %                 hm = heatmap(Percentage');
% % %                 hm.Colormap = app.GUIOPTS.redbluecmap;
% % %                 hm.XLabel = 'Region';
            end
            legend show
            
            % Make Fold change heatmap
            if values(4)==1
                
                fig = figure;
                % add an export data option
                % Create Menu bar
                ExportMenu = uimenu(fig);
                ExportMenu.Text = 'Export';
                % Create an export data option
                ExportPltDat = uimenu(ExportMenu);
                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                ExportPltDat.Text = 'Export Plot Data to .csv';

                clf
                fig.Color = app.GUIOPTS.bgclr;
                fig.InvertHardcopy = 'off';
                
                datMAT = zeros(numel(CellNms), NGroups);
                for reg_j = 1:NGroups
                    datMAT(:, reg_j) = (MeanCellsM3(:,reg_j)'./MeanCellsM2);
                end
                
                if numel(mean(Percentage))==1
                    datMAT = [datMAT; 100.*Percentage];
                else
                    datMAT = [datMAT; 100.*mean(Percentage)];
                end
                
                hm = heatmap(datMAT);
                hm.Colormap = app.GUIOPTS.redbluecmap;
                hm.XLabel = 'Region';
                
                limits = [min(min(datMAT(1:(end-1),:))), max(max(datMAT(1:(end-1),:)))]-1;
                [~,MAX] = max(abs(limits));
                [~,MIN] = min(abs(limits));
                limits(MIN) = -1*limits(MAX); % Make the color limits symetrical
                hm.ColorLimits = limits+1;
                
                % put a title on the colorbar
                axs = struct(gca); %ignore warning that this should be avoided
                cb = axs.Colorbar;
                cb.Label.String = 'Fold Change';
                
                ax = gca;
                hm.YDisplayLabels = [CellNms; {'Percentage'}];
                    
                if startsWith(DataType, 'Manually Defined Regions')
                    hm.XDisplayLabels = ManDispRegs_nms;
                else
                    if ismember('RegNames', fieldnames(app.net.(DataType)))                  
                        hm.XDisplayLabels = app.net.(DataType).RegNames;
                    else
                        
                    end
                end

            end

        end
        
        function ReduceDims_func(app, web)
            % REDUCEDIMS_FUNC Front-End of Reduce Dimensions function. It allows user to
            % run either t-SNE, U-MAP or PCA on the data, to get it to the
            % point which is able to visualize on the screen. It's useful
            % for understanding clustering results, and understanding
            % similarity between regions.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            if ~Helper.any_sample(app)
                return;
            end

            %% Build options for sorting
            [MFI_TYPE, sample] = Helper.find_MFI(app);
            if isempty(MFI_TYPE)
                MFI_TYPE = 'AllCells';
                dtype = 3;
            else
                dtype = 1;
            end
            if ~iscell(sample)
                sample = {sample};
            end
            tdata = Helper.populate_table(app, 'smpls', sample, 'MFI', MFI_TYPE, 'fill_checkbox', true);
            % For t-SNE we specifically don't want to include X,Y,Z
            tdata((end-2):end, 5) = {false};
            tdata(end + 1, :) = {[], false, 'Select All', [], false, 'Select All', false, 'Select All'}; 

            % Build the clustering options menus
            UIfig = uifigure('Name', 'Dimensionality Reduction Algorithms', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end
            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 1000 800];

            % Select Data type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Data Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+10 45 400 15];
            DataType = uidropdown(UIfig);
            DataType.Items = {'Raster Scanned Neighborhoods', 'Cell Centered Neighborhoods',  'Individual Cells'};
            DataType.Value = DataType.Items{dtype};
            DataType.Position = alpha*[10+175+10 10 175 30];
            Helper.func_SetCLR(app, DataType, 'button')

            % Create the table of options
            t = uitable(UIfig);
            t.Data = tdata;
            t.Position = alpha*[0 75 1000 725];
            t.ColumnName = {'Weight', 'Use Dimension','Phenotype (Must be in all samples)', ...
                            'Weight', 'Use Dimension', 'Channel MFI','Use for t-SNE', 'Sample'};
            t.ColumnEditable = [true true false true true false true false];
            t.ColumnWidth = {alpha*50, alpha*100, alpha*350, alpha*50, alpha*100, alpha*125, alpha*100, alpha*125};
            t.CellEditCallback = @(dd, p) switched_sample(app, DataType, dd, p);

            % Find the currently available functions
            % pull the path for the current function
            cpath = fileparts(mfilename('fullpath'));
            cpath = fullfile([cpath filesep 'Dimensionality Reduction']);
            % Get the path of the dimensionality reduction folder
            MyFolderInfo = dir(cpath);
            % remove any directory elements without dates
            MyFolderInfo = MyFolderInfo(~cellfun('isempty', {MyFolderInfo.date}));
            % remove any non-folder files
            MyFolderInfo = MyFolderInfo([MyFolderInfo.isdir]);
            % remove weird stuff
            MyFolderInfo = MyFolderInfo(~strcmp('.', {MyFolderInfo.name}));
            MyFolderInfo = MyFolderInfo(~strcmp('..', {MyFolderInfo.name}));
            
            % Select deminsoinality reduction algorithm
            lbltmp = uilabel(UIfig); lbltmp.Text = sprintf('Select Algorithm:');
            Helper.func_SetCLR(app, lbltmp, 'label')
            lbltmp.Position = alpha*[10 45 400 15];
            ClustAlgorithm = uidropdown(UIfig);
            ClustAlgorithm.Items = {MyFolderInfo.name};
            if numel(ClustAlgorithm.Items)==1
                ClustAlgorithm.Value = ClustAlgorithm.Items{1};
            else
                ClustAlgorithm.Value = ClustAlgorithm.Items{2};
            end
            ClustAlgorithm.Position = alpha*[10 10 175 30];
            Helper.func_SetCLR(app, ClustAlgorithm, 'button')

            % Select Data Preperation
            DataPrep = uidropdown(UIfig);
            DataPrep.Items = Constants.neigh_cell_norm;
            DataPrep.Value = DataPrep.Items(6);
            DataPrep.Position = alpha*[10+175+175+10+75 40 175 30];
            Helper.func_SetCLR(app, DataPrep, 'button')
                       
            % Select Data Preperation for MFI
            DataPrepMFI = uidropdown(UIfig);
            DataPrepMFI.Items = Constants.neigh_mfi_norm;
            DataPrepMFI.Value = {'Standardize: subtract Mean, divide by standard deviation'};
            DataPrepMFI.Position = alpha*[10+175+175+10+75 10 175 30];
            DataPrepMFI.BackgroundColor = app.GUIOPTS.bgclr;
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
            Helper.func_SetCLR(app, r1, 'label')
            r2 = uiradiobutton(NormPer, 'Text','Dataset',...
                'Position',alpha*[5, 5, 75, 20]);
            Helper.func_SetCLR(app, r2, 'label')
            NormPer.Visible = 'on';
            
            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap(t, app, DataPrep, DataType, ClustAlgorithm, NormPer, DataPrepMFI));
            btn.Position = alpha*[1000-100, 5, 100, 50];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')

            DataType.ValueChangedFcn = @(~, ~) ChangeDataType(app, t, DataType, DataPrep, DataPrepMFI, lbltmp);

            function ChangeDataType(app, t, DataType, DataPrep, DataPrepMFI, lbltmp)
                % changes visiblity of manual choice numerical input.
                p.Indices = [7 7];
                p.EditData = 0;
                switched_sample(app, DataType, t, p);
                if strcmp(DataType.Value, 'Individual Cells')
                    DataPrep.Visible = 'off';
                    lbltmp.Visible = 'off';
                    DataPrepMFI.Items = Constants.cell_mfi_norm;
                    DataPrepMFI.Value = DataPrepMFI.Items(1);
                else
                    DataPrep.Visible = 'on';
                    lbltmp.Visible = 'on';
                    DataPrepMFI.Items = Constants.neigh_mfi_norm;
                    DataPrepMFI.Value = DataPrepMFI.Items(3);
                end
            end

            function switched_sample(app, DataType, dd, p)
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
                    tdataTMP(end + 1, :) = {[], false, 'Select All', [], false, 'Select All', fill_select_all, 'Select All'}; 
                    dd.Data = tdataTMP;

                elseif p.Indices(1) == find(strcmp(dd.Data(:, 3), 'Select All')) && p.Indices(2) == 2
                    ind = ~cellfun('isempty', dd.Data(:, 3));
                    dd.Data(ind, 2) = {logical(p.NewData)};
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 6), 'Select All')) && p.Indices(2) == 5
                    ind = ~cellfun('isempty', dd.Data(:, 5));
                    dd.Data(ind, 5) = {logical(p.NewData)};
                end
            end

            function backend_wrap(t, app, DataPrep, DataType, ClustAlgorithm, NormPer, DataPrepMFI)
                tmp_data = t.Data(1:end - 1, :);
                INDSmpls = [tmp_data{1:numel(app.DataN.Items), 7}]==1;
                SortNames = tmp_data(1:numel([tmp_data{:,2}]),3);
                SortNames = SortNames([tmp_data{:,2}]==1,:);
                MFIList = tmp_data(1:numel([tmp_data{:,5}]),6);
                MFIList = MFIList([tmp_data{:,5}]==1,:);

                Statistics.ReduceDims_backend( ...
                    app, ...
                    tmp_data(INDSmpls, 8), ... Samples
                    SortNames, ... Phenotypes
                    [tmp_data{[tmp_data{:,2}]==1,1}], ... Phenotype weights
                    MFIList, ... MFIs
                    [tmp_data{[tmp_data{:,5}]==1,4}], ... MFI weights
                    DataPrep.Value, ...
                    DataType.Value, ...
                    ClustAlgorithm.Value, ...
                    NormPer.SelectedObject.Text, ...
                    DataPrepMFI.Value ...
                );
            end
        end

        function ReduceDims_backend(app, smplnms, SortNames, pheno_weights, MFIList, MFI_weights, DataPrep, DataType, ClustAlgorithm, NormOpt, DataPrepMFI)
            % REDUCEDIMS_BACKEND Back-End of Reduce Dimensions function. It allows user to
            % run either t-SNE, U-MAP or PCA on the data, to get it to the
            % point which is able to visualize on the screen. It's useful
            % for understanding clustering results, and understanding
            % similarity between regions.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - smplnms - cell - Samples on which data to perform
            %       dimensionality reduction.
            %   - SortNames - cell - Phenotypes on which dimensionality
            %       reduction will be performed (alongside with MFIs).
            %   - pheno_weights - array - Same size as SortNames. Array
            %       specifing how much importance to give every phenotype
            %       in the dimensionality reduction.
            %   - MFIList - cell - MFIs on which dimensionality reduction
            %       will be performed (alongside with Phenotypes).
            %   - MFI_weights - array - Same size as MFIList. Array
            %       specifing how much importance to give every MFI in the
            %       dimensionality reduction.
            %   - DataType - string - Whether to do dimensionality
            %       reduction on `Individual Cells`, or a specific model
            %       defined on Neighborhoods (has to exist as field under
            %       app.net).
            %   - DataPrep - string - How to normalize and prepare the
            %       data. Same as in Helper.Func_DataPrep.
            %   - ClustAlgorithm - string - What algorithm should be used
            %       to perform dimensionality reduction. Currently
            %       supported options are `t-SNE` and
            %       `Principal Components`.
            
            %% Prep the data for sorting
            vPD = waitbar(0, 'Preparing neighborhood data for t-SNE', ...
                'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));
            vPD.Position(4) = 85;

            switch DataType
                case 'Individual Cells'
                    type = 'Cells';
                case 'Raster Scanned Neighborhoods'
                    type = 'RSN';
                case 'Cell Centered Neighborhoods'
                    type = 'CCN';
            end

            % Define the dimesionality reduction options
%%%%%%%%%            ....
            % Pull the function name from the folder
            cpath = fileparts(mfilename('fullpath'));
            cpath = fullfile([cpath filesep 'Dimensionality Reduction']);
            funcnm = dir(fullfile([cpath filesep ClustAlgorithm], '*.m'));
            funcnm = funcnm.name(1:end-2);
            addpath(fullfile([cpath filesep ClustAlgorithm]));
            
            % name the new channel with the function used
            chname = [strrep(funcnm, 'func_', '') '_' type];
            [cls, options] = eval([funcnm '([], chname)']);

            chname = [Constants.other_tag Helper.valid_var(options{end, 2})];

            if cls==1
                Helper.func_closeVPD(app, vPD)
                return
            end
%%%%%%%%%%
            %% Load the data
            if ~isempty(MFIList)
                MFIList = Helper.valid_channel(MFIList);
            end
            % Not Individual Cells
            if ~strcmp(DataType, 'Individual Cells')
                [Dat_Pre, ~, INDDatPre, INDzrs, INDons] = Helper.func_loaddata(app, ...
                    smplnms, ...
                    SortNames, ...
                    MFIList, ...
                    pheno_weights, ...
                    MFI_weights, ...
                    DataType, ...
                    DataPrep, ...
                    NormOpt, ...
                    1, ...
                    DataPrepMFI); %<Change back to 0 maybe
            elseif strcmp(DataType, 'Individual Cells')
                [Dat_Pre, ~, INDDatPre, INDzrs, INDons] = Helper.func_loaddata(app, ...
                    smplnms, ...
                    SortNames, ...
                    MFIList, ...
                    pheno_weights, ...
                    MFI_weights, ...
                    DataType, ...
                    'Cellularity: Number of Cells / Neighborhood', ...
                    NormOpt, ...
                    1, ...
                    DataPrepMFI); 
                Dat_Pre = Dat_Pre(:, (numel(SortNames)+1):end);
            end
                        
            %% do dimensionality reduction on the neighborhoods
            if ~isvalid(vPD)
                return;
            end

            waitbar(0.5, vPD, 'Reducing Dimensions');
                       
            % run the selected dimensionality reduction function wrapper
            [tdat, ~] = eval([funcnm '(Dat_Pre, options)']);

            if isempty(tdat)
                errordlg( ...
                    "Something is wrong. nothing was returned" + newline + ...
                        "Please try different settings for this dataset, or report an issue to [WEBSITE/EMAIL LINK].", ...
                    "Empty dataset" ...
                );
                return;
            end

            %%%%%%%% Just in case it fails at the last bit
            app.CellInfo = tdat;
            %%%%%%%%
              
            %% pull the reduced dimension for each  data
            tSubDat = cell(numel(smplnms));
            for i=1:numel(smplnms)
                % Pull out the classification data for each sample
                tSubDat{i} = tdat(INDDatPre{i}{1}(1):INDDatPre{i}{1}(2), :);

                if ~strcmp(type, 'Cells')
                    % If there are were neighborhoods with 0 cells excluded
                    NNpoints = size(app.data.(smplnms{i}).(['MFI' type]), 1);
                    if NNpoints==size(tSubDat{i}, 1)
                        Class=tSubDat{i};
                    % If ROW{i} and NNeighbor do not have the same
                    % length padd with zeros where there were zero elements
                    else
                        Class = zeros(NNpoints, 2);
                        % It might be more appropriate here to use Nans but
                        % NaNs break everything
                        Class(INDzrs{i}, :) = 0.*Class(INDzrs{i}, :);
%                         Class(INDzrs{i}, :) = NaN;
                        Class(INDons{i}, :) = tSubDat{i};
                    end
                    % Redefine ROW so it is the proper length size, etc.
                    tSubDat{i} = Class;
                else
                    % If there are were neighborhoods with 0 cells excluded
                    NNpoints = size(app.data.(smplnms{i}).AllCells, 1);
                    if NNpoints==size(tSubDat{i}, 1)
                        Class=tSubDat{i};
                    % If ROW{i} and NNeighbor do not have the same
                    % length padd with zeros where there were zero elements
                    else
                        Class = zeros(NNpoints, 2);
                        % It might be more appropriate here to use Nans but
                        % NaNs break everything
                        Class(INDzrs{i}, :) = 0.*Class(INDzrs{i}, :);
%                         Class(INDzrs{i}, :) = NaN;
                        Class(INDons{i}, :) = tSubDat{i};
                    end
                    % Redefine ROW so it is the proper length size, etc.
                    tSubDat{i} = Class;
                end
            end

            %% Add data back into the dat structure
            for i=1:numel(smplnms)
                if strcmp(type, 'Cells')
                    app.data.(smplnms{i}).AllCells.([chname '_1']) = tSubDat{i}(:,1);
                    app.data.(smplnms{i}).AllCells.([chname '_2']) = tSubDat{i}(:,2);
                else
                    app.data.(smplnms{i}).(['MFI' type]).([chname '_1']) = tSubDat{i}(:,1);
                    app.data.(smplnms{i}).(['MFI' type]).([chname '_2']) = tSubDat{i}(:,2);
                end
            end

            %% Plot
            if ~isvalid(vPD)
                return;
            end
            waitbar(0.75, vPD, 'Plotting');

            Xax = [strrep(Helper.full_var(chname), Constants.other_tag, '') ' 1'];
            Yax = [strrep(Helper.full_var(chname), Constants.other_tag, '') ' 2'];

            if strcmp(type, 'Cells')
                Plotting.func_newfig(app, ...
                'X', Xax, ...
                'Y', Yax, ...
                'P', SortNames, ...
                'S', smplnms)
            else
                Plotting.func_newfig(app, ...
                'X', Xax, ...
                'Y', Yax, ...
                'P', {['Density/MFI ' type]}, ...
                'S', smplnms)
            end

            Helper.func_closeVPD(app, vPD)

            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end
        end

        function cellularity_func(app, web)
            % CELLULARITY Front-End of Plotting Cellularity. Allows user to visualize
            % sizes of the samples in a form of heatmap. It can take into
            % account volume, in order to show density etc.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            if nargin<2
                web=0;
            end
            
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            if ~Helper.any_sample(app) 
                return;
            end

            %% Make some options
            smpls = app.DataN.Items;
            
            new_dat = Helper.populate_table(app, 'smpls', smpls, 'fill_checkbox', false);
            tmp = cell(size(new_dat, 1) + 1, 8);
            tmp(1:size(new_dat, 1), 1)       = new_dat(:, 2);
            tmp(1:size(new_dat, 1), 2)       = new_dat(:, 3);
            tmp(1:size(new_dat, 1), 3)       = {false};
            % Put in sample group numbers
            tmp(1:numel(smpls), 4) = {1};
            tmp(1:numel(smpls), 5) = {100};
            tmp(1:size(new_dat, 1), 6)       = new_dat(:, 4);
            tmp(1:size(new_dat, 1), 7)       = new_dat(:, 5);
            tmp(1:numel(smpls), 8) = {0};
            tmp(end, 1:3) = {false, 'Select All', false};
            tmp(end, 6:7) = {false, 'Select All'};
            
            % Pull imaged volume and sample volume from metadata if it is
            % there
            for smp_i = 1:numel(smpls)
                if contains('MetaData', fieldnames(app.data.(smpls{smp_i})))
                    % Imaged Volume
                    if contains('ImagedVolume', fieldnames(app.data.(smpls{smp_i}).MetaData))
                        tmp{smp_i,8} = app.data.(smpls{smp_i}).MetaData.ImagedVolume;
                    end
                    % Sample Volume
                    if contains('SampleVolume', fieldnames(app.data.(smpls{smp_i}).MetaData))
                        tmp{smp_i,5} = app.data.(smpls{smp_i}).MetaData.SampleVolume;
                    end
                    if contains('Group', fieldnames(app.data.(smpls{smp_i}).MetaData))
                        grp = app.data.(smpls{smp_i}).MetaData.Group;
                        if iscell(grp)
                            tmp(smp_i,4) = grp;
                        else
                            tmp(smp_i,4) = {grp};
                        end
                    else
                        tmp(smp_i,4) = {1}; 
                    end
                else
                    tmp(smp_i,4) = {1};
                end
            end
            
            % Build the clustering options menus
            UIfig = uifigure('Name', 'Cellularity Options', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end
            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 950 800];

            % Create the table of options
            t = uitable(UIfig);
            t.Data = tmp;
            t.Position = alpha*[0 70 950 730];
            t.ColumnName = {'Plot Cellularity', 'Phenotype', 'Plot Comparisson','Group #', 'Volume', 'Include Samples', 'Sample', 'Imaged Volume'};
            t.ColumnEditable = [true false true true true true false, false];
            t.ColumnWidth = {alpha*125, alpha*350, alpha*125, alpha*50, alpha*50, alpha*125, alpha*300, alpha*50};
            t.CellEditCallback = @(dd, p) switch_phenotypes(app, dd, p);

            lbl = uilabel(UIfig);
            lbl.Text = sprintf(['Data type:' ...
                '                                                  Plot Type: ' ...
                '                                                                         Divide By:']);
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = [10 50+5 800 20];
            
            % Create a plot types button group
            plottypes = uibuttongroup(UIfig, 'Visible','off');
            plottypes.Position = [10, 10, 800, 45];
            Helper.func_SetCLR(app, plottypes, 'UICpopup')
            % Create two radio buttons in the button group.
            r1 = uicheckbox(plottypes, 'Text','Number of Cells',...
                'Position',[5, 22.5, 150, 20], 'Value', true);
            Helper.func_SetCLR(app, r1, 'table')
            r2 = uicheckbox(plottypes, 'Text','Percentage of Cells',...
                'Position',[5, 2.5, 150, 20], 'Value', false);
            Helper.func_SetCLR(app, r2, 'table')
            
            r3 = uicheckbox(plottypes, 'Text','Box Plot',...
                'Position',[5+200, 2.5, 100, 20], 'Value', false);
            Helper.func_SetCLR(app, r3, 'table')
            r4 = uicheckbox(plottypes, 'Text','Heatmap',...
                'Position',[5+200, 22.5, 100, 20], 'Value', true);
            Helper.func_SetCLR(app, r4, 'table')
            r5 = uicheckbox(plottypes, 'Text','Stacked Bar Graph',...
                'Position',[5+200+75, 22.5, 150, 20], 'Value', false);
            Helper.func_SetCLR(app, r5, 'table')
            r6 = uicheckbox(plottypes, 'Text','Box Plot By Group #',...
                'Position',[5+200+75, 2.5, 150, 20], 'Value', false);
            Helper.func_SetCLR(app, r6, 'table')
            
            r7 = uicheckbox(plottypes, 'Text','/ Median',...
            'Position',[5+200+75+200, 22.5, 150, 20], 'Value', true);
            Helper.func_SetCLR(app, r7, 'table')
            
            r8 = uicheckbox(plottypes, 'Text','/ Mean',...
            'Position',[5+200+75+200, 2.5, 150, 20], 'Value', false);
            Helper.func_SetCLR(app, r8, 'table')
            
            r9 = uicheckbox(plottypes, 'Text','/ Image Volume',...
            'Position',[5+200+75+200+75, 2.5, 150, 20], 'Value', false);
            Helper.func_SetCLR(app, r9, 'table')
            
            plottypes.Visible = 'on';
            
            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap(t, app, plottypes));
            btn.Position = alpha*[950-100, 5, 100, 50];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')

            function switch_phenotypes(app, dd, p)
                if p.Indices(2) == 6
                    fill_select_all = false;
                    if p.Indices(1) == find(strcmp(dd.Data(:, 7), 'Select All'))
                        ind = ~cellfun('isempty', dd.Data(:, 6));
                        dd.Data(ind, 6) = {logical(p.NewData)};
                    end
                    ind = ~cellfun('isempty', dd.Data(:, 6));
                    ind(end) = false;
                    ind(ind) = logical(cell2mat(dd.Data(ind, 6)));

                    new_data = cell(size(dd.Data, 1) - 1, 5);

                    ind_first_col = ~cellfun('isempty', dd.Data(:, 1));
                    ind_first_col(end) = false;
                    first_col = dd.Data(ind_first_col, 1);
                    first_col = double(cell2mat(first_col));
                    first_col = num2cell(first_col);

                    new_data(1:numel(first_col), 1) = first_col;
                    new_data(:, 2) = dd.Data(1:end - 1, 1);
                    new_data(:, 3) = dd.Data(1:end - 1, 2);
                    new_data(:, 4) = dd.Data(1:end - 1, 6);
                    new_data(:, 5) = dd.Data(1:end - 1, 7);

                    % Save
                    idx = ~cellfun('isempty', dd.Data(:, 7));
                    group_vol_keep = dd.Data(idx, 4:5);
                    img_keep = dd.Data(idx, 8);

                    smpls_tmp = dd.Data(ind, 7);
                    if isempty(smpls_tmp)
                        smpls_tmp = {app.DataN.Value};
                    end
                    new_data = Helper.populate_table(...
                        app, ...
                        'smpls', smpls_tmp, ...
                        'prev_table', new_data ...
                    );

                    dd_ph = ~cellfun('isempty', dd.Data(:, 2));
                    dd_ph(end) = false;
                    dd_ph = dd.Data(dd_ph, 2);

                    tmp_ph = ~cellfun('isempty', new_data(:, 2));
                    tmp_ph = new_data(tmp_ph, 3);

                    if ~Helper.setequal(dd_ph, tmp_ph)
                        dd.Data = cell(size(new_data, 1) + 1, 8);

                        first_col = new_data(~cellfun('isempty', new_data(:, 1)), 1);
                        first_col = logical(cell2mat(first_col));
                        first_col = num2cell(first_col);
                        dd.Data(1:numel(first_col), 1) = first_col;

                        dd.Data(1:size(new_data, 1), 2) = new_data(:, 3);
                        dd.Data(1:size(new_data, 1), 3) = new_data(:, 2);
                        dd.Data(end, 1:2) = {false, 'Select All'};

                        dd.Data(1:size(group_vol_keep, 1), 4:5) = group_vol_keep;

                        dd.Data(1:size(new_data, 1), 6) = new_data(:, 4);
                        dd.Data(1:size(new_data, 1), 7) = new_data(:, 5);
                        dd.Data(end, 6:7) = {fill_select_all, 'Select All'};
                        dd.Data(1:size(img_keep, 1), 8) = img_keep;
                    end
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
                    ind = ~cellfun('isempty', dd.Data(:, 2));
                    dd.Data(ind, 1) = {logical(p.NewData)};
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 3
                    ind = ~cellfun('isempty', dd.Data(:, 2));
                    dd.Data(ind, 3) = {logical(p.NewData)};
                end
            end

            function backend_wrap(t, app, plottypes)
                tmp_data = t.Data(1:end - 1, :);

                INDSmpls = [tmp_data{1:numel(app.DataN.Items), 6}];
                Groups = [tmp_data{INDSmpls,4}];
                smplnms = tmp_data(INDSmpls, 7);
                Volume = [tmp_data{INDSmpls, 5}]; % Recorded Volume of Sample
                Imaged_Volume = [tmp_data{INDSmpls, 8}];
                
                phenos = tmp_data(~cellfun('isempty', tmp_data(:, 2)),2)';
                phenos_cmp = phenos([tmp_data{:,3}]==1);
                phenos = phenos([tmp_data{:,1}]==1);



                Statistics.cellularity_backend(...
                    app, ...
                    smplnms, ...
                    phenos, ...
                    phenos_cmp, ...
                    Groups, ...
                    Volume, ...
                    Imaged_Volume, ...
                    [plottypes.Children.Value] ...
                );
            end
        end

        function cellularity_backend(app, smplnms, phenos, phenos_compare, Groups, Volume, Imaged_Volume, plttp)
            % CELLULARITY_BAKCEND Back-End of Plotting Cellularity. Allows user to visualize
            % sizes of the samples in a form of heatmap. It can take into
            % account volume, in order to show density etc.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - smplnms - cell - Names of samples to plot.
            %   - Groups - array - Same size as smplnms. Groupings of
            %       samples. Allows to treat multiple samples as one,
            %       without needing to merge them etc.
            %   - Volume - array - Same size as smplnms. Volume of each
            %       sample. Allows to look at density of samples, rather
            %       than the absolute values of cell numbers etc.
            
            
% %             plttp(9) % Number of Cells
% %             plttp(8) % Percentage of Cells

% %             plttp(7) % Box Plot by Cell Type
% %             plttp(6) % Heatmap
% %             plttp(5) % Stacked Bar Graph
% %             plttp(4) % Box Plot By Group Number 

% %             plttp(3) % Divide by Median
% %             plttp(2) % Divide by Mean
% %             plttp(1) % Divide by Imaged Area
            %% Initialize the number of cells
            markers = {'o', 's', 'x', '*', '.', '+', 'p'};
            % Just in case
            markers = [markers,markers,markers,markers,markers];
            markers = [markers,markers,markers,markers,markers];
            
            % re-order the samples in the order of groups tag
            [Groups, INDGroups] = sort(Groups);
            smplnms = smplnms(INDGroups);
            Volume = Volume(INDGroups);
            Imaged_Volume = Imaged_Volume(INDGroups);
            
            
            for smplj = 1:numel(smplnms)
                if contains('MetaData', fieldnames(app.data.(smplnms{smplj})))
                    if iscell(Groups(smplj))
                        app.data.(smplnms{smplj}).MetaData.Group = Groups(smplj);
                    else
                        app.data.(smplnms{smplj}).MetaData.Group = {Groups(smplj)};
                    end
                else
                    app.data.(smplnms{smplj}).MetaData = struct;
                    if iscell(Groups(smplj))
                        app.data.(smplnms{smplj}).MetaData.Group = Groups(smplj);
                    else
                        app.data.(smplnms{smplj}).MetaData.Group = {Groups(smplj)};
                    end
                end
            end
            
            
            % Pull the selected cell names and tags
            Phenotypes_vars = phenos;
            for pheno_idx = 1:numel(phenos)
                new_ph = strcmp(phenos{pheno_idx}, app.data.(smplnms{1}).GateTags{2, :});
                Phenotypes_vars{pheno_idx} = app.data.(smplnms{1}).GateTags.Properties.VariableNames{new_ph};
            end
            % This is the gate_XX type name
            Phenotypes = Phenotypes_vars;
            % This is the associated full cell type name
            CellNames = phenos;
            
            % if you are normalizing by the image volume, calculate it
            if plttp(1)
                Vol = zeros(numel(smplnms), 1);
                VolImaged = zeros(numel(smplnms), 1);
                for i=1:numel(smplnms)
                    % Find the Volume of the tissue
                    if Imaged_Volume(i) == 0
                        if contains('ImagedVolume', fieldnames(app.data.(smplnms{i}).MetaData))
                            Vol(i) = app.data.(smplnms{i}).MetaData.ImagedVolume;
                        else
                            Position = table2array(app.data.(smplnms{i}).AllCells(:, {'X', 'Y', 'Z'}));
                            % create a boundary around the tissue
                            [~, Vol(i)] = boundary(Position, 0.5);
                            Vol(i) = Vol(i)/(1000^3); % Convert to mm^3
                            if ~contains('ImagedVolume', fieldnames(app.data.(smplnms{i}).MetaData))
                                app.data.(smplnms{i}).MetaData.ImagedVolume = Vol(i);
                            end
                            if ~contains('SampleVolume', fieldnames(app.data.(smplnms{i}).MetaData))
                                app.data.(smplnms{i}).MetaData.SampleVolume = Volume(i);
                            end
                        end
                    else
                        Vol(i) = Imaged_Volume(i);
                    end
                    VolImaged(i) = Vol(i)/Volume(i);
                end
            end
            %% Load in the selected cell logic vectors
            MFIList= [];
            pheno_weights= [];
            MFI_weights= [];
            DataType= 'Individual Cells';
            DataPrep= 'Cellularity: Number of Cells / Neighborhood';
            NormOpt= 'Sample';
            rmvzeros= 0;
            DataPrepMFI = 'Raw MFI per cell';
            
            [Dat_Pre, Dat_PreALL, INDDatPre, INDzrs, INDons] = Helper.func_loaddata( ...
                app, ...
                smplnms, ...
                phenos, ...
                MFIList, ...
                pheno_weights, ...
                MFI_weights, ...
                DataType, ...
                DataPrep, ...
                NormOpt, ...
                rmvzeros, ...
                DataPrepMFI);
                    
            %% Initialize the matrices
            NCells = table;
            PCells = table;
            DCells = table;
            for j=1:numel(Phenotypes)
                NCells.(Phenotypes{j}) = zeros(numel(smplnms), 1);
                PCells.(Phenotypes{j}) = zeros(numel(smplnms), 1);
                DCells.(Phenotypes{j}) = zeros(numel(smplnms), 1);
            end
            NCells.Sample = cell(numel(smplnms), 1);
            PCells.Sample = cell(numel(smplnms), 1);
            DCells.Sample = cell(numel(smplnms), 1);
            
            %% Calculate the number of cells and percentages
            for smpl_i = 1:numel(smplnms)
                Name = smplnms{smpl_i};
                NCells.Sample{smpl_i} = Name;
                PCells.Sample{smpl_i} = Name;
                DCells.Sample{smpl_i} = Name;

                % Pull a single sample
                dattmp = Dat_Pre(INDDatPre{smpl_i}{1}(1):INDDatPre{smpl_i}{1}(2),:);
                if size(dattmp, 1)>1
                    % Calculate the total number of cells of the types selected
                    PcellsTMP = sum(dattmp, 2);
                    PcellsTMP(PcellsTMP==0) = [];
                    PcellsTMP = numel(PcellsTMP);
% %                     PcellsTMP = sum(PcellsTMP~=0)
                    % Calculate the number of cells of each cell type
                    NcellsTMP = sum(dattmp);
                end
                
                for phn_j=1:numel(Phenotypes)
                    NCells.(Phenotypes{phn_j})(smpl_i) = NcellsTMP(phn_j);
                    if plttp(8)
                        PCells.(Phenotypes{phn_j})(smpl_i) = 100.*NCells.(Phenotypes{phn_j})(smpl_i)/PcellsTMP;
                    end
                    if plttp(1)
                        DCells.(Phenotypes{phn_j})(smpl_i) = NCells.(Phenotypes{phn_j})(smpl_i)./Vol(smpl_i);
                    end
                end
            end
            %% Plot box plots of cells
            if plttp(7)
% %             plttp(3) % Divide by Median
% %             plttp(2) % Divide by Mean
% %             plttp(1) % Divide by Imaged Area
                
                % Divide by Median
                if plttp(3)
                    if plttp(9) 
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./median(Pdata);
                        lbl = 'Number of cells per image / Median';
                        plotBP(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./median(Pdata);
                        lbl = 'Percentage of cells per image / Median';
                        plotBP(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                end
                % Divide by Mean    
                if plttp(2)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./mean(Pdata);
                        lbl = 'Number of cells per image / Mean';
                        plotBP(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./mean(Pdata);
                        lbl = 'Percentage of cells per image / Mean';
                        plotBP(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                end 
                % Divide by Imaged Area
                if plttp(1)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./Vol;
                        lbl = 'Number of cells per image / Image Volume';
                        plotBP(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./Vol;
                        lbl = 'Percentage of cells per image / Image Volume';
                        plotBP(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end  
                end
                % No Normalization
                if ~plttp(3) && ~plttp(2) && ~plttp(1)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        lbl = 'Number of cells per image';
                        plotBP(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        lbl = 'Percentage of cells per image';
                        plotBP(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                end
            end
            
            function plotBP(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                % Plot the cellularity
                fig = figure;
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
                hold on
                if numel(smplnms)>1
                    boxplot(Pdata, 'Color', 'k', 'Symbol', 'w')
                end
                hold on
                for smp_j=1:numel(smplnms)
                    plt = plot(Pdata(smp_j, :), '.');
                    plt.DisplayName = strrep(strrep(smplnms{smp_j}, '_', ' '), 'Sample', '');
                    plt.Marker = markers{Groups(smp_j)};
                    plt.MarkerSize = 5;
                end

                legend show
                ax = gca;
                box off
                ax.Color = 'w';

% % %                 if i==1
                ylabel(lbl);
% % %                 elseif i==2
% % %                     ylabel('Percentage of cells per image');
% % %                 elseif i==3
% % %                     ylabel('Cells per imaged area (um^3)');
% % %                 end

                ax.XTick = 1:numel(CellNames);
                ax.XTickLabel = CellNames;
                xtickangle(45)
            end
            
            
            %% Plot the Heatmap of the number of cells
                
            if plttp(6)
% %             plttp(3) % Divide by Median
% %             plttp(2) % Divide by Mean
% %             plttp(1) % Divide by Imaged Area
                
                % Divide by Median
                if plttp(3)
                    if plttp(9) 
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./median(Pdata);
                        lbl = 'Number of cells per image / Median';
                        plotHM(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./median(Pdata);
                        lbl = 'Percentage of cells per image / Median';
                        plotHM(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                end
                % Divide by Mean    
                if plttp(2)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./mean(Pdata);
                        lbl = 'Number of cells per image / Mean';
                        plotHM(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./mean(Pdata);
                        lbl = 'Percentage of cells per image / Mean';
                        plotHM(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                end 
                % Divide by Imaged Area
                if plttp(1)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./Vol;
                        lbl = 'Number of cells per image / Image Volume';
                        plotHM(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./Vol;
                        lbl = 'Percentage of cells per image / Image Volume';
                        plotHM(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end  
                end
                % No Normalization
                if ~plttp(3) && ~plttp(2) && ~plttp(1)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        lbl = 'Number of cells per image';
                        plotHM(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        lbl = 'Percentage of cells per image';
                        plotHM(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                    end
                end
            end
            function plotHM(app, Pdata, smplnms, Groups, CellNames, markers, lbl)
                % Plot the data heatmap
                fig = figure;
                % add an export data option
                % Create Menu bar
                ExportMenu = uimenu(fig);
                ExportMenu.Text = 'Export';
                % Create an export data option
                ExportPltDat = uimenu(ExportMenu);
                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Corr_Heatmap');
                ExportPltDat.Text = 'Export Plot Data to .csv';
                
                fig.Color = app.GUIOPTS.bgclr;
                fig.InvertHardcopy = 'off';

                hm = heatmap(Pdata');
                hm.Colormap = app.GUIOPTS.redbluecmap;
                hm.YDisplayLabels = CellNames;
                hm.XDisplayLabels = strrep(strrep(smplnms, '_', ' '), 'Sample', '');
                hm.ColorScaling = 'log';
                hm.GridVisible = 'on';
% % %                 if min(min(Pdata))<0
% % %                     lm = [-1.*max(abs(max(max(Pdata))), abs(min(min(Pdata)))),max(abs(max(max(Pdata))), abs(min(min(Pdata))))];
% % %                 else
% % %                     lm = [0,max(abs(max(max(Pdata))), abs(min(min(Pdata))))];
% % %                 end
% % %                 hm.ColorLimits = ln(lm);

                title(lbl);
            end
            %% Compare populations across groups
            
            if plttp(4)
                % Divide by Median
                if plttp(3)
                    if plttp(9) 
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./median(Pdata);
                        lbl = 'Number of cells per image / Median';
                        plotBPC(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./median(Pdata);
                        lbl = 'Percentage of cells per image / Median';
                        plotBPC(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                end
                % Divide by Mean    
                if plttp(2)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./mean(Pdata);
                        lbl = 'Number of cells per image / Mean';
                        plotBPC(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./mean(Pdata);
                        lbl = 'Percentage of cells per image / Mean';
                        plotBPC(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                end 
                % Divide by Imaged Area
                if plttp(1)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./Vol;
                        lbl = 'Number of cells per image / Image Volume';
                        plotBPC(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./Vol;
                        lbl = 'Percentage of cells per image / Image Volume';
                        plotBPC(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end  
                end
                % No Normalization
                if ~plttp(3) && ~plttp(2) && ~plttp(1)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        lbl = 'Number of cells per image';
                        plotBPC(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        lbl = 'Percentage of cells per image';
                        plotBPC(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                end
                
            end
            
            
            % Plot the number of cells and then percentages
            function plotBPC(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                for ph_l=1:numel(Phenotypes)
                    % Plot the cellularity
                    fig = figure;
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
                    hold on
                    limy=0;

                    if numel(smplnms)>1
                        for grp_j=unique(Groups)
                            logic = Groups==grp_j;
                            if sum(logic) >1
                                boxplot(Pdata(logic, ph_l), 'Positions', grp_j.*ones(1, sum(logic)), 'Color', 'k')
                            end
                            for k=1:sum(logic)
                                INDi = find(logic);
                                plt = plot(grp_j, Pdata(INDi(k), ph_l));
                                plt.DisplayName = strrep(strrep(smplnms{INDi(k)}, '_', ' '), 'Sample', '');
                                plt.Marker = markers{grp_j};
                                plt.MarkerSize = 8;
                                limy = max([limy Pdata(INDi(k), ph_l)]);
                                ylim([0 limy])
                            end
                        end
                    end

                    legend show
                    ax = gca;
                    box off
                    ax.Color = 'w';

                    ylabel(lbl);

                    ax.XTick = 1:numel(unique(Groups));
                    ax.XTickLabel = 1:numel(unique(Groups));
                    xlim([min(Groups)-1, max(Groups)+1])
                    xlabel('Group Number');
                    title(CellNames{ph_l})
                end
            end
            
            %% Make stacked bar graphs
            if plttp(5)
                
                % Divide by Median
                if plttp(3)
                    if plttp(9) 
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./median(Pdata);
                        lbl = 'Number of cells per image / Median';
                        plotBG(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./median(Pdata);
                        lbl = 'Percentage of cells per image / Median';
                        plotBG(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                end
                % Divide by Mean    
                if plttp(2)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./mean(Pdata);
                        lbl = 'Number of cells per image / Mean';
                        plotBG(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./mean(Pdata);
                        lbl = 'Percentage of cells per image / Mean';
                        plotBG(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                end 
                % Divide by Imaged Area
                if plttp(1)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        Pdata = Pdata./Vol;
                        lbl = 'Number of cells per image / Image Volume';
                        plotBG(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        Pdata = Pdata./Vol;
                        lbl = 'Percentage of cells per image / Image Volume';
                        plotBG(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end  
                end
                % No Normalization
                if ~plttp(3) && ~plttp(2) && ~plttp(1)
                    if plttp(9)
                        Pdata = table2array(NCells(:, Phenotypes));
                        lbl = 'Number of cells per image';
                        plotBG(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                    if plttp(8)
                        Pdata = table2array(PCells(:, Phenotypes));
                        lbl = 'Percentage of cells per image';
                        plotBG(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                    end
                end
                
            end
            
            % Plot the number of cells and then percentages
            function plotBG(app, Pdata, smplnms, Groups, CellNames, markers, lbl,Phenotypes)
                fig = figure;
                clf
                fig.Color = app.GUIOPTS.bgclr;
                fig.InvertHardcopy = 'off';
                ax = gca;
% %                 ax.Title.String = lbl;
                hold on
                % add an export data option
                % Create Menu bar
                ExportMenu = uimenu(fig);
                ExportMenu.Text = 'Export';
                % Create an export data option
                ExportPltDat = uimenu(ExportMenu);
                ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                ExportPltDat.Text = 'Export Plot Data to .csv';
                % Deal with if only one sample is selected
                if size(Pdata,1)==1
                    plt = bar([Pdata; zeros(1, size(Pdata,2))],1,...
                        'stacked', 'EdgeColor', 'white', 'FaceColor','flat');  
                else
                    plt = bar(Pdata,1,...
                        'stacked', 'EdgeColor', 'white', 'FaceColor','flat');
                end
                % Put in the channel colors
                clrmap = jet(numel(CellNames));
                for n=1:size(Pdata, 2)
                    plt(n).CData = clrmap(n,:);
                    plt(n).DisplayName = CellNames{n};
                end

                box off
                ax.XLabel.String = 'Sample';
                ax.XLim = [0 numel(smplnms)+1];
%                 ax.YLim = [0 1];
                ax.XTick = 0:1:(numel(smplnms)+1);
                ax.XTickLabels = strrep([{''}; smplnms(:); {''}], '_', '');
                ax.XTickLabelRotation = 45;
                ax.YLabel.String = lbl;
                %%%%%
                legend off
                box off
                ax.Color = 'w';
            end
            
            
            
            %% Compare populations across timepoints using tumor volume
% % %             % Find the names of Phenotypes
% % %             Phenotypes = Phenotypes_vars;
% % %             CellNames = phenos_compare;
% % % 
% % %             % Plot the number of cells and then percentages
% % %             for i=1:2
% % %                 if i==1
% % %                     Pdata = table2array(NCells(:, Phenotypes));
% % %                 elseif i==2
% % %                     Pdata = table2array(PCells(:, Phenotypes));
% % %                 end
% % % 
% % %                 for l=1:numel(Phenotypes)
% % %                     % Normalize to the ammount of the volume imaged
% % %                     Pdata(:, l) = Pdata(:, l).*VolImaged;
% % % 
% % %                     % Plot the cellularity
% % %                     fig = figure;
% % %                     % add an export data option
% % %                     % Create Menu bar
% % %                     ExportMenu = uimenu(fig);
% % %                     ExportMenu.Text = 'Export';
% % %                     % Create an export data option
% % %                     ExportPltDat = uimenu(ExportMenu);
% % %                     ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
% % %                     ExportPltDat.Text = 'Export Plot Data to .csv';
% % %                     
% % %                     fig.Color = app.GUIOPTS.bgclr;
% % %                     fig.InvertHardcopy = 'off';
% % %                     hold on
% % %                     limy=0;
% % % 
% % %                     if numel(smplnms)>1
% % %                         for j=unique(Groups)
% % %                             logic = Groups==j;
% % %                             if sum(logic) >1
% % %                                 boxplot(Pdata(logic, l), 'Positions', j.*ones(1, sum(logic)), 'Color', 'k')
% % %                             end
% % %                             for k=1:sum(logic)
% % %                                 INDi = find(logic);
% % %                                 plt = plot(j, Pdata(INDi(k), l));
% % %                                 plt.DisplayName = strrep(strrep(smplnms{INDi(k)}, '_', ' '), 'Sample', '');
% % %                                 plt.Marker = markers{j};
% % %                                 plt.MarkerSize = 5;
% % %                                 limy = max([limy Pdata(INDi(k), l)]);
% % %                                 if limy==0
% % %                                     limy=1;
% % %                                 end
% % %                                 ylim([0 limy])
% % %                             end
% % %                         end
% % %                     end
% % % 
% % %                     ax = gca;
% % %                     box off
% % %                     ax.Color = 'w';
% % % 
% % %                     if i==1
% % %                         ylabel('(Number of cells per image)*(Percentage of tumor imaged)');
% % %                     elseif i==2
% % %                         ylabel('(Percentage of cells per image)*(Percentage of tumor imaged)');
% % %                     end
% % % 
% % %                     ax.XTick = unique(Groups);
% % %                     ax.XTickLabel = unique(Groups);
% % %                     xlim([min(Groups)-1, max(Groups)+1])
% % %                     xlabel('Group Number');
% % %                     title(CellNames{l})
% % %                 end
% % %             end

        end

        function boxplt_func(app, web)
            % BOXPLT_FUNC Front-End of Box Plot function. Allows user to plot
            % distributions of various variables among different
            % combinations of samples and phenotypes.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            if nargin<3
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

            if ~Helper.any_sample(app)
                return;
            end

            %% Pull all of the sample names
            % Load the data set
            new_data = Helper.populate_table(app, 'smpls', {app.DataN.Value}, 'MFI', 'AllCells');

            tData = cell(size(new_data, 1) + 1, 6);

            % Put in the pot options
            tData(1,1) = new_data(1,6);
% % %             tData(1,2) = {'Y-Axis'};

            tData(1:sum(~cellfun('isempty', new_data(:, 3))),2) = {true};
            tData(1:numel(new_data(:, 3)),3) = new_data(:,3);
            tData{end, 2} = false;
            tData{end, 3} = 'Select All';
            
            tData(1:numel(app.DataN.Items),4) = {'None'};
            tData(strcmp(app.DataN.Items, app.DataN.Value), 4) = {'o'};

            tData(1:size(new_data, 1), 5:6) = new_data(:, 7:8);
            tData{end, 5} = false;
            tData{end, 6} = 'Select All';

            % Build the box plot options menus
            UIfig = uifigure('Name', 'Compare Object Chennel Intensities', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end
            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 1500 800];
            ax = uiaxes(UIfig);
            ax.Position = alpha*[800 10 700 790];
            ax.BackgroundColor = app.GUIOPTS.bgclr;
            ax.TickLabelInterpreter = 'none';
            ax.XTickLabelRotation = 45;
            hold(ax, 'on')

            % Create the table of options
            markers = {'None', 'o', 's', 'x', '*', '.', '+', 'p'};
            t = uitable(UIfig);
            t.ColumnFormat = ({new_data(:, 6)', [], [], markers, []});
            t.Data = tData;
            t.Position = alpha*[0 70 790 730];
            t.ColumnName = {'Y-Axis', 'Plot','Phenotypes', 'Marker', 'Include in Plot','Sample'};
            t.ColumnEditable = [true true false true true false];
            t.ColumnWidth = {alpha*150, alpha*80, alpha*150, alpha*100, alpha*80, alpha*150};
            t.CellEditCallback = @(dd, p) smpl_change(app, markers, dd, p);

            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap(t, app, ax));
            btn.Position = alpha*[t.Position(3) - 110, 5, 100, 50];
            btn.Text = 'Plot';
            Helper.func_SetCLR(app, btn, 'button')

            function smpl_change(app, markers, dd, p)
                if p.Indices(2) == 5
                    if p.Indices(1) == find(strcmp(dd.Data(:, 6), 'Select All'))
                        ind = ~cellfun('isempty', dd.Data(:, 5));
                        dd.Data(ind, 5) = {logical(p.NewData)};
                    end
                    ind = ~cellfun('isempty', dd.Data(:, 5));
                    ind(end) = false;
                    ind(ind) = logical(cell2mat(dd.Data(ind, 5)));

                    new_dat = cell(size(dd.Data, 1) - 1, 8);
                    new_dat(:, 1) = {1}; % STUB, Corresponds to weights. Needed for populate table.
                    new_dat(:, 2) = dd.Data(1:end - 1, 2);
                    new_dat(:, 3) = dd.Data(1:end - 1, 3);
                    new_dat(:, 4) = {1}; % STUB, Corresponds to weights. Needed for populate table.
                    new_dat(:, 5) = {false}; % STUB. Corresponds to choices. Needed for populate table.
                    new_dat(:, 7) = dd.Data(1:end - 1, 4);
                    new_dat(:, 8) = dd.Data(1:end - 1, 3);

                    smpls = dd.Data(ind, 6);
                    if isempty(smpls)
                        smpls = {app.DataN.Value};
                    end

                    new_dat = Helper.populate_table(app, ...
                        'smpls', smpls, ...
                        'MFI', 'AllCells', ...
                        'prev_table', new_dat);

                    % Figure out the plotting options (keep if possible)
                    cols = new_dat(~ismember(new_dat(:, 6), Constants.ignore_names), 6);
                    dd.ColumnFormat = ({cols', [], [], markers, []});
                    prev_plt_opt = dd.Data(:, 1);
                    non_empty = ~cellfun('isempty', prev_plt_opt);
                    keep_same = ismember(prev_plt_opt(non_empty), new_dat(:, 6));
                    prev_plt_opt(~keep_same) = new_dat(1:sum(~keep_same), 6);

                    % Update the table
                    dd.Data(:, 1) = cell(numel(dd.Data(:, 1)), 1);
                    dd.Data(1:numel(prev_plt_opt), 1) = prev_plt_opt;

                    dd.Data(1:end - 1, 2) = cell(size(dd.Data, 1) - 1, 1);
                    dd.Data(1:numel(new_dat(:, 2)), 2) = new_dat(:, 2);

                    dd.Data(1:end - 1, 3) = cell(size(dd.Data, 1) - 1, 1);
                    dd.Data(1:numel(new_dat(:, 3)), 3) = new_dat(:, 3);

                    dd.Data(1:end - 1, 5) = cell(size(dd.Data, 1) - 1, 1);
                    dd.Data(1:numel(new_dat(:, 4)), 5) = new_dat(:, 7);

                    dd.Data(1:end - 1, 6) = cell(size(dd.Data, 1) - 1, 1);
                    dd.Data(1:numel(new_dat(:, 5)), 6) = new_dat(:, 8);
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 3), 'Select All')) && p.Indices(2) == 2
                    ind = ~cellfun('isempty', dd.Data(:, 3));
                    dd.Data(ind, 2) = {logical(p.NewData)};
                end
            end
            
            function backend_wrap(t, app, ax)
                tmp_data = t.Data(1:end - 1, :);
                smplnms = tmp_data([tmp_data{:, 5}]'==1, 6);
                pnms = tmp_data([tmp_data{:, 2}]'==1, 3);

                Statistics.boxplt_backend( ...
                    app, ...
                    smplnms, ...
                    pnms, ...
                    tmp_data{1, 1}, ... y label (what to make boxplot of)
                    ax ...
                );
            end
        end
        
        function boxplt_backend(app, smplnms, pnms, ynm, ax)
            % BOXPLT_BACKEND Back-End of Box Plot function. Allows user to plot
            % distributions of various variables among different
            % combinations of samples and phenotypes.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - smplnms - cell - Names of samples to be plotted. Will be
            %       combined with phenotypes to make smplnms * pnms number
            %       of things on x-axis.
            %   - pnms - cell - Names of phenotypes to be plotted. Will be
            %       combined with samples to make smplnms * pnms number of
            %       things on x-axis.
            %   - ax - axis - Axis on which to plot box plots.
            %   - ynm - string - Thing of which to make a box-plot. It's
            %       going to be visualized on plot as y-axis label.
            
            cla(ax)

            % Import the data
            ax.YLabel.String = ynm;

            ynm = Helper.valid_channel(ynm);
            if ismember(ynm, Constants.ignore_names)
                ynm = Helper.valid_var(ynm);
            end

            ax.XTick = 1:numel(pnms);

            size_pltdat = 0;

            for smpl_idx=1:numel(smplnms)
                smpl = smplnms{smpl_idx};
                size_pltdat = size_pltdat + (numel(app.data.(smpl).AllCells(:, 1)) * numel(pnms));
            end

            pltdat = zeros(1, size_pltdat);
            group = zeros(1, size_pltdat);
            index_group = 1;

            for smpl_idx=1:numel(smplnms)
                smpl = smplnms{smpl_idx};
                for pnm_i=1:numel(pnms)
                    plt_tmp = app.data.(smpl).AllCells(:, ynm);
                    indexes = app.data.(smpl).AllCells{:, Helper.gate_full2tag(app, pnms{pnm_i}, smpl)};
                    plt_tmp = plt_tmp(indexes==1, :);
                    pltdat(index_group: index_group + numel(plt_tmp) - 1) = table2array(plt_tmp)';
                    group(index_group: index_group + numel(plt_tmp) - 1) = pnm_i .* ones(1, numel(plt_tmp));
                    index_group = index_group + numel(plt_tmp);
                end
            end

            non_zero = group ~= 0;
            group = group(non_zero);
            pltdat = pltdat(non_zero);
            boxplot(ax, pltdat, group)
            ax.XLim = [0 numel(pnms)+1];
            ax.XTickLabel = Helper.full_gate(pnms);
        end

        function CellInteraction_func(app, web)
            % CELLINTERACTION_FUNC Front-End of Cell Interaction function. Allows user to plot
            % different statistics as different variations of a heatmap in
            % between different cell populations/cell regions in different
            % samples/groups. It is useful to understand how a
            % population/region relates to global average, or other
            % populations.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            if ~Helper.any_sample(app) || ~Helper.any_net(app)
                return;
            end

            % Build the clustering options menus
            UIfig = uifigure('Name', 'Cell-Cell Correlation Heatmap', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end
            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 1365 800];

            % Initialize the table of options
            t = uitable(UIfig);
            t.Position = alpha*[0 70 1360 730];

            % Select Neighborhood type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Neighborhood Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10 45 175 15];
            DataType = uidropdown(UIfig);
            DataType.Items = fieldnames(app.net);
            if ~isempty(DataType.Items)
                DataType.Value = DataType.Items{1};
            end
            DataType.Position = alpha*[10 10 175 30];
            DataType.BackgroundColor = app.GUIOPTS.bgclr;
            p = struct('Indices', [], 'EditData', 1, 'NewData', []);
            DataType.ValueChangedFcn = @(~, ~) switched_sample(app, DataType, t, p);

            %Build the table
            [MFI_TYPE, sample] = Helper.find_MFI(app);
            if isempty(MFI_TYPE)
                errordlg("In order to run this operation at least one of your samples had to be scanned beforehand.");
                return;
            end
            tmpData = Helper.populate_table(app, 'smpls', {sample}, 'MFI', MFI_TYPE, 'fill_checkbox', false);

            tData = cell(max(app.net.(DataType.Value).NR+1, size(tmpData, 1)) + 1, 10);
            % Cells
            tData(1:numel(tmpData(:, 2)), 1)= tmpData(:, 2);
            tData(1:numel(tmpData(:, 3)), 2) = tmpData(:, 3);
            tData(end, 1:2) = {false, 'Select All'};
            % MFI
            tData(1:numel(tmpData(:, 5)), 3) = tmpData(:, 5);
            tData(1:numel(tmpData(:, 6)), 4) = tmpData(:, 6);
            tData(end, 3:4) = {false, 'Select All'};
            % Sample
            tData(1:numel(tmpData(:, 7)), 5) = tmpData(:, 7);
            tData(1:sum(~cellfun(@isempty,tmpData(:, 8))), 6) = {1};
            tData(1:numel(tmpData(:, 8)), 7) = tmpData(:, 8);
            tData(end, [5, 7]) = {false, 'Select All'};
            % Regions
            tData(1, 8) = {false};
            tData(2:app.net.(DataType.Value).NR+1, 8) = {true};
            tData(2:app.net.(DataType.Value).NR+1, 9) = {1};
            tData(1:app.net.(DataType.Value).NR+1, 10) = num2cell(0:app.net.(DataType.Value).NR);

            % Populate the table of options
            t.Data = tData;
            t.ColumnName = { ...
                'Include', 'Phenotype', ...
                'Include', 'Channel MFI', ...
                'Plot HM for', 'Group',  'Sample', ...
                'Include', 'Group', 'Region' ...
            };
            t.ColumnEditable = [true false true false true true false true true false];
            t.ColumnWidth = {alpha*75, alpha*350, alpha*75, alpha*250, alpha*75, alpha*50, alpha*225, alpha*75, alpha*50, alpha*100};
            t.CellEditCallback = @(dd, p) switched_sample(app, DataType, dd, p);

            % Select Data Preperation
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Data Preperation:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175 45 175 15];
            DataPrep = uidropdown(UIfig);
            DataPrep.Items = {'Composition: Number of Cells / Total cells in Neighborhood', ...
                              'Global Composition: Number of Cells / Max Cells in Tissue Neighborhoods', ...
                              'Cellularity: Number of Cells / Neighborhood', ...
                              'Density: Number of Cells / Volume (Area if 2D) Per Neighborhood', ...
                              'Binary: If cell is in neighborhood = 1', ...
                              'Standardize: subtract Mean, divide by standard deviation', ...
                              'Corrected Density: Number of Cells / Volume (Area if 2D) of Neighborhood Inside Tissue', ...
                              'Composition+Density: Number of Cells / Total cells in Neighborhood; Density of Neighborhood'};
            DataPrep.Value = {'Cellularity: Number of Cells / Neighborhood'};
            DataPrep.Position = alpha*[10+175 10 175 30];
            DataPrep.BackgroundColor = app.GUIOPTS.bgclr;
            
            % Create a Norm per sample or per dataset label
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Normalize per:');
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175 50+5 75 20];
            % Create a Norm per sample or per dataset button
            NormPer = uibuttongroup(UIfig, 'Visible','off');
            NormPer.Position = alpha*[10+175+175, 10, 75, 45];
            Helper.func_SetCLR(app, NormPer, 'UICpopup')
            % Create two radio buttons in the button group.
            r1 = uiradiobutton(NormPer, 'Text','Sample',...
                'Position',alpha*[5, 25, 75, 20]);
            Helper.func_SetCLR(app, r1, 'table')
            r2 = uiradiobutton(NormPer, 'Text','Dataset',...
                'Position',alpha*[5, 5, 75, 20]);
            Helper.func_SetCLR(app, r2, 'table')
            NormPer.Visible = 'on';

            % Select Heatmap type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Heatmap Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+75+10 45 400 15];
            MAPType = uidropdown(UIfig);
            MAPType.Items = {'Individual Heatmap for each Sample', ...
                              'Combined Heatmap of all Samples'};
            MAPType.Value = {'Individual Heatmap for each Sample'};
            MAPType.Position = alpha*[10+175+175+75+10 10 175 30];
            MAPType.BackgroundColor = app.GUIOPTS.bgclr;

            % Select Color Scale
            lbl = uilabel(UIfig); lbl.Text = sprintf('Color Scale:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+175+75+10 45 100 15];
            MAPScale = uidropdown(UIfig);
            MAPScale.Items = {'linear', ...
                              'log'};
            MAPScale.Value = {'linear'};
            MAPScale.Position = alpha*[10+175+175+175+75+10 10 100 30];
            MAPScale.BackgroundColor = app.GUIOPTS.bgclr;

            % Select Calculation Type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Calculation:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+175+100+75+10 45 100 15];
            CalcType = uidropdown(UIfig);
            CalcType.Items = {'Pearson Correlation Coefficient', ...
                              'Spearman Correlation Coefficient', ...
                              'Kendall Correlation Coefficient', ...
                              'Covariance', ...
                              'Cross-correlation', ...
                              'Conditional Phenotype/Channel value given region', ...
                              'Conditional Phenotype/Channel value given Phenotype/Channel', ...
                              'Pair plots of raw data'};
            CalcType.Value = CalcType.Items{1};
            CalcType.Position = alpha*[10+175+175+175+100+75+10 10 175 30];
            CalcType.BackgroundColor = app.GUIOPTS.bgclr;
            
            % Select Data transform 
            lbl = uilabel(UIfig); lbl.Text = sprintf('Transform:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+175+100+175+75+10 45 100 15];
            DatScale = uidropdown(UIfig);
            DatScale.Items = {'none', ...
                              'log'};
            DatScale.Value = {'none'};
            DatScale.Position = alpha*[10+175+175+175+100+175+75+10 10 100 30];
            Helper.func_SetCLR(app, lbl, 'button')
            
            % Select Confidence Interval
            lbl = uilabel(UIfig); lbl.Text = sprintf('Confidence Interval:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+175+100+175+75+10+100 45 100 15];
            ConfInt = uieditfield(UIfig, 'numeric', ...
                'Limits', [0 1],...
                'LowerLimitInclusive', 'off', ...
                'UpperLimitInclusive', 'on' ...
            );
            ConfInt.Value = 1;
            ConfInt.Position = alpha*[10+175+175+175+100+175+75+10+100 10 100 30];
            ConfInt.Editable = 'on';
            Helper.func_SetCLR(app, ConfInt, 'UIField')

            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap( ...
                t, app, DataType, DataPrep, MAPType, MAPScale, CalcType, DatScale, ConfInt, NormPer ...
            ));
            btn.Position = alpha*[1365-100, 5, 100, 50];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')

            function switched_sample(app, DataType, dd, p)
                % if the samples index was changed
                if isempty(p.Indices) || p.Indices(2) == 5
                    % Make sure that at least one thing is selected
                    if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, p.Indices(2))), p.Indices(2))))
                        dd.Data{p.Indices(1), p.Indices(2)} = true;
                        return;
                    end
                    % Process Select All button
                    fill_select_all = false;
                    if ~isempty(p.Indices) && p.Indices(1) == find(strcmp(dd.Data(:, 7), 'Select All'))
                        ind = ~cellfun('isempty', dd.Data(:, p.Indices(2)));
                        dd.Data(ind, p.Indices(2)) = {logical(p.NewData)};
                        if ~logical(p.NewData)
                            dd.Data{strcmp(dd.Data(:, 7), app.DataN.Value), p.Indices(2)} = true;
                        end
                        fill_select_all = logical(p.NewData);
                    end
                    if strcmp(app.net.(DataType.Value).userdata.DataType, 'Cell Centered Neighborhoods')
                        scan_type = 'MFICCN';
                    else
                        scan_type = 'MFIRSN';
                    end
                    ind = ~cellfun('isempty', dd.Data(:, 5));
                    ind(end) = false;
                    ind(ind) = logical(cell2mat(dd.Data(ind, 5)));

                    %Find which samples were selected
                    fcsz = ~cellfun('isempty', dd.Data(:, 2));
                    fcsz(end) = false;
                    tmpDat = cell(size(dd.Data, 1) - 1, 6);
                    tmpDat(fcsz, 1) = {1};
                    tmpDat(fcsz, 2) = dd.Data(fcsz, 1);
                    tmpDat(:, 3) = dd.Data(1:end - 1, 2);

                    fcsz = ~cellfun('isempty', dd.Data(:, 4));
                    fcsz(end) = false;
                    tmpDat(fcsz, 4) = {1};
                    tmpDat(fcsz, 5) = dd.Data(fcsz, 3);
                    tmpDat(:, 6) = dd.Data(1:end - 1, 4);
                    tmpDat(:, 7) = dd.Data(1:end - 1, 5);
                    tmpDat(:, 8) = dd.Data(1:end - 1, 7);

                    smpls = dd.Data(ind, 7);
                    if isempty(smpls)
                        smpls = {app.DataN.Value};
                    end

                    tmpDat = Helper.populate_table(app, ...
                        'smpls', smpls, ...
                        'mfi', scan_type, ...
                        'prev_table', tmpDat ...
                    );

                    % Keep column which is not conserved in populate table
                    keep_6_col = dd.Data(:, 6);
                    prev_smpls = dd.Data(:, 7);
                    keep_regions = dd.Data(:, 8:end);

                    dd_ph = ~cellfun('isempty', dd.Data(:, 1));
                    dd_ph(end) = false;
                    dd_ph = dd.Data(dd_ph, 2);
                    dd_mfi = ~cellfun('isempty', dd.Data(:, 3));
                    dd_mfi(end) = false;
                    dd_mfi = dd.Data(dd_mfi, 4);

                    tmp_ph = ~cellfun('isempty', tmpDat(:, 2));
                    tmp_ph = tmpDat(tmp_ph, 3);
                    tmp_mfi = ~cellfun('isempty', tmpDat(:, 5));
                    tmp_mfi = tmpDat(tmp_mfi, 6);
                    if isempty(p.Indices) && app.net.(DataType.Value).NR+1 < size(dd.Data, 1)
                        dd.Data(:, 8:end) = [];
                        dd.Data(1, 8) = {false};
                        dd.Data(2:app.net.(DataType.Value).NR+1, 8) = {true};
                        dd.Data(2:app.net.(DataType.Value).NR+1, 9) = {1};
                        dd.Data(1:app.net.(DataType.Value).NR+1, 10) = num2cell(0:app.net.(DataType.Value).NR);
                    elseif ~Helper.setequal(dd_ph, tmp_ph) || ~Helper.setequal(dd_mfi, tmp_mfi) || isempty(p.Indices)
                        dd.Data = cell(max(app.net.(DataType.Value).NR + 1, size(tmpDat, 1)) + 1, 10);
                        % Cells
                        fcsz = ~cellfun('isempty', tmpDat(:, 3));
                        logic_scale = false(size(dd.Data, 1), 1);
                        logic_scale(1:size(fcsz, 1)) = fcsz;
                        dd.Data(logic_scale, 1) = tmpDat(fcsz, 2);
                        dd.Data(1:numel(tmpDat(:, 3)), 2) = tmpDat(:, 3);
                        dd.Data(end, 1:2) = {false, 'Select All'};
                        % MFI
                        fcsz = ~cellfun('isempty', tmpDat(:, 6));
                        logic_scale = false(size(dd.Data, 1), 1);
                        logic_scale(1:size(fcsz, 1)) = fcsz;
                        dd.Data(logic_scale, 3) = tmpDat(fcsz, 5);
                        dd.Data(1:numel(tmpDat(:, 6)), 4) = tmpDat(:, 6);
                        dd.Data(end, 3:4) = {false, 'Select All'};
                        % Sample
                        dd.Data(1:numel(tmpDat(:, 7)), 5) = tmpDat(:, 7);
                        if sum(~cellfun('isempty', keep_6_col)) == sum(~cellfun('isempty', tmpDat(:, 8)))
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
                        dd.Data(1:numel(tmpDat(:, 8)), 7) = tmpDat(:, 8);
                        dd.Data(end, [5, 7]) = {fill_select_all, 'Select All'};
                        % Regions
                        if ~isempty(p.Indices)
                            if size(keep_regions, 1) > size(dd.Data, 1)
                                dd.Data(:, 8:end) = keep_regions(1:size(dd.Data,1), :);
                            else
                                dd.Data(1:size(keep_regions, 1), 8:end) = keep_regions;
                            end
                        else
                            dd.Data(1,8) = {false};
                            dd.Data(2:app.net.(DataType.Value).NR+1, 8) = {true};
                            dd.Data(2:app.net.(DataType.Value).NR+1, 9) = {1};
                            dd.Data(1:app.net.(DataType.Value).NR+1, 10) = num2cell(0:app.net.(DataType.Value).NR);
                        end
                    end
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
                    ind = ~cellfun('isempty', dd.Data(:, p.Indices(2)));
                    dd.Data(ind, p.Indices(2)) = {logical(p.NewData)};
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 4), 'Select All')) && p.Indices(2) == 3
                    ind = ~cellfun('isempty', dd.Data(:, p.Indices(2)));
                    dd.Data(ind, p.Indices(2)) = {logical(p.NewData)};
                end
            end

            function backend_wrap(t, app, DataType, DataPrep, MAPType, MAPScale, CalcType, DatScale, ConfInt, NormPer)
                tmp_data = t.Data(1:end - 1, :);
                Statistics.CellInteraction_backend( ...
                    app, ...
                    cell2mat(tmp_data([tmp_data{:, 8}]==1, 9)'), ... Plotting Groups
                    cell2mat(tmp_data([tmp_data{:, 8}]==1, 10)'), ... NRegions
                    tmp_data([tmp_data{:, 5}]'==1, 7), ... Samples
                    cell2mat(tmp_data([tmp_data{:, 5}]==1, 6)'), ... Sample Groups
                    tmp_data([tmp_data{:, 1}]'==1, 2), ... Phenotypes
                    tmp_data([tmp_data{:, 3}]'==1, 4), ... MFIs
                    DataType.Value, ...
                    DataPrep.Value, ...
                    MAPType.Value, ...
                    MAPScale.Value, ...
                    CalcType.Value, ...
                    DatScale.Value, ...
                    ConfInt.Value, ...
                    NormPer.SelectedObject.Text...
                );
            end
        end

        function CellInteraction_backend(app, groups, NRegions, smplnms, smp_groups, phenos, MFIList, DataType, DataPrep, MAPType, MAPScale, CalcType, DatScale, ConfInt, NormOpt)
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
            
            % Check for valid input
            if isempty(phenos) && isempty(MFIList)
                errordlg('Please include any phenotypes or MFIs to make correlation from.');
                return;
            end
            % Make a waitbar
            vPD = waitbar(0, 'Plotting heatmaps', ...
                'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));
            vPD.Position(4) = 85;
            %% Load the data
            if ~isempty(MFIList)
                MFIListSUB = Helper.valid_channel(MFIList);
            else
                MFIListSUB = [];
            end
            [Dat_PreMAIN, ~, INDDatPre, ~, ~] = Helper.func_loaddata(app, ...
                smplnms, ...
                phenos, ...
                MFIListSUB, ...
                [], ...
                [], ...
                app.net.(DataType).userdata.DataType, ...
                DataPrep, ...
                NormOpt, ...
                0);
            %% Plot stuff
            smp_Ngroups = numel(unique(smp_groups));
            for Sgroup_i = 1:smp_Ngroups
                smplnms_i = smplnms(ismember(smp_groups,Sgroup_i));
                smplsIND = find(ismember(smp_groups,Sgroup_i));
                for i=1:numel(smplnms_i)
                    %% Prepare the data the data
                    if ~isvalid(vPD)
                        return;
                    end
                    waitbar(i/numel(smplnms_i),vPD, 'PROCESSING!');
                    
                    CellNms = phenos;

                    % Load the Region data
                    if strcmp(app.net.(DataType).userdata.DataType, 'Raster Scanned Neighborhoods')
                        ROW2 = app.data.(smplnms_i{i}).MFIRSN.([Constants.other_tag, DataType]);
                    elseif strcmp(app.net.(DataType).userdata.DataType, 'Cell Centered Neighborhoods')
                        ROW2 = app.data.(smplnms_i{i}).MFICCN.([Constants.other_tag, DataType]);
                    elseif strcmp(app.net.(DataType).userdata.DataType, 'Individual Cells')
                        ROW2 = app.data.(smplnms_i{i}).AllCells.([Constants.other_tag, DataType]);
                    end

                    if strcmp(DataPrep, 'Composition+Density: Number of Cells / Total cells in Neighborhood; Density of Neighborhood')
                        CellNms = [CellNms, 'Total Cellular Density'];
                    end

                    % Load the MFI data
                    CellNms = [CellNms', MFIList'];

                    %% Combine all data Samples from specific groups if selected
                    switch MAPType
                        case 'Individual Heatmap for each Sample'
                            Dat_Pre = Dat_PreMAIN(INDDatPre{smplsIND(i)}{1}(1):INDDatPre{smplsIND(i)}{1}(2),:);
                        case 'Combined Heatmap of all Samples'
                            if i==1
                                CombinedDat = Dat_PreMAIN(INDDatPre{smplsIND(i)}{1}(1):INDDatPre{smplsIND(i)}{1}(2),:);
                                CombinedROW2 = ROW2;
                            else
                                CombinedDat = [CombinedDat; Dat_PreMAIN(INDDatPre{smplsIND(i)}{1}(1):INDDatPre{smplsIND(i)}{1}(2),:)];
                                CombinedROW2 = [CombinedROW2; ROW2];
                            end
                            if i==numel(smplnms_i)
                                Dat_Pre = CombinedDat;
                                ROW2 = CombinedROW2;
                            end
                    end

                    % Only include selected regions
                    Ngroups = numel(unique(groups));

                    %% Exclude data outside the confidence interval
                    if ConfInt~=1
                        % Assuming a Normal Distribution of Noise
                        % 1-ConfInt
                        % find the standard deviation
                        mu = mean(Dat_Pre);
                        sigma = cov(Dat_Pre);
                        prob = mvncdf(Dat_Pre, mu, sigma);
                        keep_idxs = prob > 1 - ConfInt & prob < ConfInt;
                        if sum(keep_idxs) == 0
                            errordlg('Confidence Interval to Restrictive. Please try with higher value.');
                            return;
                        end
                        Dat_Pre = Dat_Pre(keep_idxs, :);
                        ROW2 = ROW2(keep_idxs);
                    end
                    %% Change the data scale
                    switch DatScale
                        case 'log'
                            Dat_Pre = log(Dat_Pre);
                        case 'linear'
                    end
                    %% Do the Calculations if it is the last sample or you
                    % are making individual plots
                    if strcmp(MAPType, 'Individual Heatmap for each Sample') || i==numel(smplnms_i)
                        for group_i = 1:Ngroups
                            NRegions_i = NRegions(ismember(groups,group_i));
                            pcst_i = Dat_Pre(ismember(ROW2, NRegions_i), :);
                            ROW2_i = ROW2(ismember(ROW2, NRegions_i), :);
                            YLabel = CellNms;
                            switch CalcType
                                case 'Pearson Correlation Coefficient'
                                    % Find the spatial correlation between cell-types
                                    [RatioMatrix, P_Matrix] = corr(pcst_i, 'Type', 'Pearson');
                                    ttla = 'Pearson Correlation coefficients between cells';
                                    ttlb = 'p-values of Pearson Correlation';
                                case 'Spearman Correlation Coefficient'
                                    % Find the spatial correlation between cell-types
                                    [RatioMatrix, P_Matrix] = corr(pcst_i, 'Type', 'Spearman');
                                    ttla = 'Spearman Correlation coefficients between cells';
                                    ttlb = 'p-values of Spearman Correlation';
                                case 'Kendall Correlation Coefficient'
                                    % Find the spatial correlation between cell-types
                                    [RatioMatrix, P_Matrix] = corr(pcst_i, 'Type', 'Kendall');
                                    ttla = 'Kendall Correlation coefficients between cells';
                                    ttlb = 'p-values of Kendall Correlation';
                                case 'Covariance'
                                    % Find the spatial correlation between cell-types
                                    RatioMatrix = cov(pcst_i);
                                    P_Matrix = 0.*RatioMatrix+1;
                                    ttla = 'Covariance matrix between neighborhoods';
                                case 'Cross-correlation'
                                    % Find the spatial correlation between cell-types
                                    RatioMatrix = xcorr(pcst_i);
                                    P_Matrix = 0.*RatioMatrix+1;
                                    ttla = 'Cross-correlation matrix between neighborhoods';
                                case 'Conditional Phenotype/Channel value given region'
                                    % Phenotypes are first in the pcst
                                    global_mean = mean(pcst_i);
                                    RatioMatrix = zeros(numel(NRegions_i) + 1, size(pcst_i, 2));
                                    YLabel = cell(1, size(RatioMatrix, 1));
                                    YLabel(1) = {'Mean of Selected Regions'};
                                    RatioMatrix(1, :) = global_mean;
                                    for reg_idx=1:numel(NRegions_i)
                                        reg_i = NRegions_i(reg_idx);
                                        mean_reg_i = mean(pcst_i(ROW2_i == reg_i, :));
                                        RatioMatrix(reg_idx + 1, :) = mean_reg_i - global_mean;
                                        YLabel(reg_idx + 1) = {['Region ' char(num2str(reg_i))]};
                                    end
                                    P_Matrix = 0.*RatioMatrix+1;
                                    ttla = 'Conditional Phenotype/Channel value given region in relation to global mean';
                                case 'Conditional Phenotype/Channel value given Phenotype/Channel'
                                    % Phenotypes are first in the pcst
                                    global_mean = mean(pcst_i);
                                    RatioMatrix = zeros(size(pcst_i, 2));
                                    % loop through the cell types
                                    for pcst_idx=1:size(pcst_i, 2)
                                        % find the normalized difference from the global mean for cell_type_i
                                        mean_pcst_idx = (pcst_i(:, pcst_idx)- global_mean(pcst_idx))./global_mean(pcst_idx);
                                        % Multiply the fold change for cell_i by the number of cellj
                                        mean_pcst_idx = mean_pcst_idx.*pcst_i;
                                        % The 'interaction' is the the average of this value
                                        RatioMatrix(pcst_idx, :) = mean(mean_pcst_idx);
                                    end
                                    P_Matrix = 0.*RatioMatrix+1;
                                    ttla = 'mean[ (Fold change of cell_i(row) * number of cells_j(column) per neighborhood]';
                                case 'Pair plots of raw data'
                                    region = cell(numel(ROW2_i), 1);
                                    for k = 1:numel(ROW2_i)
                                        region{k} = ['R' num2str(ROW2_i(k))];
                                    end
                                    fig = figure;
                                    fig.Color = app.GUIOPTS.bgclr;
                                    fig.InvertHardcopy = 'off';
                                    clf
                                    hold on
                                    pairplot(pcst_i, YLabel, region, app.net.(DataType).cmap((unique(ROW2_i, 'stable')+1), :))
                            end

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

                                fig.Color = app.GUIOPTS.bgclr;
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

                                limits = [-1, 1];
                                hm.ColorLimits = limits;

                                if strcmp(MAPType, 'Individual Heatmap for each Sample')
                                    title(sprintf([ttla '\n' strrep(smplnms_i{i}, '_', ' ') '\n Regions: ' num2str(NRegions_i)]))
                                elseif strcmp(MAPType, 'Combined Heatmap of all Samples')
                                    if smp_Ngroups==1
                                        title(sprintf([ttla ' \n All Samples' '\n Regions: ' num2str(NRegions_i)]))
                                    else
                                        title(sprintf([ttla ' \n Group ' num2str(Sgroup_i) ' Samples' '\n Regions: ' num2str(NRegions_i)]))
                                    end
                                end

                                % Plot the p-Values, if there are any
                                if contains(CalcType, 'Correlation Coefficient')
                                    % Plot the p-values heatmap
                                    fig = figure;
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
                                        title(sprintf([ttlb '\n' strrep(smplnms_i{i}, '_', ' ') '\n Regions: ' num2str(NRegions_i)]))
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
                                    % Plot the difference between the correlation heatmaps
                                    fig = figure;
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
                                        title(sprintf(['Difference of the ' ttla 'between groups\n' strrep(smplnms_i{i}, '_', ' ')]))
                                    elseif strcmp(MAPType, 'Combined Heatmap of all Samples')
                                        title(sprintf(['Difference of the ' ttla ' between groups\n All Samples']))
                                    end
                                end

                                % if the samples have groups
                                if Sgroup_i == 1 && smp_Ngroups > 1
                                    RatioMatrix_G1 = RatioMatrix;
                                end
                                if Sgroup_i == smp_Ngroups && smp_Ngroups > 1
                                    % Plot the difference between the correlation heatmaps
                                    fig = figure;
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
                                        title(sprintf(['Difference of the ' ttla 'between groups\n' strrep(smplnms_i{i}, '_', ' ')]))
                                    elseif strcmp(MAPType, 'Combined Heatmap of all Samples')
                                        title(sprintf(['Difference of the ' ttla ' between groups\n All Samples']))
                                    end
                                end
                            end % end if not pair plots
                        end % End if plot
                    end % end region groups loop

                end
            end
            Helper.func_closeVPD(app, vPD)

            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end
        end

        function func_Pseudospace(app, web)
            % FUNC_PSEUDOSPACE Front-End of Pseudospce Plot function. It allows user to
            % reduce and sort data in 1-D to show how neighborhoods
            % interact with each other. It is useful for showing relations
            % between cells/neighborhood on human-readable scale.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            
            %% Intro function handling
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

            if ~Helper.any_sample(app)
                return;
            end

            [MFI_TYPE, smpls] = Helper.find_MFI(app);
            if isempty(MFI_TYPE)
                errordlg("In order to run this operation at least one of your samples had to be scanned beforehand.");
                return;
            end

            if ~iscell(smpls)
                smpls = {smpls};
            end

            %% Pull all of the fieldnames from the MFIRSN table
            dispnames = app.data.(smpls{1}).(MFI_TYPE).Properties.VariableNames';
            dispnames = [dispnames(startsWith(dispnames, Constants.gate_tag)); dispnames(~startsWith(dispnames, Constants.gate_tag))];
            dispnames(startsWith(dispnames, Constants.gate_tag)) = Helper.gate_tag2full(app, dispnames(startsWith(dispnames, Constants.gate_tag)), smpls{1});
            dispnames = Helper.full_channel(dispnames);
            tData = cell(max(numel(smpls), numel(dispnames)) + 1, 8);

            %% Populate initial table
            % Plot (T/F)
            tData(1:numel(dispnames), 1)= {false};
            tData(1, 1)= {true};
            % Use in Sort (T/F)
            tData(1:numel(dispnames), 2)= {false};
            tData(1, 2)= {true};
            % All Variable Names From Either CCN or RSN
            tData(1:numel(dispnames), 3) = dispnames;
            tData(end, 1:3) = {false, false, 'Select All'};
            % Order that neighborhood parameters will be sorted in
            tData(1:numel(dispnames), 4) = {1};
            % Weight to use in sorting data
            tData(1:numel(dispnames), 5) = {1};
            % Ammount to smooth each cell/MFI parameter
            tData(1:numel(dispnames), 6) = {5000};
            % Include Sample (T/F)
            tData(1:numel(fieldnames(app.data)), 7) = num2cell(strcmp(fieldnames(app.data), smpls));
            % List of sample names
            tData(1:numel(fieldnames(app.data)), 8) = fieldnames(app.data);
            tData(end, 7:8) = {false, 'Select All'};

            % 'Order (if same num use in same sort, if 1,2,3 sort one then 2 then 3)'
            % Weight; use positive numbers to put on right of graph and negative to put on left of graph

            %% Define additional fields around table
            % Build the clustering options menus
            UIfig = uifigure('Name', 'Pseudospace Options', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end
            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 1000 800];
            t = uitable(UIfig);

            % Select Data Preperation
% %             lbl = uilabel(UIfig); lbl.Text = sprintf('Select Data Preperation:');
% %             Helper.func_SetCLR(app, lbl, 'label')
% %             lbl.Position = alpha*[10+175+75 45 400 15];

            DataPrep = uidropdown(UIfig);
            DataPrep.Items = Constants.neigh_cell_norm;
            DataPrep.Value = {'Composition: Number of Cells / Total cells in Neighborhood'};
            DataPrep.Position = alpha*[10+175+75 10 175 30];
            DataPrep.BackgroundColor = app.GUIOPTS.bgclr;
            
            % Select Data Preperation for MFI
            DataPrepMFI = uidropdown(UIfig);
            DataPrepMFI.Items = Constants.neigh_mfi_norm;
            DataPrepMFI.Value = {'Sum MFI per neighborhood'};
            DataPrepMFI.Position = alpha*[10+175+75 40 175 30];
            DataPrepMFI.BackgroundColor = app.GUIOPTS.bgclr;

            % Create a Norm per sample or per dataset label
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Normalize per:');
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175 50+5 75 20];
            % Create a Norm per sample or per dataset button
            NormPer = uibuttongroup(UIfig, 'Visible','off');
            NormPer.Position = alpha*[10+175, 10, 75, 45];
            Helper.func_SetCLR(app, NormPer, 'UICpopup')
            % Create two radio buttons in the button group.
            r1 = uiradiobutton(NormPer, 'Text','Sample',...
                'Position',alpha*[5, 25, 75, 20]);
            Helper.func_SetCLR(app, r1, 'table')
            r2 = uiradiobutton(NormPer, 'Text','Dataset',...
                'Position',alpha*[5, 5, 75, 20]);
            Helper.func_SetCLR(app, r2, 'table')
            NormPer.Visible = 'on';
            
            % Select Neighborhood type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Neighborhood Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10 45 175 15];
            DataType = uidropdown(UIfig);
            DataType.Items = {'Raster Scanned Neighborhoods', ...
                              'Cell Centered Neighborhoods', ...
                              'Individual Cells'};
            DataType.Value = {'Raster Scanned Neighborhoods'};
            DataType.Position = alpha*[10 10 175 30];
            DataType.BackgroundColor = app.GUIOPTS.bgclr;

            % Select Heatmap type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Plot Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+75+10 45 175 15];
            MAPType = uidropdown(UIfig);
            MAPType.Items = {'Individual Plot for each Sample', ...
                              'Combined Plot of all Samples'};
            MAPType.Value = {'Individual Plot for each Sample'};
            MAPType.Position = alpha*[10+175+175+75+10 10 175 30];
            MAPType.BackgroundColor = app.GUIOPTS.bgclr;

            % Select Color Scale
            lbl = uilabel(UIfig); lbl.Text = sprintf('Scale:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+175+75+10+175 45 100 15];

            MAPScale = uidropdown(UIfig);
            MAPScale.Items = {'linear', ...
                              'log'};
            MAPScale.Value = {'linear'};
            MAPScale.Position = alpha*[10+175+175+75+10+175 10 100 30];
            MAPScale.BackgroundColor = app.GUIOPTS.bgclr;

            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap(t, app, DataType, DataPrep, MAPType, MAPScale,NormPer,DataPrepMFI));
            btn.Position = alpha*[1000-100, 5, 100, 50];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')


            %% Creation of the table GUI
            t.Data = tData;
            t.Position = alpha*[0 70 1000 730];
            t.ColumnName = {'Plot (T/F)','Use in Sort (T/F)', 'Variable Names', 'Order', 'Weight', 'Smoothing', 'Include (T/F)', 'Sample List'};

            t.ColumnEditable = [true true false true true true true false];
            t.ColumnWidth = {alpha*75, alpha*100, alpha*200, alpha*75, alpha*75, alpha*75, alpha*75, alpha*200};
            t.CellEditCallback = @(dd, p) switched_sample(app, DataType, dd, p);

            function switched_sample(app, DataType, dd, p)
                if p.Indices(2) == 7
                    ind = ~cellfun('isempty', dd.Data(:, 7));
                    ind(end) = false;
                    ind(p.Indices(1)) = logical(p.PreviousData);
                    ind(ind) = logical(cell2mat(dd.Data(ind, 7)));
                    old_samples = dd.Data(ind, 8);

                    % Make sure that at least one thing is selected
                    if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, p.Indices(2))), p.Indices(2))))
                        dd.Data{p.Indices(1), p.Indices(2)} = true;
                        return;
                    end
                    % Process Select All button
                    fill_select_all = false;
                    if p.Indices(1) == find(strcmp(dd.Data(:, 8), 'Select All'))
                        ind = ~cellfun('isempty', dd.Data(:, p.Indices(2)));
                        dd.Data(ind, p.Indices(2)) = {logical(p.NewData)};
                        if ~logical(p.NewData)
                            dd.Data{strcmp(dd.Data(:, 8), app.DataN.Value), p.Indices(2)} = true;
                        end
                        fill_select_all = logical(p.NewData);
                    end
                    if strcmp(DataType.Value, 'Cell Centered Neighborhoods')
                        scan_type = 'MFICCN';
                    else
                        scan_type = 'MFIRSN';
                    end

                    % Get current and old samples.
                    ind = ~cellfun('isempty', dd.Data(:, 7));
                    ind(end) = false;
                    ind(ind) = logical(cell2mat(dd.Data(ind, 7)));
                    samples = dd.Data(ind, 8);

                    % Get old phenotypes shorts
                    old_short = {};
                    for old_smpl_idx=1:numel(old_samples)
                        [~, old_short_i] = Helper.get_gates(app, old_samples{old_smpl_idx});
                        old_short = union(old_short, old_short_i,'stable');
                    end

                    % Group column into gates and everything else (so it corresponds to populate table format)
                    ind = dd.Data(:, 2);
                    ind = ~cellfun('isempty', ind);
                    ind(end) = false;
                    gate_idx = ismember(dd.Data(:, 3), old_short);
                    gates = dd.Data(gate_idx, 3);
                    rest = dd.Data(~gate_idx & ind, 3);

                    % Populate table
                    % first cols are used for changing between checkboxes and '1's (again compatibility with populate table)
                    tmpDat = cell(max(numel(fieldnames(app.data), max(numel(gates), numel(rest)))), 8);
                    first_col = dd.Data(gate_idx, 4);
                    first_col = double(cell2mat(first_col));
                    first_col = num2cell(first_col);
                    tmpDat(1:numel(gates), 1) = first_col;
                    tmpDat(1:numel(gates), 2) = dd.Data(gate_idx, 2);
                    tmpDat(1:numel(gates), 3) = gates;

                    first_col = dd.Data(~gate_idx & ind, 4);
                    first_col = double(cell2mat(first_col));
                    first_col = num2cell(first_col);
                    tmpDat(1:numel(rest), 4) = first_col;
                    tmpDat(1:numel(rest), 5) = dd.Data(~gate_idx & ind, 2);
                    tmpDat(1:numel(rest), 6) = rest;

                    tmpDat(1:numel(fieldnames(app.data)), 7) = dd.Data(1:numel(fieldnames(app.data)), 7);
                    tmpDat(1:numel(fieldnames(app.data)), 8) = fieldnames(app.data);

                    % Keep Order, Weight and Smoothing columns in place.
                    keys = dd.Data(ind, 3);
                    order = dd.Data(ind, 4);
                    weight = dd.Data(ind, 5);
                    smoothing = dd.Data(ind, 6);

                    % Make new table with current samples
                    tmpDat = Helper.populate_table(app, ...
                        'smpls', samples, ...
                        'MFI', scan_type, ...
                        'prev_table', tmpDat ...
                    );

                    dd_ph = ~cellfun('isempty', dd.Data(:, 1));
                    dd_ph(end) = false;
                    dd_ph = dd.Data(dd_ph, 3);

                    tmp_ph = ~cellfun('isempty', tmpDat(:, 2));
                    tmp_ph = tmpDat(tmp_ph, 3);
                    tmp_mfi = ~cellfun('isempty', tmpDat(:, 5));
                    tmp_mfi = tmpDat(tmp_mfi, 6);

                    if ~Helper.setequal(dd_ph, union(tmp_ph, tmp_mfi,'stable'))
                        idx_gate = ~cellfun('isempty', tmpDat(:, 1));
                        idx_rest = ~cellfun('isempty', tmpDat(:, 4));
                        newDat = cell(max(sum(idx_gate) + sum(idx_rest), numel(fieldnames(app.data))) + 1, 8);
                        newDat(1:sum(idx_gate), 1) = tmpDat(idx_gate, 1);
                        newDat(1:sum(idx_gate), 2) = tmpDat(idx_gate, 2);
                        newDat(1:sum(idx_gate), 3) = tmpDat(idx_gate, 3);
                        newDat(end, 1:3) = {0, false, 'Select All'};
                        newDat(sum(idx_gate) + 1:sum(idx_gate) + sum(idx_rest), 1) = tmpDat(idx_rest, 4);
                        newDat(sum(idx_gate) + 1:sum(idx_gate) + sum(idx_rest), 2) = tmpDat(idx_rest, 5);
                        newDat(sum(idx_gate) + 1:sum(idx_gate) + sum(idx_rest), 3) = tmpDat(idx_rest, 6);

                        first_col_idx = ~cellfun('isempty', newDat(:, 1));
                        first_col = newDat(first_col_idx, 1);
                        first_col = logical(cell2mat(first_col));
                        first_col = num2cell(first_col);
                        newDat(first_col_idx, 1) = first_col;

                        for row_idx=1:size(newDat, 1) - 1
                            if isempty(newDat(row_idx, 1))
                                continue;
                            end
                            idx = find(strcmp(keys, newDat(row_idx, 3)));
                            if any(idx)
                                newDat(row_idx, 4) = order(idx);
                                newDat(row_idx, 5) = weight(idx);
                                newDat(row_idx, 6) = smoothing(idx);
                            else
                                newDat(row_idx, 4) = {1};
                                newDat(row_idx, 5) = {1};
                                newDat(row_idx, 6) = {5000};
                            end
                        end

                        newDat(1:numel(fieldnames(app.data)), 7) = dd.Data(1:numel(fieldnames(app.data)), 7);
                        newDat(1:numel(fieldnames(app.data)), 8) = dd.Data(1:numel(fieldnames(app.data)), 8);
                        newDat(end, 7:8) = {fill_select_all, 'Select All'};
                        dd.Data = newDat;
                    end
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 3), 'Select All')) && p.Indices(2) <= 2
                    ind = ~cellfun('isempty', dd.Data(:, p.Indices(2)));
                    dd.Data(ind, p.Indices(2)) = {logical(p.NewData)};
                end
            end
            
            function backend_wrap(t, app, DataType, DataPrep, MAPType, MAPScale,NormPer,DataPrepMFI)
                tmp_data = t.Data(1:end - 1, :);
                
                Statistics.Pseudospace_backend( ...
                    app, ...
                    tmp_data([tmp_data{:, 7}]'==1, 8), ... Samples
                    tmp_data([tmp_data{:, 1}]'==1, 3), ... Phenotypes Plot
                    tmp_data([tmp_data{:, 2}]'==1, 3), ... Phenotypes Sort
                    DataType.Value, ...
                    DataPrep.Value, ...
                    MAPType.Value, ...
                    MAPScale.Value, ...
                    [tmp_data{[tmp_data{:, 2}]'==1, 4}], ... order
                    [tmp_data{[tmp_data{:, 2}]'==1, 5}], ... weight
                    [tmp_data{[tmp_data{:, 1}]'==1, 6}], ... smoothing
                    NormPer.SelectedObject.Text, ...
                    DataPrepMFI.Value... 
                );
            end
        end
        
        function Pseudospace_backend(app, smplnms, phenos_plot, phenos_sort, DataType, DataPrep, MAPType, MAPScale, order, weight, smoothing, NormOpt, DataPrepMFI)
            % PSEUDOSPACE_BACKEND Back-End of Pseudospce Plot function. It allows user to
            % reduce and sort data in 1-D to show how neighborhoods
            % interact with each other. It is useful for showing relations
            % between cells/neighborhood on human-readable scale.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - smplnms - cell - Samples from which data to make
            %       pseudospace plot will be pulled.
            %   - phenos_plot - cell - Phenotypes which will be plotted on
            %       pseudospace.
            %   - phenos_sort - cell - Phenotypes on which pseudospace will
            %       be sorted. Those are the ones, which are used to
            %       distinguish between the regions.
            %   - DataType - string - Name of the Model to use in creating
            %       the pseudospace. Has to exist as field under app.net.
            %   - DataPrep - string - How to normalize and prepare the
            %       data. Same as in Helper.Func_DataPrep.
            %   - MAPType - string - Whether to make one `Combined Heatmap
            %       of all Samples` or multiple `Individual Heatmap for
            %       each Sample`.
            %   - MAPScale - string - Whether to apply apply `log` or
            %       `linear` scale to pseudospace in the end.
            %   - order - array - Same size as phenos_sort. Order of sorted
            %       phenotypes. Smaller the number, more important the
            %       phontype will be.
            %   - weight - array - Same size as phenos_sort. Second-order
            %       importance metric for phenotype. Will be multiplied
            %       before sort to all the phenotype. Higher the number
            %       more important the phenotype. Mainly used with negative
            %       numbers (if not default), to revert the sort order.
            %   - smoothing - array - Same size as phenos_plot. How much a
            %       specific phenotype should be smoothened. It increases
            %       readability of the plot, while decrease the detail of
            %       it.
            
            
            % Initialize some variables
            INDDatPrePlot = cell(numel(smplnms), 1);
            INDDatPreSort = cell(numel(smplnms), 1);
            
            for i=1:numel(smplnms)
                
                %% Load the data

                % Make a waitbar
                vPD = waitbar(0, 'Preparing Data', ...
                    'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));
                vPD.Position(4) = 85;
                
                dispnamesFull = phenos_plot;

                PlotNames_i = Helper.gate_full2tag(app, phenos_plot, smplnms{i});
                Plot_idx = ~startsWith(PlotNames_i, Constants.other_tag) & ~startsWith(PlotNames_i, Constants.neigh_tag);
                PlotNames_i(Plot_idx) = Helper.valid_channel(PlotNames_i(Plot_idx));
                PlotNames_i(ismember(PlotNames_i, Constants.ignore_names)) = Helper.valid_var(PlotNames_i(ismember(PlotNames_i, Constants.ignore_names)));
                                
                SortNames_i = Helper.gate_full2tag(app, phenos_sort, smplnms{i});
                Sort_idx = ~startsWith(SortNames_i, Constants.other_tag) & ~startsWith(SortNames_i, Constants.neigh_tag);
                SortNames_i(Sort_idx) = Helper.valid_channel(SortNames_i(Sort_idx));
                SortNames_i(ismember(SortNames_i, Constants.ignore_names)) = Helper.valid_var(SortNames_i(ismember(SortNames_i, Constants.ignore_names)));
                
                if i==1
                    if isempty(SortNames_i)
                        errordlg("Nothing was found in the 'Use in Sort' options. Select at least one parameter to sort the neighborhoods with.");
                    end
                end
                
                % Load the Cell data
                if strcmp(DataType, 'Raster Scanned Neighborhoods')
                    Dat_PreALL = app.data.(smplnms{i}).MFIRSN;
                    PlotDat_i = table2array(app.data.(smplnms{i}).MFIRSN(:, PlotNames_i));
                    SortDat_i = table2array(app.data.(smplnms{i}).MFIRSN(:, SortNames_i));
                elseif strcmp(DataType, 'Cell Centered Neighborhoods')
                    Dat_PreALL = app.data.(smplnms{i}).MFICCN;
                    PlotDat_i = table2array(app.data.(smplnms{i}).MFICCN(:, PlotNames_i));
                    SortDat_i = table2array(app.data.(smplnms{i}).MFICCN(:, SortNames_i));
                elseif strcmp(DataType, 'Individual Cells')
                    DataPrep = 'Cellularity: Number of Cells / Neighborhood';
                    Dat_PreALL = app.data.(smplnms{i}).AllCells;
                    PlotDat_i = table2array(app.data.(smplnms{i}).AllCells(:, PlotNames_i));
                    SortDat_i = table2array(app.data.(smplnms{i}).AllCells(:, SortNames_i));
                end
                
                % remove any cells, or other extranious data from Dat_PreAll                
                IND = contains(Dat_PreALL.Properties.VariableNames, Constants.gate_tag);
                Dat_PreALL(:,IND) = [];
                IND = contains(Dat_PreALL.Properties.VariableNames, Constants.channel_tag);
                Dat_PreALL(:,IND) = [];
                IND = contains(Dat_PreALL.Properties.VariableNames, Constants.neigh_tag);
                Dat_PreALL(:,IND) = [];
                IND = contains(Dat_PreALL.Properties.VariableNames, Constants.other_tag);
                Dat_PreALL(:,IND) = [];

                if strcmp(NormOpt, 'Sample')
                    
                    % Prepare the data as selected by the user
                    PlotDat_i = Helper.Func_DataPrep(Dat_PreALL, PlotDat_i, DataPrep);
                    SortDat_i = Helper.Func_DataPrep(Dat_PreALL, SortDat_i, DataPrep);
                end

                if i==1
                    INDDatPrePlot{1} = {[1, size(PlotDat_i, 1)], smplnms{1}};
                    INDDatPreSort{1} = {[1, size(SortDat_i, 1)], smplnms{1}};
                    CombinedDatPLT = PlotDat_i;
                    CombinedDatSRT = SortDat_i;
                    CombinedDatALL = Dat_PreALL;
                else
                    INDDatPrePlot{i} = {[size(CombinedDatPLT, 1)+1, size(CombinedDatPLT, 1)+size(PlotDat_i, 1)], smplnms{i}};
                    INDDatPreSort{i} = {[size(CombinedDatSRT, 1)+1, size(CombinedDatSRT, 1)+size(SortDat_i, 1)], smplnms{i}};
                    CombinedDatPLT = [CombinedDatPLT; PlotDat_i];
                    CombinedDatSRT = [CombinedDatSRT; SortDat_i];
                    CombinedDatALL = [CombinedDatALL; Dat_PreALL];
                end
                if i==numel(smplnms) && strcmp(NormOpt, 'Dataset')
                    % Prepare the data as selected by the user
                    CombinedDatPLT = Helper.Func_DataPrep(CombinedDatALL, CombinedDatPLT, DataPrep);
                    CombinedDatSRT = Helper.Func_DataPrep(CombinedDatALL, CombinedDatSRT, DataPrep);
                end
                waitbar(1,vPD, 'Done!');
                Helper.func_closeVPD(app, vPD)
            end
            for i=1:numel(smplnms)
                % Make a waitbar
                vPD = waitbar(0, 'Plotting Data', ...
                    'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));
                vPD.Position(4) = 85;
                
                % Do the plot if its the last sample or you are making
                % plots for every sample
                if strcmp(MAPType, 'Individual Plot for each Sample') || i==numel(smplnms)
                    switch MAPType
                        case 'Individual Plot for each Sample'
                            PlotDat = CombinedDatPLT(INDDatPrePlot{i}{1}(1):INDDatPrePlot{i}{1}(2),:);
                            SortDat = CombinedDatSRT(INDDatPreSort{i}{1}(1):INDDatPreSort{i}{1}(2),:);
                        case 'Combined Plot of all Samples'
                            PlotDat = CombinedDatPLT;
                            SortDat = CombinedDatSRT;
                    end
                
                    if strcmp(MAPScale, 'log')
                        PlotDat = real(log(PlotDat));
                        SortDat = real(log(SortDat));
                    end
                    
                    if strcmp(DataPrep, 'Composition+Density: Number of Cells / Total cells in Neighborhood; Density of Neighborhood')
                        dispnamesFull = [dispnamesFull, 'Total Cellular Density'];
                    end
                    % remove any zero neighborhoods from the plot
                    IND_Zero = any(PlotDat.*(PlotDat>10e-3),2); % I hate rounding errors... can't just ask ==0
                    PlotDat = PlotDat(IND_Zero,:);
                    SortDat = SortDat(IND_Zero,:);
                    % Don't sort the neighborhoods if nothing was selected
                    if ~isempty(SortNames_i)
                        for sort_i = unique(order)
                            if size(SortDat,1)>1
                                sortdat_i = sum(weight(order==sort_i).*SortDat(:,order==sort_i),2);
                            else
                                sortdat_i = weight(order==sort_i).*SortDat(:,order==sort_i);
                            end

                            [~, IND_i] = sort(sortdat_i);
                            PlotDat = PlotDat(IND_i,:);
                            SortDat = SortDat(IND_i,:);
                        end
                    end
                    % TODO; find a way to put the indices of the sorted neighborhoods back into the neighborhoods table

                    for ph_i = 1:size(PlotDat,2)
                        PlotDat(:,ph_i) = smoothdata(PlotDat(:,ph_i),'gaussian', smoothing(ph_i));
                    end

                    % Noralize the lines
                    PlotDat = PlotDat./max(PlotDat);

                    % Plot the data heatmap
                    fig = figure;
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
                    clf
                    hold on

                    % add an export data option
                    % Create Menu bar
                    ExportMenu = uimenu(fig);
                    ExportMenu.Text = 'Export';

                    % Create an export data option
                    ExportPltDat = uimenu(ExportMenu);
                    ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                    ExportPltDat.Text = 'Export Plot Data to .csv';

                    for Ph_i=1:size(PlotDat,2)
                        plot(PlotDat(:,Ph_i), '.','DisplayName', dispnamesFull{Ph_i})
                    end

                    axis tight
                    ax = gca;
                    ax.XColor = app.GUIOPTS.txtclr;
                    ax.YColor = app.GUIOPTS.txtclr;
                    ax.XLim = [0 size(PlotDat, 1)];
                    ax.YDir = 'Normal';
                    ax.XDir = 'Normal';
                    ax.XLabel.String = 'Pseudo-space';
                    ax.YLabel.String = DataPrep;

                    if strcmp(MAPType, 'Individual Plot for each Sample')
                        if strcmp(DataType, 'Individual Cells')
                            title(sprintf(['Sorted Cells \n' strrep(smplnms{i}, '_', ' ')]))
                        else
                            title(sprintf(['Sorted Neighborhoods \n' strrep(smplnms{i}, '_', ' ')]))
                        end
                    elseif strcmp(MAPType, 'Combined Plot of all Samples')
                        if strcmp(DataType, 'Individual Cells')
                            title(sprintf('Sorted Cells \n All Samples'))
                        else
                            title(sprintf('Sorted Neighborhoods \n All Samples'))
                        end
                    end

                    box off
                    hold on
                end % End if plot

                waitbar(1,vPD, 'Done!');
                Helper.func_closeVPD(app, vPD)
            end
            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end
        end
        
    end
end
