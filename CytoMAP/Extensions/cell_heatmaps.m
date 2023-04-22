function cell_heatmaps(app)
    %% Build a user interface
    % HEATMAP_FUNC Front-End of making heatmap function. Allows user to create
    % heatmap of Phenotypes/Channels within specific samples
    % (either seperate or combined throughout samples).
    %
    % Input:
    %   - app - Instance of CytoMAP
    
    alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

% % %     if ~Helper.any_sample(app) || ~Helper.any_net(app)
% % %         return;
% % %     end

    [~, sample] = Helper.find_MFI(app);
    
    tmpData = Helper.populate_table(app, 'MFI', 'AllCells', 'smpls', {sample}, 'fill_checkbox', false);

    tData = cell(size(tmpData, 1) + 1, 6);
    tData(1:size(tmpData, 1), 1)= tmpData(:, 2);
    tData(1:size(tmpData, 1), 2) = tmpData(:, 3);

    tData(1:size(tmpData, 1), 3) = tmpData(:, 5);
    tData(1:size(tmpData, 1), 4) = tmpData(:, 6);

    tData(1:size(tmpData, 1), 5) = tmpData(:, 7);
    tData(1:size(tmpData, 1), 6) = tmpData(:, 8);
    tData(end, :) = {false, 'Select All', false, 'Select All', false, 'Select All'}; 

    % Build the main figure
    UIfig = uifigure('Name', 'Plot Cell MFI Heatmap Options', 'Scrollable', 'on');
    Helper.func_SetCLR(app, UIfig, 'UIfigure');
    UIfig.Position = alpha*[10 10 1000 800];

    % Select model type
    lbl = uilabel(UIfig); lbl.Text = sprintf('Select What To Compare:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = alpha*[10+175+150 45 400 15];
    DataType = uidropdown(UIfig);
    DataType.Items = [fieldnames(app.net); {'Phenotypes'}];
% % %     DataType.Items = fieldnames(app.net);
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
    t.CellEditCallback = @(dd, p) switched_sample(app, dd, p);

    % Select Data Preperation for MFI
    DataPrepMFI = uidropdown(UIfig);
    DataPrepMFI.Items = Constants.cell_mfi_norm;
    DataPrepMFI.Value = DataPrepMFI.Items{1};
    DataPrepMFI.Position = alpha*[10+75 10 175 30];
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
    MAPScale.Value = {'linear'};
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
    btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap(t, app, DataType, MAPType, MAPScale, NormPer, DataPrepMFI));
    btn.Position = alpha*[1000-150, 10, 100, 50];
    btn.Text = 'Ok';
    Helper.func_SetCLR(app, btn, 'button')

    function switched_sample(app, dd, p)
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
                'MFI', 'AllCells', ...
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

    function backend_wrap(t, app, DataType, MAPType, MAPScale, NormPer, DataPrepMFI)
        tmp_data = t.Data(1:end-1, :);
        heatmap_backend( ...
            app, ...
            tmp_data([tmp_data{:, 5}]'==1, 6), ... Samples
            tmp_data([tmp_data{:, 1}]'==1, 2), ... Phenotypes
            tmp_data([tmp_data{:, 3}]'==1, 4), ... MFIs
            DataType.Value, ...
            MAPType.Value, ...
            MAPScale.Value, ...
            NormPer.SelectedObject.Text, ...
            DataPrepMFI.Value ...
        );
    end

    function heatmap_backend(app, smplnms, phenos, MFIList, DataType, MAPType, MAPScale, NormOpt, DataPrepMFI, web)
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
            
            
            if ~strcmp(DataType, 'Phenotypes')
                is_in = ones(numel(smplnms), 1);
                for smpl_idx=1:numel(smplnms)
                    if ~ismember(...
                            [Constants.other_tag, DataType], ...
                            app.data.(smplnms{smpl_idx}).AllCells.Properties.VariableNames...
                    )
                        is_in(smpl_idx) = 0;
                    end
                end
                is_in = logical(is_in);
                if ~all(is_in)
                    if any(is_in)
                        warndlg( ...
                            "Samples: " + join(smplnms(~is_in), ", ") + ...
                            " do not contain the chosen model." + newline, ...
                            "However samples: " + join(smplnms(is_in), ",") + ...
                            " contain this model, and process will continue with only these samples."...
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
            end
            %% Load the data                
            rmvzeros = 1;
            
            if ~isempty(MFIList)
                MFInms = MFIList;
                MFIList = Helper.valid_channel(MFIList);
            end
            [Dat_PreMAIN, ~, INDDatPre, INDzrs, INDons] = Helper.func_loaddata( ...
                app, ...
                smplnms, ...
                phenos, ...
                MFIList, ...
                [], ...
                [], ...
                'Individual Cells', ...
                'Cellularity: Number of Cells / Neighborhood', ...
                NormOpt, ...
                rmvzeros, ...
                DataPrepMFI);

% % %                 type = 'AllCells';
                if strcmp(DataType, 'Phenotypes')
                    % Pull the logicals
                    Dat_PreLOGICMain = Dat_PreMAIN(:, 1:numel(phenos));
                end
                % remove the binary cell type classifiers from the data table
                Dat_PreMAIN = Dat_PreMAIN(:, (numel(phenos)+1):end);
            
            switch MAPType
                case 'Combined Heatmap of all Samples'
                    Dat_Pre = Dat_PreMAIN;
                    if strcmp(DataType, 'Phenotypes')
                        Dat_PreLOGIC = Dat_PreLOGICMain;
                    end
            end
            
            %% Plot the heatmap
            for i=1:numel(smplnms)
                
                % Make a waitbar
                vPD = waitbar(0, 'Plotting heatmap', ...
                    'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));

                switch MAPType
                    case 'Individual Heatmap for each Sample'
                        Dat_Pre = Dat_PreMAIN(INDDatPre{i}{1}(1):INDDatPre{i}{1}(2),:);
                        if strcmp(DataType, 'Phenotypes')
                            Dat_PreLOGIC = Dat_PreLOGICMain(INDDatPre{i}{1}(1):INDDatPre{i}{1}(2),:);
                        end
                end

                if ~strcmp(DataType, 'Phenotypes')
                    % Load the Model data
                    if rmvzeros==1
                        ROW2 = app.data.(smplnms{i}).AllCells.([Constants.other_tag, DataType])(INDons{i},:);
                    else
                        ROW2 = app.data.(smplnms{i}).AllCells.([Constants.other_tag, DataType]);
                    end
                end
% % %                 RowMAP = app.net.(DataType).cmap;

                % Load the MFI data
                switch MAPType
                    case 'Combined Heatmap of all Samples'
                        
                        if ~strcmp(DataType, 'Phenotypes')
                            if i==1
                                CombinedROW2 = ROW2;
                            else
                                CombinedROW2 = [CombinedROW2; ROW2];
                            end
                        end
                        
                        if i==numel(smplnms)
                            if strcmp(MAPScale, 'log')
                                Dat_Pre = real(log(Dat_Pre));
                            end
                            
                            if ~strcmp(DataType, 'Phenotypes')
                                [CombinedROW2, IND2] = sort(CombinedROW2);
                                Dat_Pre = [CombinedROW2, CombinedROW2, Dat_Pre(IND2, :)];
                                Dat_Pre(CombinedROW2==0, :) = []; 
                                ROW2 = CombinedROW2;
                            else
                                Dat_Pre = [zeros(size(Dat_Pre, 1), 1), Dat_Pre]; %add empty column for cell percentages
                            end
                        end
                    case 'Individual Heatmap for each Sample'
                        % This will only import the data that was used by the NN to
                        % cluster the neighborhoods, it should also sort them into the
                        % same order as was used in the NN
                        if strcmp(MAPScale, 'log')
                            Dat_Pre = real(log(Dat_Pre));
                        end
                        
                        if ~strcmp(DataType, 'Phenotypes')
                            %Sort the data
                            [ROW2, IND] = sort(ROW2);
                            Dat_Pre = [ROW2, ROW2, Dat_Pre(IND, :)]; %add information of color
                            Dat_Pre(ROW2==0, :) = []; % remove neighborhoods with no cells
                        else
                            Dat_Pre = [zeros(size(Dat_Pre, 1), 1), Dat_Pre]; %add empty column for cell percentages
                        end
                end

                % Do the plot if its the last sample or you are making
                % plots for every sample
                if strcmp(MAPType, 'Individual Heatmap for each Sample') || i == numel(smplnms)
                    % Plot the data heatmap
                    fig = figure;
                    fig.Color = app.GUIOPTS.bgclr;
                    fig.InvertHardcopy = 'off';
                    if web==1
                        fig.Visible = 'off';
                    end
                    clf
                    % add an export data option
                    % Create Menu bar
                    ExportMenu = uimenu(fig);
                    ExportMenu.Text = 'Export';

                    % Create an export data option
                    ExportPltDat = uimenu(ExportMenu);
                    ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Corr_Heatmap');
                    ExportPltDat.Text = 'Export Plot Data to .csv';
                    
                    % average over the groups (clusters or cells)
                    if strcmp(DataType, 'Phenotypes')
                        % Cells
                        ncell = numel(phenos);
                        Npts = size(Dat_Pre,1);
                        for cel_i = 1:ncell
                                N_sub = sum(Dat_PreLOGIC(:, cel_i));
                                Dat_Pre(Dat_PreLOGIC(:, cel_i)==1, 1) = 100.*(N_sub/Npts);
                            if size(Dat_Pre(Dat_PreLOGIC(:, cel_i)==1, :), 1) ==0
                                % there are no cells in this group
                                regmean = zeros(1, size(Dat_Pre, 2));
                                Dat_PreLOGIC = [zeros(1, size(Dat_PreLOGIC, 2)); Dat_PreLOGIC];
                                Dat_Pre = [regmean; Dat_Pre];
                            elseif size(Dat_Pre(Dat_PreLOGIC(:, cel_i)==1, :), 1) >1
                                regmean = mean(Dat_Pre(Dat_PreLOGIC(:, cel_i)==1, :));
                                Dat_PreLOGIC = [zeros(1, size(Dat_PreLOGIC, 2)); Dat_PreLOGIC];
                                Dat_Pre = [regmean; Dat_Pre];
                            else
                                Dat_PreLOGIC = [zeros(1, size(Dat_PreLOGIC, 2)); Dat_PreLOGIC];
                                Dat_Pre = [Dat_Pre(Dat_PreLOGIC(:, cel_i)==1, :); Dat_Pre];
                            end
                        end
                        Dat_Pre = Dat_Pre(1:ncell, :);

                    else
                        % Clusters or regions
                        regs = unique(Dat_Pre(:,1), 'stable');
                        nreg = numel(regs);
                        Npts = size(Dat_Pre,1);
                        for reg_i = 1:nreg
                                N_sub = sum(Dat_Pre(:,1)==regs(reg_i));
                                Dat_Pre(Dat_Pre(:,1)==regs(reg_i), 2) = 100.*(N_sub/Npts);
                            if size(Dat_Pre(Dat_Pre(:,1)==regs(reg_i), :), 1) >1
                                regmean = mean(Dat_Pre(Dat_Pre(:,1)==regs(reg_i), :));
                                Dat_Pre(Dat_Pre(:,1)==regs(reg_i), :) = [];
                                Dat_Pre = [regmean; Dat_Pre];
                            end
                        end
                    end
                    % Create the image
% %                     imagesc(Dat_Pre);
                    if strcmp(DataType, 'Phenotypes')
                        ax = heatmap(Dat_Pre(end:-1:1, :));
                    else
                        ax = heatmap(Dat_Pre(:, 2:end));
                    end

                    cmap = app.GUIOPTS.redbluecmap;
                    cmap(1:round(size(cmap,1)/2), :) = [];
                    ax.Colormap = cmap;
                    
                    
                    
                    if strcmp(DataType, 'Phenotypes')
                        ax.XDisplayLabels =[{'% of Cells in Phenotype'}; MFInms];
                        ax.YDisplayLabels =  strrep(phenos', '_', ' ');
                        
                        ax.YLabel = 'Phenotype';
                    else
                        if strcmp(app.net.(DataType).userdata.DataType, 'Raster Scanned Neighborhoods') || ...
                                strcmp(app.net.(DataType).userdata.DataType, 'Cell Centered Neighborhoods')
                            ax.XDisplayLabels =[{'% of Cells in Region'}; MFInms];
                            ax.YDisplayLabels =  strcat('R', num2str(Dat_Pre(:, 1)));
                            ax.YLabel = 'Region';
                        elseif strcmp(app.net.(DataType).userdata.DataType, 'Individual Cells') 
                            ax.XDisplayLabels =[{'% of Cells in Cluster'}; MFInms];
                            ax.YDisplayLabels =  strcat('C', num2str(Dat_Pre(:, 1)));
                            ax.YLabel = 'Cluster';
                        end
                        % add annotations if you got em
                        if ismember('RegNames', fieldnames(app.net.(DataType)))   
                            ax.YDisplayLabels = app.net.(DataType).RegNames(Dat_Pre(:, 1));
                        end     
                    end
                    
                    ax.XLabel = 'Channel';
                    
                    if strcmp(DataPrepMFI, 'Standardize: subtract Mean, divide by standard deviation')
                        limits = [min(min(Dat_Pre(:, 3:end))), max(max(Dat_Pre(:, 3:end)))]-1;
                        [~,MAX] = max(abs(limits));
                        [~,MIN] = min(abs(limits));
                        limits(MIN) = -1*limits(MAX); % Make the color limits symetrical
                        colormap(app.GUIOPTS.redbluecmap)
                        caxis(limits);
                        cttl = 'Standardized MFI';
                    else
                        limits = [min(min(Dat_Pre(:, 3:end))), max(max(Dat_Pre(:, 3:end)))];
                        if limits(1) == -Inf
                            limits(1) = -6;
                        end
                        if limits(2) > 0
                            limits(1) = 0;
                        end
                        caxis(limits);
                        cttl = DataPrepMFI;
                    end
                    strtmp = strcat(phenos, ',');
                    if strcmp(MAPType, 'Individual Heatmap for each Sample')
                        title(sprintf(['Cellular Mean Fluorescent Intensity \n' ...
                            strrep(smplnms{i}, '_', ' ') '\n' ...
                            'Cells: ' strtmp{:}]))
                    elseif strcmp(MAPType, 'Combined Heatmap of all Samples')
                        title(sprintf(['Cellular Mean Fluorescent Intensity \n ' ...
                            'All Samples \n' ...
                            'Cells: ' strtmp{:}]))
                    end
                    
                    % put a title on the colorbar
                    axs = struct(gca); %ignore warning that this should be avoided
                    cb = axs.Colorbar;
                    cb.Label.String = cttl;
                    
                end % End if plot

                waitbar(1,vPD, 'Done!');
                close(vPD)
            end

            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end

        end
   
    
end