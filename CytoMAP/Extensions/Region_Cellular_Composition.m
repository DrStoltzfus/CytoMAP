function Region_Cellular_Composition(app)
    %% Build a user interface
    % HEATMAP_FUNC Front-End of making heatmap function. Allows user to create
    % heatmap of Phenotypes/Channels within specific samples
    % (either seperate or combined throughout samples).
    %
    % Input:
    %   - app - Instance of CytoMAP
    
    alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

    [~, sample] = Helper.find_MFI(app);

    tmpData = Helper.populate_table(app, 'MFI', 'AllCells', 'smpls', {sample}, 'fill_checkbox', false);

    tData = cell(size(tmpData, 1) + 1, 7);
    % Put in phenotypes
    tData(1:size(tmpData, 1), 1)= tmpData(:, 2);
    tData(1:size(tmpData, 1), 2) = tmpData(:, 3);
    % Put in sample names
    tData(1:size(tmpData, 1), 3) = tmpData(:, 7);
    tData(1:size(tmpData, 1), 4) = tmpData(:, 8);
    
    % column for region groups
    smpls = fieldnames(app.data);
    groups = cell(1, numel(smpls));
    for smpli = 1:numel(smpls)
        if contains('MetaData', fieldnames(app.data.(smpls{smpli})))
            if contains('Group', fieldnames(app.data.(smpls{smpli}).MetaData))
                grp = app.data.(smpls{smpli}).MetaData.Group;
                if iscell(grp)
                    groups{smpli} = grp{1};
                else
                    groups{smpli} = grp;
                end
            else
                groups{smpli} = 0; 
            end
        else
            groups{smpli} = 0;
        end
    end
    tData(1:numel(smpls), 5) = groups;

    % column for include region
    
    % column for region numbers or names 
    nnms = fieldnames(app.net);
    if contains('RegNames', fieldnames(app.net.(nnms{1})))
        rnms = app.net.(nnms{1}).RegNames;
    else
        rnms = num2cell(1:app.net.(nnms{1}).NR);
    end
    tData(1:app.net.(nnms{1}).NR, 7) = rnms;
    tData(1:app.net.(nnms{1}).NR, 6) = {false};
    
    
    tData(end, :) = {false, 'Select All', false, 'Select All', [], [], []}; 

    % Build the main figure
    UIfig = uifigure('Name', 'Plot Region Cellular Composition Options', 'Scrollable', 'on');
    Helper.func_SetCLR(app, UIfig, 'UIfigure');
    UIfig.Position = alpha*[10 10 1100 800];

    % Select model type
    lbl = uilabel(UIfig); lbl.Text = sprintf('Select Model To Analyze:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = alpha*[10+200 45 400 15];
    DataType = uidropdown(UIfig);
    DataType.Items = fieldnames(app.net);
    DataType.Value = DataType.Items{1};
    DataType.Position = alpha*[10+200 15 175 30];
    DataType.BackgroundColor = app.GUIOPTS.bgclr;

    % Create the table of options
    t = uitable(UIfig);
    t.Data = tData;
    t.Position = alpha*[0 70 1100 730];
    t.ColumnName = {'Include','Phenotype (Must be in all samples)', ...
                    'Include', 'Sample', 'Sample Group', ...
                    'Analyze Cells From', 'Region'};
    t.ColumnEditable = [true false true false true true false];
    t.ColumnWidth = {75, 350, 75, 200, 100, 140, 100};
    t.CellEditCallback = @(dd, p) switched_sample(app, dd, p);

    % Select Color Scale
    lbl = uilabel(UIfig); lbl.Text = sprintf('Scale:');
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = alpha*[10+200+175 45 100 15];
    MAPScale = uidropdown(UIfig);
    MAPScale.Items = {'linear', ...
                      'log'};
    MAPScale.Value = {'linear'};
    MAPScale.Position = alpha*[10+200+175 15 100 30];
    MAPScale.BackgroundColor = app.GUIOPTS.bgclr;

    % Create a Norm per sample or per dataset label
    lbl = uilabel(UIfig);
    lbl.Text = sprintf('Type of Plot:');
    lbl.HorizontalAlignment = 'left';
    Helper.func_SetCLR(app, lbl, 'label')
    lbl.Position = alpha*[10 50 150 20];
    % Create a Norm per sample or per dataset button
    NormPer = uibuttongroup(UIfig, 'Visible','off');
    NormPer.Position = alpha*[10, 7, 200, 45];
    Helper.func_SetCLR(app, NormPer, 'UICpopup')
    % Create two radio buttons in the button group.
    r1 = uicheckbox(NormPer, 'Text','Number of Cells in Cluster',...
        'Position',alpha*[5, 25, 200, 20], 'Value', true);
    Helper.func_SetCLR(app, r1, 'table')
    r2 = uicheckbox(NormPer, 'Text','Percentage of total cells in Cluster',...
        'Position',alpha*[5, 5, 200, 20], 'Value', true);
    Helper.func_SetCLR(app, r2, 'table')
    NormPer.Visible = 'on';

    % Create a push button
    btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap(t, app, DataType, MAPScale, NormPer));
    btn.Position = alpha*[1000-150, 10, 100, 50];
    btn.Text = 'Ok';
    Helper.func_SetCLR(app, btn, 'button')
    
    DataType.ValueChangedFcn = @(btn, event) ChangeDataType(app, t, DataType.Value);

    function ChangeDataType(app, dd, netnm)
        % change available regions to edit
    
        % column for region numbers or names 
        if contains('RegNames', fieldnames(app.net.(netnm)))
            regnms = app.net.(netnm).RegNames;
        else
            regnms = num2cell(1:app.net.(netnm).NR);
        end
        dd.Data(:, 6) = {[]};
        dd.Data(:, 7) = {[]};
        dd.Data(1:app.net.(netnm).NR, 7) = regnms;
        dd.Data(1:app.net.(netnm).NR, 6) = {false};
        
    end

    
    function switched_sample(app, dd, p)
        if p.Indices(2) == 3
            % Make sure that at least one thing is selected
            if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 3)), 3)))
                dd.Data{p.Indices(1), 3} = true;
                return;
            end
            % Process Select All button
            fill_select_all = false;
            if p.Indices(1) == find(strcmp(dd.Data(:, 4), 'Select All'))
                ind = ~cellfun('isempty', dd.Data(:, 3));
                dd.Data(ind, 3) = {logical(p.NewData)};
                if ~logical(p.NewData)
                    dd.Data{strcmp(dd.Data(:, 4), app.DataN.Value), 3} = true;
                end
                fill_select_all = logical(p.NewData);
            end
            ind = dd.Data(:, 3);
            ind = ind(~cellfun('isempty', ind));
            ind = cell2mat(ind);

            if isempty(dd.Data(ind, 4))
                dd.Data(p.Indices(1), p.Indices(2)) = {true};
                return;
            end

            tmpDat = cell(size(dd.Data, 1) - 1, 7);

            tmpDat(1:end, 1) = {1}; % Stub. Needed for populate table
            tmpDat(1:end, 2) = dd.Data(1:end-1, 1);
            tmpDat(1:end, 3) = dd.Data(1:end-1, 2);
            tmpDat(1:end, 4) = {1}; % Stub. Needed for populate table
            tmpDat(1:end, 7) = dd.Data(1:end-1, 3);
            tmpDat(1:end, 8) = dd.Data(1:end-1, 4);

            tmpDat = Helper.populate_table(app, ...
                'smpls', dd.Data(ind, 4), ...
                'MFI', 'AllCells', ...
                'prev_table', tmpDat ...
            );

            % Only re-make the data table if you have to
            INDddCL = ~cellfun('isempty', dd.Data(:, 2));
            INDddCL(end) = false;
            INDtmpCL = ~cellfun('isempty', tmpDat(:, 3));

            if ~Helper.setequal(dd.Data(INDddCL, 2), tmpDat(INDtmpCL, 3))
                dd.Data(1:numel(tmpDat(:, 2)), 1) = tmpDat(:, 2);
                dd.Data(1:numel(tmpDat(:, 3)), 2) = tmpDat(:, 3);
                dd.Data(1:numel(tmpDat(:, 7)), 3) = tmpDat(:, 7);
                dd.Data(1:numel(tmpDat(:, 8)), 4) = tmpDat(:, 8);
                dd.Data(end, :) = {false, 'Select All', fill_select_all, 'Select All', [], [], []}; 
            end
        elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
            ind = ~cellfun('isempty', dd.Data(:, 2));
            dd.Data(ind, 1) = {logical(p.NewData)};
        end
    end

    function backend_wrap(t, app, DataType, MAPScale, NormPer)
        tmp_data = t.Data(1:end-1, :);
        plot_backend( ...
            app, ...
            tmp_data([tmp_data{:, 3}]'==1, 4), ... Samples
            tmp_data([tmp_data{:, 1}]'==1, 2), ... Phenotypes
            DataType.Value, ...
            MAPScale.Value, ...
            cell2mat(tmp_data([tmp_data{:, 3}]'==1, 5)), ... % Group Number
            find([tmp_data{:, 6}]'==1), ...        % region numbers
            tmp_data([tmp_data{:, 6}]'==1, 7), ...  % Region Names
            [NormPer.Children.Value] ...
        );
    end

    function plot_backend(app, smplnms, phenos, DataType, MAPScale, grpnum, regnums, regnms, plottypes)
            % Back-End of making plots function. Allows user to create
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
            if isempty(phenos)
                return
            end
            
            % re-order things if they are groups are wierdly ordered
            [~, INDgrp]= sort(grpnum);
            grpnum = grpnum(INDgrp);
            smplnms = smplnms(INDgrp);
            regnms(1:end-1) = strcat(regnms(1:end-1), ' and_');
            
            for smplj = 1:numel(smplnms)
                if contains('MetaData', fieldnames(app.data.(smplnms{smplj})))
                    if iscell(grpnum(smplj))
                        app.data.(smplnms{smplj}).MetaData.Group = grpnum(smplj);
                    else
                        app.data.(smplnms{smplj}).MetaData.Group = {grpnum(smplj)};
                    end
                else
                    app.data.(smplnms{smplj}).MetaData = struct;
                    if iscell(grpnum(smplj))
                        app.data.(smplnms{smplj}).MetaData.Group = grpnum(smplj);
                    else
                        app.data.(smplnms{smplj}).MetaData.Group = {grpnum(smplj)};
                    end
                end
            end
        
            
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
            %% Load the data                
            rmvzeros = 1;

            [Dat_PreMAIN, ~, INDDatPre, INDzrs, INDons] = Helper.func_loaddata( ...
                app, ...
                smplnms, ...
                phenos, ...
                [Constants.other_tag, DataType], ...
                [], ...
                [], ...
                'Individual Cells', ...
                'Cellularity: Number of Cells / Neighborhood', ...
                'Sample', ...
                rmvzeros, ...
                'Raw MFI per cell');

                % Pull the logicals
                Dat_PreLOGICMain = Dat_PreMAIN(:, 1:numel(phenos));
                % remove the binary cell type classifiers from the data table
                Dat_PreMAIN = Dat_PreMAIN(:, (numel(phenos)+1):end);
            

            %% Analyze the data
            datplt = zeros(numel(smplnms), 2*numel(phenos));
            Dat_PlotXBox = zeros(numel(smplnms), numel(grpnum)*numel(phenos));
            Dat_PlotYBox = zeros(numel(smplnms), numel(grpnum)*numel(phenos));
            
            for i=1:numel(smplnms)
                
                % Make a waitbar
                vPD = waitbar(0, 'Plotting heatmap', ...
                    'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));

                % Pull the ith data
                Dat_Pre = Dat_PreMAIN(INDDatPre{i}{1}(1):INDDatPre{i}{1}(2),:);
                Dat_PreLOGIC = Dat_PreLOGICMain(INDDatPre{i}{1}(1):INDDatPre{i}{1}(2),:);
                
                %% Build the table
                % Total Number of cells
                datplt(i, 1:numel(phenos)) = sum(Dat_PreLOGIC);
                if numel(regnums)>1
                    INDreg = any(Dat_Pre==regnums', 2);
                else
                    INDreg = Dat_Pre==regnums;
                end
                datplt(i, numel(phenos)+1:end) = sum(Dat_PreLOGIC(INDreg, :));
                
                % Do the plot if its the last sample or you are making
                if i == numel(smplnms)
                    for phenoi = 1:numel(phenos)
                        %% Plot the Number of cells
                        if plottypes(2)==1
                            
                            fig = figure;
                            fig.Color = app.GUIOPTS.bgclr;
                            fig.InvertHardcopy = 'off';

                            clf
                            % add an export data option
                            % Create Menu bar
                            ExportMenu = uimenu(fig);
                            ExportMenu.Text = 'Export';

                            % Create an export data option
                            ExportPltDat = uimenu(ExportMenu);
                            ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                            ExportPltDat.Text = 'Export Plot Data to .csv';
                            hold on
                            
                            for smplj = 1:numel(smplnms)
                                plot( ...
                                    grpnum(smplj), ... % X data
                                    datplt(smplj, numel(phenos)+phenoi), ...  % Y data
                                    'o', 'DisplayName', smplnms{smplj} ...
                                    )
                            end

                            % Label Stuff
                            ax = gca;
                            ax.YLabel.String = 'Number of Cells';
                            ax.XLim = [min(grpnum)-1, max(grpnum+1)];

                            % add the box and wiskers
                            boxplot( ...
                                    datplt(:, numel(phenos)+phenoi), ...
                                    grpnum, ...
                                    'Colors', [0.6, 0.6, 0.6], 'Symbol', '.k', 'Widths', 0.2 ...
                                    );

                            box off
                            ax = gca;
                            ax.XTick = min(grpnum):max(grpnum);
                            ax.XTickLabel = unique(grpnum);
                            ax.XLabel.String = 'Sample Group Number';
                            ax.Title.String = [phenos{phenoi} ' in regions ' regnms{:}];
                            ax.YLim(1) = 0;
                            ax.YScale = MAPScale;

                        end
                        %% Plot the Percentage of cells
                        if plottypes(1)==1
                            fig = figure;
                            fig.Color = app.GUIOPTS.bgclr;
                            fig.InvertHardcopy = 'off';

                            clf
                            % add an export data option
                            % Create Menu bar
                            ExportMenu = uimenu(fig);
                            ExportMenu.Text = 'Export';

                            % Create an export data option
                            ExportPltDat = uimenu(ExportMenu);
                            ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                            ExportPltDat.Text = 'Export Plot Data to .csv';
                            hold on

                            
                            for smplj = 1:numel(smplnms)
                                % Create the plots
                                numinreg = datplt(smplj, numel(phenos)+phenoi);
                                numtotal = datplt(smplj, phenoi);
                                plot( ...
                                    grpnum(smplj), ... % X data
                                    (numinreg./numtotal)*100, ...                % Y data
                                    'o', 'DisplayName', smplnms{i} ...
                                    )
                            end

                            % Label Stuff
                            ax = gca;
                            ax.YLabel.String = 'Percentage of Cells';
                            ax.XLim = [min(grpnum)-1, max(grpnum+1)];
                            % add the box and wiskers
                            boxplot( ...
                                    (datplt(:, numel(phenos)+phenoi)./datplt(:, phenoi))*100, ...
                                    grpnum, ...
                                    'Colors', [0.6, 0.6, 0.6], 'Symbol', '.k', 'Widths', 0.2 ...
                                    );

                            box off
                            ax = gca;
                            ax.XTick = min(grpnum):max(grpnum);

                            ax.XTickLabel = unique(grpnum);
                            ax.XLabel.String = 'Sample Group Number';
                            ax.Title.String = [phenos{phenoi} ' in regions ' regnms{:}];
                            ax.YLim(1) = 0;
                            ax.YScale = MAPScale;
                        end
                        
                    end
                end % End if plot

                waitbar(1,vPD, 'Done!');
                close(vPD)
            end

            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end
            

        end
   
    
end