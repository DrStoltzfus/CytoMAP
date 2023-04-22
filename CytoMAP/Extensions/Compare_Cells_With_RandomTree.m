function Compare_Cells_With_RandomTree(app)
    %% Build a user interface
    % Input:
    %   - app - Instance of CytoMAP
    %
    
    if ~Helper.any_sample(app) || ~Helper.any_net(app)
        return;
    end
    
    alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
    
    % Build the clustering options menus
    UIfig = uifigure('Name', 'Random Tree', 'Scrollable', 'on');
    Helper.func_SetCLR(app, UIfig, 'UIfigure');
    UIfig.Position = alpha*[10 10 1225 700];

    % Initialize the table of options
    t = uitable(UIfig);
    t.Position = alpha*[0 70 1215 630];
    
    %Build the table
    [~, sample] = Helper.find_MFI(app);
    tmpData = Helper.populate_table(app, 'smpls', {sample}, 'MFI', 'AllCells', 'fill_checkbox', false);

    tData = cell(size(tmpData, 1) + 1, 7);
    % Cells
    tData(1:sum(~cellfun(@isempty,tmpData(:, 2))), 1) = {1};
    tData(1:numel(tmpData(:, 2)), 2)= tmpData(:, 2);
    tData(1:numel(tmpData(:, 3)), 3) = tmpData(:, 3);
    tData(end, 2:3) = {false, 'Select All'};
    % MFI
    tData(1:numel(tmpData(:, 5)), 4) = tmpData(:, 5);
    tData(1:numel(tmpData(:, 6)), 5) = tmpData(:, 6);
    tData(end, 4:5) = {false, 'Select All'};
    % Sample
    tData(1:numel(tmpData(:, 7)), 6) = tmpData(:, 7);
    tData(1:numel(tmpData(:, 8)), 7) = tmpData(:, 8);
    tData(end, 6:7) = {false, 'Select All'};
    
    % Populate the table of options
    t.Data = tData;
    t.ColumnName = { ...
        'Group', 'Include', 'Phenotype', ...
        'Include', 'Channel MFI', ...
        'Include', 'Sample', ...
    };
    t.ColumnEditable = [true true false true false true false];
    t.ColumnWidth = {alpha*50, alpha*75, alpha*350, alpha*75, alpha*350, alpha*75, alpha*200};
    t.CellEditCallback = @(dd, p) switched_sample(app, dd, p);

    % Create a push button
    btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap( ...
        t, app ...
    ));
    btn.Position = alpha*[1215-100, 5, 100, 50];
    btn.Text = 'Ok';
    Helper.func_SetCLR(app, btn, 'button')
    
    function switched_sample(app, dd, p)
        % if the samples index was changed
        if p.Indices(2) == 6
            % Make sure that at least one thing is selected

            % Process Select All Button
            fill_select_all = false;
            if p.Indices(1) == find(strcmp(dd.Data(:, 7), 'Select All'))
                ind = ~cellfun('isempty', dd.Data(:, 6));
                dd.Data(ind, 6) = {logical(p.NewData)};
                if ~logical(p.NewData)
                    dd.Data{strcmp(dd.Data(:, 7), app.DataN.Value), 6} = true;
                end
                fill_select_all = logical(p.NewData);
            end

            ind = ~cellfun('isempty', dd.Data(:, 6));
            ind(end) = false;
            ind(ind) = logical(cell2mat(dd.Data(ind, 6)));

            smpls = dd.Data(ind, 7);
            if isempty(smpls)
                [~, smpls] = Helper.find_MFI(app);
                smpls = {smpls};
            end
            
            tmpDat = cell(size(dd.Data, 1) - 1, 8);

            tmpDat(1:end, 1) = {1}; % Stub. Needed for populate table
            tmpDat(1:end, 2) = dd.Data(1:end-1, 2);
            tmpDat(1:end, 3) = dd.Data(1:end-1, 3);
            tmpDat(1:end, 4) = {1}; % Stub. Needed for populate table
            tmpDat(1:end, 5) = dd.Data(1:end-1, 4);
            tmpDat(1:end, 6) = dd.Data(1:end-1, 5);

            tmpDat = Helper.populate_table(app, ...
                'smpls', smpls, ...
                'MFI', 'AllCells', ...
                'prev_table', tmpDat);

            dd_ph = ~cellfun('isempty', dd.Data(:, 2));
            dd_ph(end) = false;
            dd_ph = dd.Data(dd_ph, 3);
            dd_mfi = ~cellfun('isempty', dd.Data(:, 4));
            dd_mfi(end) = false;
            dd_mfi = dd.Data(dd_mfi, 5);

            tmp_ph = ~cellfun('isempty', tmpDat(:, 2));
            tmp_ph = tmpDat(tmp_ph, 3);
            tmp_mfi = ~cellfun('isempty', tmpDat(:, 5));
            tmp_mfi = tmpDat(tmp_mfi, 6);

            if ~Helper.setequal(dd_ph, tmp_ph) || ~Helper.setequal(dd_mfi, tmp_mfi)
                dd.Data = cell(size(tmpDat, 1) + 1, 7);
                dd.Data(1:numel(tmpDat(:, 2)), 1) = {1};
                dd.Data(1:numel(tmpDat(:, 2)), 2) = tmpDat(:, 2);
                dd.Data(1:numel(tmpDat(:, 3)), 3) = tmpDat(:, 3);
                dd.Data(1:numel(tmpDat(:, 5)), 4) = tmpDat(:, 5);
                dd.Data(1:numel(tmpDat(:, 6)), 5) = tmpDat(:, 6);
                dd.Data(1:numel(tmpDat(:, 7)), 6) = tmpDat(:, 7);
                dd.Data(1:numel(tmpDat(:, 8)), 7) = tmpDat(:, 8);
                
                dd.Data(end, 2:end) = {false, 'Select All', false, 'Select All', fill_select_all, 'Select All'};
            end                            
        elseif p.Indices(1) == find(strcmp(dd.Data(:, 3), 'Select All')) && p.Indices(2) == 2
            ind = ~cellfun('isempty', dd.Data(:, p.Indices(2)));
            dd.Data(ind, p.Indices(2)) = {logical(p.NewData)};
        elseif p.Indices(1) == find(strcmp(dd.Data(:, 5), 'Select All')) && p.Indices(2) == 4
            ind = ~cellfun('isempty', dd.Data(:, p.Indices(2)));
            dd.Data(ind, p.Indices(2)) = {logical(p.NewData)};
        end
    end

    function backend_wrap(t, app)
        tmp_data = t.Data(1:end-1, :); %%remove the select all 
        
        INDSmpls = [tmp_data{1:numel(app.DataN.Items), 6}];

        smplnms = tmp_data(INDSmpls, 7);

        phenos = tmp_data(~cellfun('isempty', tmp_data(:, 3)),3)';
        phenos = phenos([tmp_data{:,2}]==1);
        
        Groups = [tmp_data{[tmp_data{:,2}]==1, 1}];
        if numel(unique(Groups))==1
            Groups = 1:numel(Groups);
        end
        
        MFIList = tmp_data(~cellfun('isempty', tmp_data(:, 5)), 5);
        MFIList = MFIList([tmp_data{:, 4}]==1);
        MFIList = Helper.valid_channel(MFIList);

        random_tree_backend( ...
            app, ...
            Groups, ...
            smplnms, ... Sample
            phenos, ... Phenotypes
            MFIList ...
        );
    end

    function random_tree_backend(app, Groups, smplnms, phenos, MFIList)
        
        %% Use random forest tree to determine predictor importance
        % RANDOM_TREE_BACKEDN Back_End of random tree. Allow user to
        % visulize 
        %
        % Input:
        %   - app - Instance of CytoMAP
        %   - smplnms - cell - Names of samples to plot
        %   - Groups - array - Same size as smplnms. Groupings of samples.
        %       Allows to treat multiple samples as one, without needing
        %       to merge them etc.
        %   - phenos - 
        %   - MFIList - 
        
        % order the samples in the order of groups
        
        [Dat_Pre, ~, INDDatPre, INDzrs, INDons] = Helper.func_loaddata(app, ...
            smplnms, ...
            phenos, ...
            MFIList, ...
            [], ...
            [], ...
            'Individual Cells', ...
            'Cellularity: Number of Cells / Neighborhood', ...
            'Sample', ...
            1); 
        
        Group_id = ones(size(Dat_Pre,1),1);
        
        for idx = 1: length(Groups)
            Group_id(Dat_Pre(:,idx)==1) = Groups(idx);
        end
        
        Dat_Channel = Dat_Pre(:, (numel(phenos)+1):end);
        Dat_Channel = mat2dataset(Dat_Channel);
        Dat_Channel.Properties.VarNames = MFIList;
        Dat_Channel.Group_id = Group_id;
        rng(1); % For reproducibility
        MdlDefault = fitctree(Dat_Channel,'Group_id','CrossVal','on');
        view(MdlDefault.Trained{1},'Mode','graph');
        Mdl = fitctree(Dat_Channel,'Group_id','PredictorSelection', 'curvature');
        imp = predictorImportance(Mdl);
        
        fig = figure;
        
        fig.Color = app.GUIOPTS.bgclr;
        fig.InvertHardcopy = 'off';
        
        barh(imp);
        title(sprintf('Predictor Importance Estimates\n (estimated using fitctree)'));
        xlabel('Estimates');
        ylabel('Predictors');
        h = gca;
        h.YTick = 1:numel(Mdl.PredictorNames);
        h.YTickLabel = Mdl.PredictorNames;
        h.TickLabelInterpreter = 'none';
    end
end