classdef Plt_Helper

% Plt_Helper defines a suite of functions which help plot and visualize data

    methods (Static)
        function Import_Definitions_Func
% % %             import Helper.*;
% % %             import Constants.*;
            set(groot, 'DefaultTextInterpreter', 'none')
            set(groot, 'DefaultLegendInterpreter', 'none')
        end

        function [dat, vrnms, pnms] = UpdateFigMenu(app, num, varargin)
            % UPDATEFIGMENU Returns all of the points which should be
            % currently plotted in the given figure based on chices on
            % samples and phenotypes. It also provides an easy access to
            % variable names on all of the axes (up to special cases of
            % names).
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - num - string, char - Identifier for a figure to pull data
            %       from
            %
            % Key-word Input:
            %   - Combine - default: true - Whether to combine all of data
            %       into one table, or to make a table for each of the
            %       sample/phenotype combinations currently chosen.
            %
            % Output:
            %   - dat - table, cell - table if combine was true, cell if
            %       combine was false. It has all of the data currently
            %       plotted with all rows being unique (if Combine was
            %       true). If it's cell, then it has size of:
            %       (#samples, #phenotypes), with each entry being a table
            %       corresponding to data from that sample filtered through
            %       that phenotype.
            %   - vrnms - cell - Cell of strings/chars which should be
            %       contents of x-axis. They can be also used on y and z
            %       axes, since everything additional on those 2 axes are
            %       in the end of options. However, IT CANNOT be used on
            %       the color axis, since it can lead to index mismatch.
            %   - pnms - cell - Cell of string/chars of all phenotypes
            %       which are currently chosen in a given figure.
            
            % Define defualt arguments for varargin
            defaults = struct( ...
                ... % For whether to combine data from all phenotypes and samples in one table,
                ... % or return it as a cell of tables with shape (num_smpl, num_phenos)
                'Combine', true ...
            );
            % Check that varargin is of even length
            if mod(length(varargin), 2) ~= 0
                error('Uneven number of input arguments! func_newfig needs propertyName/propertyValue pairs, after app argument')
            end
            % Process key, i.e. change default values.
            for pair = reshape(varargin, 2, [])
                if isfield(defaults, pair{1})
                    defaults.(pair{1}) = pair{2}; % Change default axis value to the given value
                else  % Key is not in the defaults
                    error('%s is not a recognized parameter name', pair{1})
                end
            end

            % Update Sample drop down in case you uploaded more files.
            app.figOpts.sample.(num).String = fieldnames(app.data);
            % Pull the selected sample name
            smpls = app.figOpts.smpls.(num);  % app.DataN.Items{app.figOpts.sample.(num).Value};

            % Returns the number of gates
            pnms = Helper.get_gate_tags(app, smpls)';
            % If the there is an AllCells "0" level gate, remove it (this is the same as AllCells for plotting)
            pnms = pnms(~strcmp(pnms, [Constants.gate_tag '0']));

            phnm = app.figOpts.phenos.(num);
            if ~iscell(phnm)  % Sanity check, if someone decides to change field to not cell
                app.figOpts.phenos.(num) = {app.figOpts.phenos.(num)};
                phnm = app.figOpts.phenos.(num);
            end
            
            %% Pull Data to be plotted
            if defaults.Combine
                %% If the user wants to combine all of the data
                if numel(phnm) == 1 && strcmp(phnm, 'Density/MFI RSN')
                    dat = Helper.merge_MFI(app, app.figOpts.smpls.(num), 'MFIRSN');
                elseif numel(phnm) == 1 && strcmp(phnm, 'Density/MFI CCN')
                    dat = Helper.merge_MFI(app, app.figOpts.smpls.(num), 'MFICCN');
                elseif numel(phnm) == 1 && strcmp(phnm, 'All Cells')
                    [ch_oth_table, gate_table, ~] = Helper.merge_AllCells(app, app.figOpts.smpls.(num));
                    allcells = horzcat(ch_oth_table, gate_table);
                    [dat, ~] = unique(allcells, 'rows', 'stable');
                elseif all(startsWith(phnm, Constants.neigh_tag))
                    dat = Helper.merge_MFI(app, app.figOpts.smpls.(num), 'MFIRSN');
                    if isempty(dat) || ~ismember(phnm(1), dat.Properties.VariableNames)
                        dat = Helper.merge_MFI(app, app.figOpts.smpls.(num), 'MFICCN');
                        if isempty(dat) || ~ismember(phnm(1), dat.Properties.VariableNames)
                            errordlg('Error occured when loading neighborhood');
                            error('Neither MFIRSN or MFICCN is in all samples, or first neighborhood is not in all of the samples.');
                        end
                    end
                    logic = zeros(size(dat, 1), 1);
                    for ph_idx=1:numel(phnm)
                        logic = logic | table2array(dat(:, app.figOpts.phenos.(num)(ph_idx)));
                        if isempty(dat) || ~ismember(phnm(1), dat.Properties.VariableNames)
                            errordlg('Error occured when loading neighborhood');
                            error('%d neighborhood is not in all of the samples.', ph_idx);
                        end
                    end
                    dat = dat(logic, :);
                else
                    % Pull the data that is to be plotted
                    [ch_oth_table, gate_table, gate_tags] = Helper.merge_AllCells(app, app.figOpts.smpls.(num)); % SLOW
                    allcells = horzcat(ch_oth_table, gate_table);
                    [allcells, ~] = unique(allcells, 'rows', 'stable');
                    logic = zeros(size(allcells, 1), 1);
                    % If any of the phenotypes aren't in AllCells
                    if any(~ismember(phnm, allcells.Properties.VariableNames))
                        for ph_idx=1:numel(phnm)
                            tag_idx = find(strcmp(gate_tags{2, :}, phnm(ph_idx)));
                            if ~isempty(tag_idx)
                                phnm(ph_idx) = gate_tags.Properties.VariableNames(tag_idx);
                            end
                        end
                    end
                    for ph_idx=1:numel(phnm)
                        if ismember(app.figOpts.phenos.(num)(ph_idx), {'Density/MFI RSN', 'Density/MFI CCN', 'All Cells'}) || ...
                            startsWith(app.figOpts.phenos.(num)(ph_idx), Constants.neigh_tag)
                            continue;
                        end
                        logic = logic | table2array(allcells(:, phnm(ph_idx)));
                    end
                    dat = allcells(logic, :);
                end
                %% Pull out the available variable names
                dat = Helper.reorder_cols(dat);
                % Move any Model Names to the end
                modelnames = strcat(Constants.other_tag, fieldnames(app.net)); 
                for m_i = 1:numel(modelnames)
                    if ismember(modelnames{m_i}, dat.Properties.VariableNames)
                        dat = movevars(dat, modelnames{m_i}, 'After', size(dat, 2));
                    end
                end
                vrnms = dat.Properties.VariableNames;
                gate_neigh_idxs = startsWith(vrnms, Constants.gate_tag) | startsWith(vrnms, Constants.neigh_tag);
                tmp = vrnms(~gate_neigh_idxs);
                tmp(end + 1: end + sum(gate_neigh_idxs)) = vrnms(gate_neigh_idxs);
                vrnms = tmp;
            else % If combine is not true
                dat = cell(numel(smpls), numel(phnm));
                if numel(phnm) == 1 && strcmp(phnm, 'Density/MFI RSN')
                    for smpl_idx = 1:numel(smpls)
                        smpl = smpls{smpl_idx};
                        dat(smpl_idx, 1) = {app.data.(smpl).MFIRSN};
                    end
                elseif numel(phnm) == 1 && strcmp(phnm, 'Density/MFI CCN')
                    for smpl_idx = 1:numel(smpls)
                        smpl = smpls{smpl_idx};
                        dat(smpl_idx, 1) = {app.data.(smpl).MFICCN};
                    end
                elseif numel(phnm) == 1 && strcmp(phnm, 'All Cells')
                    for smpl_idx = 1:numel(smpls)
                        smpl = smpls{smpl_idx};
                        dat(smpl_idx, 1) = {app.data.(smpl).AllCells};
                    end
                elseif all(startsWith(phnm, Constants.neigh_tag))
                    for phn_idx = 1:numel(phnm)
                        phn = phnm{phn_idx};
                        in_rsn = true;
                        % Determine if it's MFI RSN or CCN. Consistently for all Samples
                        for smpl_idx = 1:numel(smpls)
                            smpl = smpls{smpl_idx};
                            if ~isfield(app.data.(smpl), 'MFIRSN') || ...
                                    ~ismember({phn}, app.data.(smpl).MFIRSN.Properties.VariableNames)
                                in_rsn = false;
                            end
                        end
                        % Load all the neighborhoods into dat.
                        for smpl_idx = 1:numel(smpls)
                            smpl = smpls{smpl_idx};
                            if in_rsn
                                dat{smpl_idx, phn_idx} = app.data.(smpl).MFIRSN(app.data.(smpl).MFIRSN.(phn)==1, :);
                            else
                                dat{smpl_idx, phn_idx} = app.data.(smpl).MFICCN(app.data.(smpl).MFICCN.(phn)==1, :);
                            end
                        end
                    end
                else
                    % Pull the data that is to be plotted
                    for smpl_idx = 1:numel(smpls)
                        smpl = smpls{smpl_idx};
                        for phn_idx = 1:numel(phnm)
                            phn = phnm{phn_idx};
                            if ismember({phn}, {'Density/MFI RSN', 'Density/MFI CCN', 'All Cells'}) || ...
                                startsWith({phn}, Constants.neigh_tag)
                                % Just a sanity check it shouldn't be possible to come in here.
                                warndlg("Something went wrong. Neighborhoods and phenotypes shouldn't be mixed. You might need to restart CytoMap.");
                                continue;
                            end
                            tag = Helper.gate_full2tag(app, phn, smpl);
                            if iscell(tag)
                                tag = tag{1};
                            end
                            dat{smpl_idx, phn_idx} = app.data.(smpl).AllCells(app.data.(smpl).AllCells.(tag)==1, :);
                        end
                    end
                end
                dat{1, 1} = Helper.reorder_cols(dat{1, 1});
                % Move any Model Names to the end
                modelnames = strcat(Constants.other_tag, fieldnames(app.net)); 
                for m_i = 1:numel(modelnames)
                    if ismember(modelnames{m_i}, dat{1, 1}.Properties.VariableNames)
                        dat{1, 1} = movevars(dat{1, 1}, modelnames{m_i}, 'After', size(dat{1, 1}, 2));
                    end
                end
                vrnms = dat{1, 1}.Properties.VariableNames;
                gate_neigh_idxs = startsWith(vrnms, Constants.gate_tag) | startsWith(vrnms, Constants.neigh_tag);
                tmp = vrnms(~gate_neigh_idxs);
                tmp(end + 1: end + sum(gate_neigh_idxs)) = vrnms(gate_neigh_idxs);
                vrnms = tmp;
                for v_i= 1:numel(vrnms)
                    if isstring(vrnms{v_i})
                        vrnms{v_i} = char(vrnms{v_i});
                    end
                end
                for d_i = 2:numel(dat)
                    dat{d_i} = Helper.reorder_cols(dat{d_i});
                    % Move any Model Names to the end
                    modelnames = strcat(Constants.other_tag, fieldnames(app.net)); 
                    for m_i = 1:numel(modelnames)
                        if ismember(modelnames{m_i}, dat{d_i}.Properties.VariableNames)
                            dat{d_i} = movevars(dat{d_i}, modelnames{m_i}, 'After', size(dat{d_i}, 2));
                        end
                    end
                
                    vrnms_i = dat{d_i}.Properties.VariableNames;
                    for v_i= 1:numel(vrnms_i)
                        if isstring(vrnms_i{v_i})
                            vrnms_i{v_i} = char(vrnms_i{v_i});
                        end
                    end
                    vrnms = intersect(vrnms, vrnms_i, 'stable');
                end
                vrnms = vrnms(~ismember(vrnms, {'X', 'Y', 'Z'}));
                vrnms = [{'X', 'Y', 'Z'}, vrnms];
            end % end if combine statement
            %% Now make the variable names menus
            vrnmsREP = vrnms;

            % Gate tag to Cell name could be specific per sample
            % dat will always have tags from 1st sample
            vrnmsREP(startsWith(vrnmsREP, Constants.gate_tag)) = Helper.gate_tag2full(app, vrnmsREP(startsWith(vrnmsREP, Constants.gate_tag)), smpls{1});
            % remove the channel tag
            vrnmsREP(startsWith(vrnmsREP, Constants.channel_tag)) = Helper.full_channel(vrnmsREP(startsWith(vrnmsREP, Constants.channel_tag)));
            % remove the 'other' tag
            vrnmsREP(startsWith(vrnmsREP, Constants.other_tag)) = Helper.full_other(vrnmsREP(startsWith(vrnmsREP, Constants.other_tag)));
            % remove neighborhood tag
            vrnmsREP(startsWith(vrnmsREP, Constants.neigh_tag)) = Helper.full_neigh(vrnmsREP(startsWith(vrnmsREP, Constants.neigh_tag)));
            % Update drop down menus
            app.figOpts.xaxIF.(num).String = vrnmsREP;
            % These are the same for all cases
            app.figOpts.zaxIF.(num).String = app.figOpts.xaxIF.(num).String;
            app.figOpts.yaxIF.(num).String = [app.figOpts.xaxIF.(num).String', {'HistogramNum', 'HistogramFreq'}];
            if app.figOpts.D3.(num)
                if any(strcmp({'Point Density Grid'}, app.figOpts.caxIF.(num).String)) && app.figOpts.caxIF.(num).Value > 3
                    app.figOpts.caxIF.(num).Value = app.figOpts.caxIF.(num).Value - 2;
                end
                app.figOpts.caxIF.(num).String = [{'None'}, {'Point Density'}, app.figOpts.xaxIF.(num).String'];
            else
                if ~any(strcmp({'Point Density Grid'}, app.figOpts.caxIF.(num).String)) && app.figOpts.caxIF.(num).Value > 3
                    app.figOpts.caxIF.(num).Value = app.figOpts.caxIF.(num).Value + 1;
                end
                app.figOpts.caxIF.(num).String = [{'None'}, {'Number of cells / Neighborhood'}, {'Point Density Grid'}, {'Point Density Contour'}, {'Point Density'}, app.figOpts.xaxIF.(num).String'];
            end
            app.figOpts.sldrIF.(num).String = app.figOpts.xaxIF.(num).String;
        end

        function func_poly(app, name0, web)
            % FUNC_POLY Creates a polygon on a current plot, and allows
            % user to modify it
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            %   - name0 - Name of the polygon
            % 
            % Modifies:
            %   - app - In CytoMap it modifies polygons field, either
            %       adding to it, or removing from it.
            if nargin<3
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            %% Build the options menus
            options = cell(2,4);

            if isempty(name0)
                %Draw polygons
                POLY = drawpolygon(gca, 'Color', 'r');
                app.polygons.TMP = struct;
                app.polygons.TMP.P1.obj = POLY;
                app.polygons.TMP.P1.pos = POLY.Position;
                app.polygons.TMP.P1.clr = POLY.Color;
               
                options(:,1) = {'r', 'New_polygon'};
            else
                options(:,1) = {'r', name0};
                polysubnames = fieldnames(app.polygons.(name0));
                for poly_i=1:numel(polysubnames)
                    POLY = app.polygons.(name0).(polysubnames{poly_i});
                    iPOLY = drawpolygon(gca, 'Position', POLY.pos, 'Color', POLY.clr);
                    iPOLY.InteractionsAllowed = 'all';
                    app.polygons.TMP.(['P' num2str(poly_i)]).obj = iPOLY;
                    app.polygons.TMP.(['P' num2str(poly_i)]).pos = iPOLY.Position;
                    app.polygons.TMP.(['P' num2str(poly_i)]).clr = iPOLY.Color;
                   
                end
            end
            options(:,2) = {'Color', 'Polygon Name'};

            UIfig = uifigure('Name', 'Polygon Options');
            if web==1
                UIfig.Visible='OFF';
            end

            UIfig.Position(3:4) = alpha*[350 250];
            % Create the table of options
            t = uitable(UIfig);
            t.Data = options;
            t.Position = alpha*[0 40 350 210];
            t.ColumnName = {'Value', 'Option', 'Statistics'};
            t.ColumnEditable = [true, false, false, false];
            t.ColumnWidth = {alpha*100, alpha*210, alpha*100, alpha*210};

            % Add polygon
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) func_poly_add(t,app));
            btn.Position = alpha*[10, 5, 100, 30];
            btn.Text = 'Add Polygon';
            Helper.func_SetCLR(app, btn, 'button')

            %save/exit
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) func_poly_save(t,app,UIfig));
            btn.Position = alpha*[340-75, 5, 75, 30];
            btn.Text = 'Save';
            Helper.func_SetCLR(app, btn, 'button')

            % delete polygon
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) func_poly_deleteLast(app));
            btn.Position = alpha*[110, 5, 75, 30];
            btn.Text = 'Del.';
            Helper.func_SetCLR(app, btn, 'button')

            % delete all polygons
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) func_poly_delete(app));
            btn.Position = alpha*[185, 5, 75, 30];
            btn.Text = 'Del. All';
            Helper.func_SetCLR(app, btn, 'button')

        % % %
        % % %     % Gate on cells within surface
        % % %     btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) ExportPopulation(t,app, num, smpl, Type));
        % % %     btn.Position = [10, 5, 115, 30];
        % % %     btn.Text = 'Export Population';


            function func_poly_add(~, app)
                % Draw polygons
                POLY2 = drawpolygon(gca, 'Color', 'r');
                Npoly = numel(fieldnames(app.polygons.TMP))+1;
                app.polygons.TMP.(['P' num2str(Npoly)]).obj = POLY2;
                app.polygons.TMP.(['P' num2str(Npoly)]).pos = POLY2.Position;
                app.polygons.TMP.(['P' num2str(Npoly)]).clr = POLY2.Color;
            end
            function func_poly_deleteLast(app)
                delete(app.polygons.TMP.(['P' num2str(numel(fieldnames(app.polygons.TMP)))]).obj)
                app.polygons.TMP = rmfield(app.polygons.TMP, ['P' num2str(numel(fieldnames(app.polygons.TMP)))]);
            end
            function func_poly_delete(app)
                pnms = fieldnames(app.polygons.TMP);
                for i=1:numel(fieldnames(app.polygons.TMP))
                    delete(app.polygons.TMP.(pnms{i}).obj)
                    app.polygons.TMP = rmfield(app.polygons.TMP, pnms{i});
                end
                app.polygons.TMP = struct;
            end
            function func_poly_save(t,app,UIfig)
                name = Helper.valid_var(t.Data{2,1});
                app.polygons.(name) = app.polygons.TMP;
                app.polygons = rmfield(app.polygons, 'TMP');
                close(UIfig)
            end
        end

        function update_plot_order(app, num)
            % UPDATE_PLOT_ORDER Rearranges currently plotted things in a
            % given figure in such a way that they correspond to the order
            % user set them in
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - num - char, string - Identifier of figure on which order
            %       of children will be changed.
            %
            % Modifies:
            %   - app - Specifically
            %       app.figOpts.order_stack/order_items/order_list
            %       By rearranging, truncating, adding new items to it,
            %       while trying to preserve maximal amount of information.
            
            ax_child = app.figOpts.plt_ax.(num).Children;
            combined_list = cell(numel(ax_child), 1);
            for child_idx = 1:numel(ax_child)
                combined_list{child_idx} = ax_child(child_idx).DisplayName;
            end
            old_out = ~ismember(app.figOpts.order_items.(num), combined_list);
            new_in  = ~ismember(combined_list, app.figOpts.order_items.(num));
            if any(old_out)
                app.figOpts.order_stack.(num) = app.figOpts.order_stack.(num)(~old_out);
                [~, idxs] = sort(app.figOpts.order_stack.(num));
                app.figOpts.order_stack.(num)(idxs) = 1:numel(app.figOpts.order_stack.(num));
                app.figOpts.order_items.(num) = app.figOpts.order_items.(num)(~old_out);
            end
            if any(new_in)
                app.figOpts.order_stack.(num)(end+1:end+sum(new_in)) = numel(app.figOpts.order_stack.(num)) + (1:sum(new_in));
                app.figOpts.order_items.(num)(end+1:end+sum(new_in)) = combined_list(new_in);
            end
            if isfield(app.figOpts.order_list, num) && isvalid(app.figOpts.order_list.(num))
                app.figOpts.order_list.(num).Items = app.figOpts.order_items.(num);
            end
        end
        
        function func_EditColormap(app, web)
            % FUNC_EDITCOLORMAP Lets the user change the color-map for the
            % selected clustering model,
            %
            % Input:
            %   - app - Instance of CytoMAP
            % 
            % Modifies:
            %   - Colors of the regions corresponding to the regions of a
            %       chosen model.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            if ~Helper.any_sample(app)
                return;
            end
            %% Make table options
            % make sure there are actually models to rename
            if isempty(fieldnames(app.net))
                return;
            end
            models = fieldnames(app.net);

            %just in case, make sure the model names is a cell
            if ~iscell(models)
                models = {models};
            end
            % Pull region colors
            colors =  app.net.(models{1}).cmap;

            % build region names
            regions = cell(size(colors,1),1);
            for region_i = 1:size(colors,1)
%                 regions{region_i} = ['R' num2str(region_i-1)];
                regions{region_i} = region_i-1;
            end

            tData = cell(size(colors,1),5);
            tData(1:end,1) = {false};
            tData(1:end,2) = regions;
            tData(1:end,3) = mat2cell(colors(:,1),ones(size(colors,1),1));
            tData(1:end,4) = mat2cell(colors(:,2),ones(size(colors,1),1));
            tData(1:end,5) = mat2cell(colors(:,3),ones(size(colors,1),1));
            %% make the table edit function
            function switched_sample(app, t, DataType)
                % Pull region colors
                colors_i =  app.net.(DataType.Value).cmap;
                % build region names
                regions_j = cell(size(colors_i,1),1);
                for r_i = 1:size(colors_i,1)
%                     regions_j{r_i} = ['R' num2str(r_i-1)];
                    regions_j{r_i} = r_i-1;
                end
                t_TMP = cell(size(colors_i,1),5);
                t_TMP(1:end,1) = {false};
                t_TMP(1:end,2) = regions_j;
                t_TMP(1:end,3) = mat2cell(colors_i(:,1),ones(size(colors_i,1),1));
                t_TMP(1:end,4) = mat2cell(colors_i(:,2),ones(size(colors_i,1),1));
                t_TMP(1:end,5) = mat2cell(colors_i(:,3),ones(size(colors_i,1),1));
                t.Data = t_TMP;
            end


            %% Make a table figure
            UIfig = uifigure('Name', 'Edit Regions and Colors');
            if web==1
                UIfig.Visible='OFF';
            end
            UIfig.Position = alpha*[10 10 550 800];
            t = uitable(UIfig);
            t.Data = tData;
            t.Position = alpha*[0 70 550 730];
            t.ColumnName = {'Edit','Region', ...
                            'Red', 'Green','Blue'};
            t.ColumnEditable = [true true true true true];
            t.ColumnWidth = {alpha*100, alpha*100, alpha*100, alpha*100, alpha*100};
            %t.ColumnWidth = cellfun(@(x) alpha*x,{100, 100, 100, 100, 100});

            % Select Neighborhood type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Clustering Model:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10 40 400 15];

            DataType = uidropdown(UIfig);
            DataType.Items = fieldnames(app.net);
            DataType.Value = DataType.Items{1};
            DataType.Position = alpha*[10 5 235 30];
            DataType.BackgroundColor = app.GUIOPTS.bgclr;
            DataType.ValueChangedFcn = @(~, ~) switched_sample(app, t, DataType);

            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) func_Edit(t));
            btn.Position = alpha*[10+235, 5, 100, 50];
            btn.Text = 'Edit Selected';
            Helper.func_SetCLR(app, btn, 'button')

            % Create a confirm colors push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Plt_Helper.edit_colormap_backend( ...
                app, DataType.Value, cell2mat(t.Data(1:end,2)), cell2mat(t.Data(1:end,3:5)) ...
            ));
            btn.Position = alpha*[10+235+100+100, 5, 100, 50];
            btn.Text = 'Save Changes';
            Helper.func_SetCLR(app, btn, 'button')
            
            % Create a reset colors push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) func_Reset(t));
            btn.Position = alpha*[10+235+100, 5, 100, 50];
            btn.Text = 'Reset Colors';
            Helper.func_SetCLR(app, btn, 'button')
            

            function func_Edit(t)
                IND = [t.Data{1:end,1}];
                if sum(IND) == 0
                    return;
                end
                selected = cell2mat(t.Data(IND,3:5));
                c = uisetcolor(selected(1,:));
                t.Data(IND,3) = {c(1)};
                t.Data(IND,4) = {c(2)};
                t.Data(IND,5) = {c(3)};
            end

            function func_Reset(t)
                % Generalize new colors to any number of regions
                new_color=jet(size(t.Data,1));
                % Make the first region white
                t.Data(1,3:5) = {1,1,1};
                I = 2:size(t.Data,1);
                % Add the data back to the table
                t.Data(I,3:5) = num2cell(new_color(I,:));
                %run apply changes
%                 Plt_Helper.edit_colormap_backend(app, DataType.Value, cell2mat(t.Data(1:end,3:5)));
            end
        end
        
        function edit_colormap_backend(app, DataType, newORDER, rgb)
            % ask the user to name the model
            prompt = {'Enter name of new model:'};
            title = 'Name new model';
            dims = [1 35];
            definput = {[DataType '_New' ]};
            ModelName = inputdlg(prompt,title,dims,definput);
            if isempty(ModelName)
                return;
            end
            % Make sure the user actually wants to overwrite an existing
            % model
            if any(ismember(fieldnames(app.net), ModelName))
                answer = questdlg(['Are you certain you want to overwrite' ModelName], ...
                    'Dessert Menu', ...
                    'Yes', 'No', 'No');
                if strcmp(answer, 'No') || isempty(answer)
                    % User decided not to override data, abort 
                    return;
                end
                % If the user chose to override an existing model do
                % that
                app.net.(ModelName{1}) = app.net.(DataType);
            else
                % If the model does not exist add it now
                app.net.(ModelName{1}) = app.net.(DataType);
            end
            
            % Edit the region colormap
            app.net.(ModelName{1}).cmap = rgb;
            % Pull all smaple names
            smpls = fieldnames(app.data);
            % Define the original numbering scheme
            oldORDER = 0:app.net.(DataType).NR;
            %% Check to see if the region numbers have been chaged
            if any(oldORDER~=newORDER)
                % Now that we are sure the user really wants to do this,
                % rearrange the numbers in all of the data
                for sample_i = 1:numel(smpls)
                    %% initialize some stuff
                    RowOrder = zeros(1,2);
                    ORDER = newORDER;
                    if ismember([Constants.other_tag DataType], ...
                            fieldnames(app.data.(smpls{sample_i}).AllCells))
                        ROWAC = app.data.(smpls{sample_i}).AllCells.([Constants.other_tag DataType]);
                        IND_AC = zeros(numel(ROWAC), numel(ORDER));
                    end
                    if ismember('MFIRSN', fieldnames(app.data.(smpls{sample_i})))
                        if ismember([Constants.other_tag DataType], ...
                                fieldnames(app.data.(smpls{sample_i}).MFIRSN))
                        ROWRSN = app.data.(smpls{sample_i}).MFIRSN.([Constants.other_tag DataType]);
                        IND_RSN = zeros(numel(ROWRSN), numel(ORDER));
                        end
                    end
                    if ismember('MFICCN', fieldnames(app.data.(smpls{sample_i})))
                        if ismember([Constants.other_tag DataType], ...
                                fieldnames(app.data.(smpls{sample_i}).MFICCN))
                        ROWCCN = app.data.(smpls{sample_i}).MFICCN.([Constants.other_tag DataType]);
                        IND_CCN = zeros(numel(ROWCCN), numel(ORDER));
                        end
                    end
                    %% Loop through the region numbers and re-define them
                    for i=1:numel(ORDER)
                        % If the model is in this AllCells in sample
                        if ismember([Constants.other_tag DataType], ...
                                fieldnames(app.data.(smpls{sample_i}).AllCells))
                            if ORDER(i) ~= i-1
                                RowOrder = [RowOrder; i-1, ORDER(i)];
                                % Pull the indeces of the elements
                                IND_AC(:, i) = ROWAC==(i-1);
                            end
                        end
                        % If the model is in MFIRSN
                        if ismember('MFIRSN', fieldnames(app.data.(smpls{sample_i})))
                            if ismember([Constants.other_tag DataType], ...
                                    fieldnames(app.data.(smpls{sample_i}).MFIRSN))
                                if ORDER(i) ~= i-1
                                    IND_RSN(:, i) = ROWRSN==(i-1);
                                end
                            end
                        end
                        % If the model is in MFICCN
                        if ismember('MFICCN', fieldnames(app.data.(smpls{sample_i})))
                            if ismember([Constants.other_tag DataType], ...
                                    fieldnames(app.data.(smpls{sample_i}).MFICCN))
                                if ORDER(i) ~= i-1
                                    IND_CCN(:, i) = ROWCCN==(i-1);
                                end
                            end
                        end
                    end % end loop through region numbers
                    for i=1:numel(ORDER)
                        % If the model is in this AllCells in sample
                        if ismember([Constants.other_tag DataType], ...
                                fieldnames(app.data.(smpls{sample_i}).AllCells))
                            if ORDER(i) ~= i-1
                                % Change the elements the indeces of the elements
                                ROWAC(IND_AC(:, i)==1) = ORDER(i);
                            end
                            app.data.(smpls{sample_i}).AllCells.([Constants.other_tag ModelName{1}]) = ROWAC;
                        end
                        % If the model is in MFIRSN
                        if ismember('MFIRSN', fieldnames(app.data.(smpls{sample_i})))
                            if ismember([Constants.other_tag DataType], ...
                                    fieldnames(app.data.(smpls{sample_i}).MFIRSN))
                                if ORDER(i) ~= i-1
                                    ROWRSN(IND_RSN(:, i)==1) = ORDER(i);
                                end
                                app.data.(smpls{sample_i}).MFIRSN.([Constants.other_tag ModelName{1}]) = ROWRSN;
                            end
                        end
                        % If the model is in MFICCN
                        if ismember('MFICCN', fieldnames(app.data.(smpls{sample_i})))
                            if ismember([Constants.other_tag DataType], ...
                                    fieldnames(app.data.(smpls{sample_i}).MFICCN))
                                if ORDER(i) ~= i-1
                                    ROWCCN(IND_CCN(:, i)==1) = ORDER(i);
                                end
                                app.data.(smpls{sample_i}).MFICCN.([Constants.other_tag ModelName{1}]) = ROWCCN;
                            end
                        end
                    end % end second loop through region numbers
                end % end loop through samples
            end % end if there were redefined numbers
        end

        function func_EditChannels(app, web)
            % FUNC_EDITCHANNELS It is a function which allows the user
            % to edit the previous channels and then turn the result to
            % a new channel with given name. 
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            if nargin<2
                web=0;
            end
            alpha = app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            %% Make table options
            % make sure there are actually channels to modify
            % Pull channels
            smples = fieldnames(app.data);
            % build region names
            tDatasm = cell(size(smples,1),2);
            tDatasm(1:end,1) = {false};
            tDatasm(1:end,2) = smples;
            %% Make a table figure
            UIfigsm = uifigure('Name', 'Choose Samples to Edit Channels');
            if web==1
                UIfigsm.Visible='OFF';
            end
            UIfigsm.Position = alpha*[10 10 240 270];

            t = uitable(UIfigsm);
            t.Data = tDatasm;
            t.Position = alpha.*[0 70 240 200];
            t.ColumnName = {'Choose','SamplesName'};
            t.ColumnEditable = [true true true true true];
            t.ColumnWidth = {alpha*100, alpha*100, alpha*100, alpha*100, alpha*100};
            
            % Create a push button
            btn = uibutton(UIfigsm,'push', 'ButtonPushedFcn', @(btn,event) edit_channel(app, smples, t.Data));
            btn.Position = alpha*[30, 10, 100, 50];
            btn.Text = 'Edit Selected';
            Helper.func_SetCLR(app, btn, 'button');
            
            function edit_channel(app,smples,tbl)
                IND = [tbl{1:end,1}];
                if sum(IND) == 0
                    return;
                end
                for i=1:size(IND,2)
                    if IND(i)==1
                        frountend_Edit_channels(app,smples,i);
                    end
                end
            end
            
            function frountend_Edit_channels(app,smples, smp_i, web)
                if nargin<4
                    web=0;
                end
                alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
%                 channels = fieldnames(app.data.(smples{smp_i}).AllCells);
                datTMP = app.data.(smples{smp_i}).AllCells;
                channels = datTMP.Properties.VariableNames';
                %% Make table options
                % make sure there are actually channels to modify
                % build channels names
                tDatasub = cell(numel(channels),8);
                tDatasub(:,1) = channels;
                tDatasub(:,2:end) =  num2cell(datTMP{1:7, :}');
                
%                 for i=1:size(datTMP,2)-3
%                     for j=2:8
% %                         tDatasub(i,j) = (num2cell(table2array(datTMP(j-1,i+3))));
%                         tDatasub(i,j) = table2cell(datTMP(j-1,i+3));
%                     end
%                 end
                %% Make a table figure
                UIfigsub = uifigure('Name',['Channels Details for sample: ' smples{smp_i}]);
                if web==1
                    UIfigsub.Visible='OFF';
                end

                UIfigsub.Position = alpha*[0 0 710 800];
                tch = uitable(UIfigsub);
                tch.Data = tDatasub;
                tch.Position = alpha*[0 70 710 730];
                tch.ColumnName = {'Channel Names','1','2','3','4','5','6','7'};
                tch.ColumnWidth = {alpha*100, alpha*80, alpha*80, alpha*80, alpha*80, alpha*80, alpha*80, alpha*80};
               
                label = uilabel(UIfigsub);
                label.Text = 'New Channel Name';
                label.FontSize = 14;
                label.Position = alpha*[20 35 180 35];
                
                edt = uieditfield(UIfigsub,'text');
                edt.Value = 'NewChannelName';
                edt.Position = alpha*[20 12 180 30];
                edt.Visible = 'ON';

                label2 = uilabel(UIfigsub);
                label2.Text = 'Equation, example: C.Gate_1+C.Gate_2';
                label2.FontSize = 14;
                label2.Position = alpha*[210 35 380 35];
                
                edt2 = uieditfield(UIfigsub,'text');
                edt2.Value = 'C.Gate_1+C.Gate_2';
                edt2.Position = alpha*[210 12 380 30];
                edt2.Visible = 'ON';
                
                % Create a confirm colors push button
                btnsub = uibutton(UIfigsub,'push', 'ButtonPushedFcn', @(btn,event) Plt_Helper.Edit_channels_backend(app,smples{smp_i},edt.Value,edt2.Value,tch));
                btnsub.Position = alpha*[600, 10, 100, 50];
                btnsub.Text = 'Confirm';
                Helper.func_SetCLR(app, btnsub, 'button') 
            end
        end
        
        function Edit_channels_backend(app, smp, NewChannelName, equation, tch)
            %% Pull all of the data
            NewChannelName = Helper.valid_channel(NewChannelName);
            C = app.data.(smp).AllCells;
            channels = C.Properties.VariableNames';
            % equation = 'C.Gate_1+C.Gate_2';
            if ismember(NewChannelName, channels)
                m = questdlg([NewChannelName ' already exists. Do you want to replace this existing channel?'] , '', 'Yes', 'No', 'No');
                if strcmp(m, 'No')
                    return
                end
                app.data.(smp).AllCells.(NewChannelName) = eval(equation);
            else
                app.data.(smp).AllCells.(NewChannelName) = eval(equation);
            end
            % re-populate the table
            C = app.data.(smp).AllCells;
            channels = C.Properties.VariableNames';
            % make sure there are actually channels to modify
            % build channels names
            tData = cell(numel(channels),8);
            tData(:,1) = channels;
            tData(:,2:end) =  num2cell(C{1:7, :}');
            tch.Data = tData;
        end
        
        function func_EditNeighbors(app, web)
            % FUNC_EDITCHANNELS It is a function which allows the user
            % to edit the previous channels and then turn the result to
            % a new channel with given name. 
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            if nargin<2
                web=0;
            end
            alpha = app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            %% Make table options
            % make sure there are actually channels to modify
            % Pull channels
            smples = fieldnames(app.data);
            % build region names
            tDatasm = cell(size(smples,1),2);
            tDatasm(1:end,1) = {false};
            tDatasm(1:end,2) = smples;
            %% Make a table figure
            UIfigsm = uifigure('Name', 'Choose Samples to Edit Neighborhoods');
            if web==1
                UIfigsm.Visible='OFF';
            end
            UIfigsm.Position = alpha*[10 10 240 270];

            t = uitable(UIfigsm);
            t.Data = tDatasm;
            t.Position = alpha.*[0 70 240 200];
            t.ColumnName = {'Choose','SamplesName'};
            t.ColumnEditable = [true true true true true];
            t.ColumnWidth = {alpha*100, alpha*100, alpha*100, alpha*100, alpha*100};
            
            % Create a push button
            btn = uibutton(UIfigsm,'push', 'ButtonPushedFcn', @(btn,event) edit_channel(app, smples, t.Data));
            btn.Position = alpha*[30, 10, 100, 50];
            btn.Text = 'Edit Selected';
            Helper.func_SetCLR(app, btn, 'button');
            
            function edit_channel(app,smples,tbl)
                IND = [tbl{1:end,1}];
                if sum(IND) == 0
                    return;
                end
                for i=1:size(IND,2)
                    if IND(i)==1
                        frountend_Edit_neighbors(app,smples,i);
                    end
                end
            end
            
            function frountend_Edit_neighbors(app,smples, smp_i, web)
                if nargin<4
                    web=0;
                end
                alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
%                 channels = fieldnames(app.data.(smples{smp_i}).AllCells);
                datTMP = app.data.(smples{smp_i}).MFIRSN;
                channels = datTMP.Properties.VariableNames';
                %% Make table options
                % make sure there are actually channels to modify
                % build channels names
                tDatasub = cell(numel(channels),8);
                tDatasub(:,1) = channels;
                tDatasub(:,2:end) =  num2cell(datTMP{1:7, :}');
                %% Make a table figure
                UIfigsub = uifigure('Name',['Neighborhood Details for sample: ' smples{smp_i}]);
                if web==1
                    UIfigsub.Visible='OFF';
                end

                UIfigsub.Position = alpha*[0 0 710 800];
                tch = uitable(UIfigsub);
                tch.Data = tDatasub;
                tch.Position = alpha*[0 70 710 730];
                tch.ColumnName = {'Channel Names','1','2','3','4','5','6','7'};
                tch.ColumnWidth = {alpha*100, alpha*80, alpha*80, alpha*80, alpha*80, alpha*80, alpha*80, alpha*80};
               
                label = uilabel(UIfigsub);
                label.Text = 'New Channel Name';
                label.FontSize = 14;
                label.Position = alpha*[20 35 180 35];
                
                edt = uieditfield(UIfigsub,'text');
                edt.Value = 'NewChannelName';
                edt.Position = alpha*[20 12 180 30];
                edt.Visible = 'ON';

                label2 = uilabel(UIfigsub);
                label2.Text = 'Equation, example: C.Gate_1+C.Gate_2';
                label2.FontSize = 14;
                label2.Position = alpha*[210 35 380 35];
                
                edt2 = uieditfield(UIfigsub,'text');
                edt2.Value = 'C.Gate_1+C.Gate_2';
                edt2.Position = alpha*[210 12 380 30];
                edt2.Visible = 'ON';
                
                % Create a confirm colors push button
                btnsub = uibutton(UIfigsub,'push', 'ButtonPushedFcn', @(btn,event) Plt_Helper.Edit_neighbors_backend(app,smples{smp_i},edt.Value,edt2.Value,tch));
                btnsub.Position = alpha*[600, 10, 100, 50];
                btnsub.Text = 'Confirm';
                Helper.func_SetCLR(app, btnsub, 'button') 
            end
        end
        
        function Edit_neighbors_backend(app, smp, NewChannelName, equation, tch)
            %% Pull all of the data
            NewChannelName = Helper.valid_channel(NewChannelName);
            C = app.data.(smp).MFIRSN;
            channels = C.Properties.VariableNames';
            % equation = 'C.Gate_1+C.Gate_2';
            if ismember(NewChannelName, channels)
                m = questdlg([NewChannelName ' already exists. Do you want to replace this existing channel?'] , '', 'Yes', 'No', 'No');
                if strcmp(m, 'No')
                    return
                end
                app.data.(smp).MFIRSN.(NewChannelName) = eval(equation);
            else
                app.data.(smp).MFIRSN.(NewChannelName) = eval(equation);
            end
            % re-populate the table
            C = app.data.(smp).MFIRSN;
            channels = C.Properties.VariableNames';
            % make sure there are actually channels to modify
            % build channels names
            tData = cell(numel(channels),8);
            tData(:,1) = channels;
            tData(:,2:end) =  num2cell(C{1:7, :}');
            tch.Data = tData;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function func_Spots(app, web)
            % FUNC_SPOTS Creates interface to look at the spots on a
            % current figure. It allows to both add, plot, and delete
            % spots.
            % 
            % Input:
            %   - app - Instance of CytoMAP
            %
            % Modifies:
            %   - app - Adds app.points.PNTS(x, y, z) field, if spots are
            %       saved.
            %
            % Note:
            %   - Currently spots are held globally across all figures, and
            %       saving overwrites all of the previous spots.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            hold on
            [PNTSx, PNTSy] = getpts(gca);
            plot(PNTSx, PNTSy, '.r', 'MarkerSize', 20, 'DisplayName', 'Selected Points')
            options = zeros(numel(PNTSx),3);
            options(:, 1) = PNTSx;
            options(:, 2) = PNTSy;
            % Replace this with the actual Z range/plane these points were place in
            if iscell(options(:, 3))
                options(:, 3) = {0};
            else
                 options(:, 3) = 0;
            end
            UIfig = uifigure('Name', 'Points Options');
            if web==1
                UIfig.Visible='OFF';
            end

            UIfig.Position(3:4) = [350 250];
            % Create the table of options
            t = uitable(UIfig);
            t.Data = options;
            t.Position = alpha*[0 40 350 210];
            t.ColumnName = {'X', 'Y', 'Z'};
            t.ColumnEditable = [true, true, true];
            t.ColumnWidth = {alpha*100, alpha*100, alpha*100};
            % Save/Exit
            btn = uibutton(UIfig, 'push', 'ButtonPushedFcn', @(~,~) Plt_Helper.spots_save(t.Data(:, 1), t.Data(:, 2), t.Data(:, 3), app, UIfig));
            btn.Position = alpha*[340-100-100, 5, 100, 30];
            btn.Text = 'Save and exit';
            Helper.func_SetCLR(app, btn, 'button')
            % Exit
            btn = uibutton(UIfig, 'push', 'ButtonPushedFcn', @(~,~) close(UIfig));
            btn.Position = alpha*[340-100, 5, 100, 30];
            btn.Text = 'Cancel';
            Helper.func_SetCLR(app, btn, 'button')
            % Plot
            btn = uibutton(UIfig, 'push', 'ButtonPushedFcn', @(~,~) Plt_Helper.spots_plot(t.Data(:, 1), t.Data(:, 2), t.Data(:, 3)));
            btn.Position = alpha*[340-300, 5, 100, 30];
            btn.Text = 'Plot';
            Helper.func_SetCLR(app, btn, 'button')
        end

        function spots_plot(x, y, z)
            % SPOTS_PLOT Plots given points, which then can be defined as spots.
            
            plot3( ...
                x, y, z, '.r', ...
                'MarkerSize', 20, ...
                'DisplayName', 'Selected Points' ...
            );
        end

        function spots_save(x, y, z, app, UIfig)
            % SPOTS_SAVE Saves spots in current CytoMAP instance, and
            % closes spots figure
            
            app.points.PNTSx = x;
            app.points.PNTSy = y;
            app.points.PNTSz = z;
            if isa(UIfig, 'matlab.ui.Figure') && isvalid(UIfig)
                close(UIfig)
            end
        end

        function plotpts(app, num)
            % PLOTPTS Plots spots currently existing in the given CytoMAP instance
            
            hold on
            if app.figOpts.D3.(num)
                plot3(app.points.PNTSx,app.points.PNTSy,app.points.PNTSz, '.r', 'MarkerSize', 20, 'DisplayName', 'Selected Points')
            else
                plot(app.points.PNTSx,app.points.PNTSy, '.r', 'MarkerSize', 20, 'DisplayName', 'Selected Points')
            end
            hold off
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot settings changes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plt3D(btn, app, num)
            % PLT3D Toggles the view in a current figure either to 3D or 2D.
            
            if strcmp(btn.State, 'on')
                app.figOpts.D3.(num) = true;
                view(3)
            else
                app.figOpts.D3.(num) = false;
                view(2)
            end

            app.figOpts.xaxIFTMP.(num) = 'NewFig';
            app.figOpts.yaxIFTMP.(num) = 'NewFig';
            app.figOpts.zaxIFTMP.(num) = 'NewFig';
            app.figOpts.caxIFTMP.(num) = 'NewFig';
            Plotting.func_plot(app, num)
        end

        function func_refresh(app, num)
            % FUNC_REFRESH Completely refreshes the current figure. May be used to
            % remove all of drawing on top of the figure (such as spots,
            % polygons etc.)
            
            Plotting.func_plot(app, num);
            if isfield(app.figOpts, 'phn_smpl_tab') && isfield(app.figOpts.phn_smpl_tab, num)
                app.figOpts.phn_smpl_tab.(num).Data = ...
                    Plotting.update_phn_smpl_table(app, num, app.figOpts.phn_smpl_tab.(num).Data);
            end
            % Fix the figure size bob
            Plotting.resize_figure(gcf, app)
        end

        function pltax(btn, app, num)
            % PLTAX Toggles the view in a current figure to square axis (i.e.
            % axes are unscaled and true to the data), or normalized axis
            % (i.e. axis are scaled to fit plot nicely).
            
            if strcmp(btn.State, 'on')
                axis equal
                app.figOpts.eqax.(num) = true;
            else
                axis normal
                app.figOpts.eqax.(num) = false;
            end
        end

        function pltIC(btn, app, num)
            % PLTIC Inverts theme of plot. Can be either dark or light themed.
            
            if strcmp(btn.State, 'on')
                app.figOpts.IC.(num) = true;
                app.figOpts.bgclr.(num) = 'k';
                app.figOpts.txtclr.(num) = 'w';
            else
                app.figOpts.IC.(num) = false;
                app.figOpts.bgclr.(num)  = 'w';
                app.figOpts.txtclr.(num) = 'k';
            end
            figNUM = gcf;
            figNUM.Color = app.figOpts.bgclr.(num);
            figNUM = gca;
            figNUM.Color = app.figOpts.bgclr.(num);
            figNUM.XColor = app.figOpts.txtclr.(num);
            figNUM.YColor = app.figOpts.txtclr.(num);
            figNUM.ZColor = app.figOpts.txtclr.(num);
        end

        function plot_options(app, num, web)
            % PLOT_OPTIONS Allows user to change different things about
            % given plot, such as axis scale, ranges etc.
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            %   - num - char, string - Identifier of figure which will be
            %       changed
            %
            % Modifies:
            %   - app - More specifically multiple fields in app.figOpts,
            %       which describe axes of the plot.
            %   - given plot - How it looks, not the data itself.
            
            % Build the options menus
            if nargin<3
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            ax2 = gca;
            app.figOpts.limits.(num) = cell(4,4);
            if ~isempty(ax2.Colorbar)
                app.figOpts.limits.(num)(1:4, 1) = {ax2.XLim(1), ax2.YLim(1), ax2.ZLim(1), ax2.Colorbar.Limits(1)};
                app.figOpts.limits.(num)(1:4, 2)= {ax2.XLim(2), ax2.YLim(2), ax2.ZLim(2), ax2.Colorbar.Limits(2)};
            else
                app.figOpts.limits.(num)(1:3, 1) = {ax2.XLim(1), ax2.YLim(1),ax2.ZLim(1)};
                app.figOpts.limits.(num)(1:3, 2)= {ax2.XLim(2), ax2.YLim(2),ax2.ZLim(2)};
            end
            app.figOpts.limits.(num)(:, 3) = {false};
            app.figOpts.limits.(num)(:, 4) = {false};
            app.figOpts.limits.(num)(4, 5) = {false};
            app.figOpts.limits.(num)(1:3, 5) = num2cell(strcmp({ax2.XDir, ax2.YDir, ax2.ZDir}, 'reverse'));
            UIfig = uifigure('Name', 'Plot Options');
            if web==1
                UIfig.Visible='OFF';
            end

            UIfig.Position(3:4) = [370 200];
            % Create the table of options
            t = uitable(UIfig);
            t.Data = app.figOpts.limits.(num);
            t.Position = alpha*[0 40 370 160];
            t.ColumnName = {'Min','Max','Auto', 'Log', 'Invert'};
            t.RowName = {'X-Axis', 'Y-Axis', 'Z-Axis', 'Color-Axis'};
            t.ColumnEditable = [true, true, true, true, true];
            t.ColumnWidth = {75 75 50, 50, 50};
            
            % Add change colormap dropdown
            clrmapopt = uidropdown(UIfig);
            clrmapopt.Position = alpha*[5+100+5 5 150 30];
            clrmapopt.Items = {'Default Colormap', ...
                               'jet', ...
                               'redblue', ...
                               'red', ...
                               'bone', ...
                               'autumn', ...
                               'prism', ...
                               'hsv', ...
                               'hot'};
            clrmapopt.ValueChangedFcn = @(~, ~) Plt_Helper.plot_options_backend( ...
                app, UIfig, num, ...
                cell2mat(t.Data(1:4, 3)), ... Whether ranges are set to automode, or manual
                cell2mat(t.Data(1:4, 1:2)), ... Manual ranges
                cell2mat(t.Data(1:4, 4)), ... Whether scale if log or linear
                cell2mat(t.Data(1:4, 5)), ... Invert the selected axis
                [], clrmapopt);
            clrmapopt.Value = clrmapopt.Items{1};
            Helper.func_SetCLR(app, clrmapopt, 'button');
            
            
            t.CellEditCallback = @(dd, p) Plt_Helper.plot_options_backend( ...
                app, UIfig, num, ...
                cell2mat(dd.Data(1:4, 3)), ... Whether ranges are set to automode, or manual
                cell2mat(dd.Data(1:4, 1:2)), ... Manual ranges
                cell2mat(dd.Data(1:4, 4)), ... Whether scale if log or linear
                cell2mat(dd.Data(1:4, 5)), ... Invert the selected axis
                p, clrmapopt);
            
            % Create an edit pot points button
            btn = uibutton(UIfig,'push');
            btn.Position = alpha*[5+100+5+150+5, 5, 100, 30];
            btn.Text = 'Edit Points';
            btn.ButtonPushedFcn = @(btn,event) Plt_Helper.edit_point_colors(app);
            Helper.func_SetCLR(app, btn, 'button')
            
            % Create an apply options to any selected figure button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Plt_Helper.plot_options_backend( ...
                app, UIfig, num, ...
                cell2mat(t.Data(1:4, 3)), ... Whether ranges are set to automode, or manual
                cell2mat(t.Data(1:4, 1:2)), ... Manual ranges
                cell2mat(t.Data(1:4, 4)), ... Whether scale if log or linear
                cell2mat(t.Data(1:4, 5)), ... Invert the selected axis
                [], clrmapopt));
            btn.Position = alpha*[5, 5, 100, 30];
            btn.Text = 'Apply To Figure';
            Helper.func_SetCLR(app, btn, 'button')
        end

        function plot_options_backend(app, UIfig, num, automodes, ranges, scales, invert, p, clrmapopt)
            % PLOT_OPTIONS_BACKEND Changes what are the ranges, modes, and
            % scales of axes on the given figure
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            %   - UIfig - Handle to plot options front-end which will be
            %       closed after this function executes.
            %   - num - char, string - Identifier of figure which will be
            %       changed.
            %   - autmodes - logical, numerical - Array of size (4, 1) of
            %       truths/falses, on whether to apply automatic or manual
            %       ranges of axes (x, y, z, c).
            %   - ranges - numerical - Array of size (4, 2) of min, max
            %       rows which defines manual ranges of axes.
            %   - scales - logical, numerical - Array of size (4, 1) of
            %       truths/falses, on whether to apply logarithmic or
            %       linear scale on axes (x, y, z, c).
            %
            % Modifies:
            %   - app - More specifically multiple fields in app.figOpts,
            %       which describe axes of the plot.
            %   - given plot - How it looks, not the data itself.
            
            % This should set all of the plot options that were
            % selected
            ax3 = gca;
            % if the user selected the dropdown
            if isempty(p)
                % change colormap scheme
                switch clrmapopt.Value
                    case 'Default Colormap'

                    case 'redblue'
                        ax3.Colormap = app.GUIOPTS.redbluecmap;
                    case 'red'
                        cmap = app.GUIOPTS.redbluecmap;
                        cmap(1:round(size(cmap,1)/2), :) = [];
                        cmap = [1,1,1; cmap];
                        ax3.Colormap = cmap;
                    otherwise
                        ax3.Colormap = eval([clrmapopt.Value '();']);
                end
            end
            
            if automodes(1)
                ax3.XLimMode = 'auto';
            else
                ax3.XLim = ranges(1, :);
            end
            if automodes(2)
                ax3.YLimMode = 'auto';
            else
                ax3.YLim = ranges(2, :);
            end
            if automodes(3)
                ax3.ZLimMode = 'auto';
            else
                ax3.ZLim = ranges(3, :);
            end

            if scales(1)
                ax3.XScale= 'log';
            else
                ax3.XScale= 'linear';
            end
            if scales(2)
                ax3.YScale= 'log';
            else
                ax3.YScale= 'linear';
            end
            if scales(3)
                ax3.ZScale= 'log';
            else
                ax3.ZScale= 'linear';
            end

            if ~isempty(ax3.Colorbar)
                if automodes(4)
                    ax3.Colorbar.LimitsMode = 'auto';
                else
                    limits = ranges(4, :);
                    caxis(limits);
                    ax3.Colorbar.Limits = limits;
                    ax3.Colorbar.Ticks = limits(1):((limits(2)-limits(1))/10):limits(2);
                end
                if scales(4)
                    set(ax3,'colorscale','log')
                else
                    set(ax3,'colorscale','linear')
                end
            end
            if invert(1)
                set(ax3, 'xdir', 'reverse');
            else
                set(ax3, 'xdir', 'normal');
            end
            if invert(2)
                set(ax3, 'ydir', 'reverse');
            else
                set(ax3, 'ydir', 'normal');
            end
            if invert(3)
                set(ax3, 'zdir', 'reverse');
            else
                set(ax3, 'zdir', 'normal');
            end
            if ~isempty(p)
                if p.Indices == [4, 5] 
                    oldcmap = ax3.Colormap;
                    ax3.Colormap = flipud(oldcmap);
                    
                end                
            end
        end

        function edit_point_colors(app)
            axh = gca;
            Dispnms = {axh.Children.DisplayName};
            clrs = cell(numel(Dispnms), 3);
            for ch_i = 1:numel(Dispnms)
                switch axh.Children(ch_i).CDataMode
                    case 'auto'
                        clrs(ch_i, :) = num2cell(axh.Children(ch_i).CData);
                    case 'manual'
                        %need to put something in here
                end
            end            
            %% Build a table with each plotted element as a row
            dat = cell(numel(Dispnms),7);
            dat(:, 1) = {true};
            dat(:, 2) = {axh.Children.SizeData};
            dat(:, 3) = {axh.Children.Marker};
            dat(:, 4) = {axh.Children.MarkerFaceAlpha};
            dat(:, 5:7) = clrs;
            
            
            UIfig_points = uifigure('Name', 'Edit Points');
            UIfig_points.Position(3:4) = [700 200];
            % Create the table of options
            t = uitable(UIfig_points);
            t.Data = dat;
            t.Position = [0 0 700 200];
            t.ColumnName = {'Visible', 'Size','Symbol', 'Opacity', 'Color R', 'Color G', 'Color B'};
            t.RowName = Dispnms;
            t.ColumnEditable = true;
            t.ColumnWidth = {55, 50, 60, 75,75,75,75};
            
            t.CellEditCallback = @(dd, p) Plt_Helper.plot_points_backend(...
                app, ...
                cell2mat(dd.Data(:, 1)), ... visible points
                cell2mat(dd.Data(:, 2)), ... point sizes
                cell2mat(dd.Data(:, 3)), ... what symbol to use
                cell2mat(dd.Data(:, 4)), ... marker opacity
                cell2mat(dd.Data(:, 5:7)) ...what colors
                );
        end
        
        function plot_points_backend(app, visible, size, symbol, opacity, color)
            axh = gca;
            for ch_j = 1:numel({axh.Children.DisplayName})
                switch visible(ch_j)
                    case 1
                        axh.Children(ch_j).Visible = 'on';
                    case 0
                        axh.Children(ch_j).Visible = 'off';
                end
                axh.Children(ch_j).SizeData = size(ch_j);
                axh.Children(ch_j).Marker = symbol(ch_j);
                switch axh.Children(ch_j).CDataMode
                    case 'auto'
                        axh.Children(ch_j).CData = color(ch_j, :);
                    case 'manual'
                        %need to put something in here
                        try
                            axh.Children(ch_j).CData = color(ch_j, :); 
                        catch
                            
                        end
                end
                
                axh.Children(ch_j).MarkerEdgeAlpha = opacity(ch_j);
                axh.Children(ch_j).MarkerFaceAlpha = opacity(ch_j);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rotation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function frontend_rotate(app, web)
            % FRONTEND_ROTATE Front-End to rotating plot around an axis, by some user
            % specified angle.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            options = cell(1,4);
            options(1,1:2) ={0};
            options(1,3) = {1};
            options(1,4) = {5};

            UIfig = uifigure('Name', 'Rotate Plot');
            if web==1
                UIfig.Visible='OFF';
            end

            UIfig.Position(3:4) = [375 250];

            % Create the table of options
            t = uitable(UIfig);
            t.Data = options;
            t.Position = alpha*[0 40 375 210];
            t.ColumnName = {'X axis', 'Y axis', 'Z axis', 'Rotation Angle'};
            t.ColumnEditable = [true, true, true, true];
            t.ColumnWidth = {alpha*50, alpha*50, alpha*50, alpha*150};

            % Plot/Edit Surface
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Plt_Helper.rotate_plot(t.Data{1,4}, cell2mat(t.Data(1, 1:3))));
            btn.Position = alpha*[10, 5, 175, 30];
            btn.Text = 'Rotate about selected axis';
            Helper.func_SetCLR(app, btn, 'button')
        end
        
        function rotate_plot(alpha, azel)
            % ROTATE_PLOT Look at docs below. azel = [THETA PHI]
            %
            % We are using built-in matlab function here, but with support
            % for scatter plots. Also we fixed their indetation. Like who
            % uses a mix of 2 and 3 spaces?
            %
            % ROTATE Rotate objects about specified origin and direction.
            %   ROTATE(H,[THETA PHI],ALPHA) rotates the objects with handles H
            %   through angle ALPHA about an axis described by the 2-element
            %   direction vector [THETA PHI] (spherical coordinates).
            %   All the angles are in degrees.  The handles in H must be children
            %   of the same axes.
            %
            %   THETA is the angle in the xy plane counterclockwise from the
            %   positive x axis.  PHI is the elevation of the direction vector
            %   from the xy plane (see also SPH2CART).  Positive ALPHA is defined
            %   as the righthand-rule angle about the direction vector as it
            %   extends from the origin.
            %
            %   ROTATE(H,[X Y Z],ALPHA) rotates the objects about the direction
            %   vector [X Y Z] (Cartesian coordinates). The direction vector
            %   is the vector from the center of the plot box to (X,Y,Z).
            %
            %   ROTATE(...,ORIGIN) uses the point ORIGIN = [x0,y0,z0] as
            %   the center of rotation instead of the center of the plot box.
            %
            %   See also SPH2CART, CART2SPH.
            %
            %   Copyright 1984-2017 The MathWorks, Inc.
            
            h = gca;
            h = h.Children;
            % Determine the default origin (center of plot box).
            if ~isgraphics(h)
                error(message('MATLAB:rotate:InvalidHandle'));
            end
            ax = ancestor(h(1),'axes');
            if isempty(ax) || ax==0
                error(message('MATLAB:rotate:InvalidHandle'));
            end
            origin = sum([get(ax,'xlim')' get(ax,'ylim')' get(ax,'zlim')'])/2;

            % find unit vector for axis of rotation
            if numel(azel) == 2 % theta, phi
                theta = pi*azel(1)/180;
                phi = pi*azel(2)/180;
                u = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
            elseif numel(azel) == 3 % direction vector
                u = azel(:)/norm(azel);
            end

            alph = alpha*pi/180;
            cosa = cos(alph);
            sina = sin(alph);
            vera = 1 - cosa;
            x = u(1);
            y = u(2);
            z = u(3);
            rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
                   x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
                   x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';

            for i=1:numel(h)
                t = get(h(i),'type');
                skip = 0;
                if strcmp(t,'surface') || strcmp(t,'line') || strcmp(t,'patch') || strcmp(t,'scatter')
                    % If patch, rotate vertices
                    if strcmp(t,'patch')
                        verts = get(h(i),'Vertices');
                        x = verts(:,1); y = verts(:,2);
                        if size(verts,2)>2
                            z = verts(:,3);
                        else
                            z = [];
                        end

                    % If surface or line, rotate {x,y,z}data
                    else
                        x = get(h(i),'xdata');
                        y = get(h(i),'ydata');
                        z = get(h(i),'zdata');
                    end

                    if isempty(z)
                        z = -origin(3)*ones(size(y));
                    end
                    [m,n] = size(z);
                    if numel(x) < m*n
                        [x,y] = meshgrid(x,y);
                    end
                elseif strcmp(t,'text')
                    p = get(h(i),'position');
                    x = p(1); y = p(2); z = p(3);
                elseif strcmp(t,'image')
                    x = get(h(i),'xdata');
                    y = get(h(i),'ydata');
                    z = zeros(size(x));
                else
                    skip = 1;
                end

                if ~skip
                    [m,n] = size(x);
                    newxyz = [x(:)-origin(1), y(:)-origin(2), z(:)-origin(3)];
                    newxyz = newxyz*rot;
                    newx = origin(1) + reshape(newxyz(:,1),m,n);
                    newy = origin(2) + reshape(newxyz(:,2),m,n);
                    newz = origin(3) + reshape(newxyz(:,3),m,n);

                    if strcmp(t,'surface') || strcmp(t,'line') || strcmp(t,'scatter')
                        set(h(i),'xdata',newx,'ydata',newy,'zdata',newz);
                    elseif strcmp(t,'patch')
                        set(h(i),'Vertices',[newx,newy,newz]);
                    elseif strcmp(t,'text')
                        set(h(i),'position',[newx newy newz])
                    elseif strcmp(t,'image')
                        set(h(i),'xdata',newx,'ydata',newy)
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Random Points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function frontend_random_pts(app, num, web)
            % FRONTEND_RANDOM_PTS It is a user interface for generating,
            % plotting and saving random points, which can be used as a
            % metric to compare existing distributions to random ones.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - num - char, str - Figure indetifier. Decides where to
            %       pull data from.
            %
            % Modifies:
            %   - app - Can modify currently chosen samples by adding new
            %       gate corresponding to random points in it.
            %   - current plot - Can make a new plot of random points.
            if nargin<3
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

            %% Build the options menus
            options = cell(5,3);

            % Pull the currently plotted data points
            %% Get plotted data
            dat = app.figOpts.dat.(num);
            vrnms = app.figOpts.vrnms.(num);
            % Do not use UpdateFigMenu unless you really need to
            if ~istable(dat)
                if numel(dat) == 1
                    dat = dat{1};
                    app.figOpts.dat.(num) = dat;
                else
                    [dat, vrnms, ~] = Plt_Helper.UpdateFigMenu(app, num, 'Combine', true);
                    app.figOpts.dat.(num) = dat;
                    app.figOpts.vrnms.(num) = vrnms;
                end
            end
            XRange = xlim; YRange = ylim; ZRange = zlim;

            options(:,1) = {XRange(1), YRange(1), ZRange(1), '', ''};
            options(:,2) = {XRange(2), YRange(2), ZRange(2), size(dat, 1), 'New_Random_Points'};
            options(:,3) = {'X range', 'Y Range', 'Z Range', 'Number of Cells', 'Export Name'};
            UIfig = uifigure('Name', 'Randomize Location');
            if web==1
                UIfig.Visible='OFF';
            end
            UIfig.Position(3:4) = alpha*[375 250];
            % Create the table of options
            t = uitable(UIfig);
            t.Data = options;
            t.Position = alpha*[0 40 375 210];
            t.ColumnName = {'Min','Max', 'Option'};
            t.ColumnEditable = [true, true, false];
            t.ColumnWidth = {alpha*115, alpha*125, alpha*120};

            %% Edit number of points
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Plt_Helper.backend_random_pts( ...
                app, vrnms, num, ...
                t.Data{4,2}, ... Number of random points
                cell2mat(t.Data(1:3, 1:2)), ... Ranges on which to apply the random points
                false, ... Whether to add those random points to sample or just plot them
                t.Data{5,2} ... Name of the gate representing random points
            ));
            btn.Position = alpha*[375-85, 5, 75, 30];
            btn.Text = 'Plot';
            Helper.func_SetCLR(app, btn, 'button')

            %% add random points to data
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Plt_Helper.backend_random_pts( ...
                app, vrnms, num, ...
                t.Data{4,2}, ... Number of random points
                cell2mat(t.Data(1:3, 1:2)), ... Ranges on which to apply the random points
                true, ... Whether to add those random points to sample or just plot them
                t.Data{5,2} ... Name of the gate representing random points
            ));
            btn.Position = alpha*[10, 5, 115, 30];
            btn.Text = 'Export Population';
            Helper.func_SetCLR(app, btn, 'button')
        end

        function backend_random_pts(app, vrnms, num, num_rand_pts, ranges, export, gate_name)
            % BACKEND_RANDOM_PTS It is backend for generating, plotting and
            % saving random points, which can be used as metric to compare
            % existing distributions to random ones.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - dat - table - Table of currently plotted things
            %   - vrnms - cell of char/string - Things which are legal on
            %       current axes.
            %   - num - char, string - Figure indetifier. Decides where to
            %       pull data from.
            %   - num_rand_pts - int - Number of random points to be
            %       generated. Has to be <= size(dat, 1)
            %   - ranges - numerical matrix of size 3x2 - Contains ranges
            %       on which random points should be generated. Every row
            %       will correspond to 1st, 2nd, 3rd axis.
            %   - export - boolean - Whether to export random points to all
            %       currently plotted samples (true) or just plot them
            %       (false).
            %   - gate_name - cell, string, char - If export is true, then
            %       it's name of gate which will be created in all
            %       currently plotted samples.
            %
            % Modifies:
            %   - app.data.(smpl) - Where smpl is every sample in
            %       app.figOpts.samples.(num). It will be added a new gate
            %       in all things keeping track of gates. If export is
            %       true.
            %   - current plot - Can make a new plot of random points. If
            %       export is false.
            
            % Zero out all elements of the selected table
            % (This is my sloppy way of making a zero element table with the correct dimensions)

            % Find the Sample Name
            smpls = app.figOpts.smpls.(num);
            if ~iscell(smpls)
                smpls = {smpls};
            end
            vrnmstmp = app.data.(smpls{1}).AllCells.Properties.VariableNames;
            % Initialize the random points table
% % %             typ = cell(numel(vrnmstmp),1);
% % %             typ(:) = {'double'};
            typ = varfun(@class,app.data.(smpls{1}).AllCells,'OutputFormat','cell');
            dat = table('Size',[num_rand_pts, numel(vrnmstmp)],'VariableTypes',typ,'VariableNames', vrnmstmp);

                % X axis
                dat.(vrnms{app.figOpts.xaxIF.(num).Value}) = ...
                    rand(num_rand_pts, 1) .* (ranges(1, 2) - ranges(1, 1)) + ranges(1, 1);
                % Y axis
                dat.(vrnms{app.figOpts.yaxIF.(num).Value}) = ...
                    rand(num_rand_pts, 1) .* (ranges(2, 2) - ranges(2, 1)) + ranges(2, 1);
                if app.figOpts.D3.(num) % 3D plot
                    % Z axis
                    dat.(vrnms{app.figOpts.zaxIF.(num).Value}) = ...
                        rand(num_rand_pts, 1) .* (ranges(3, 2) - ranges(3, 1)) + ranges(3, 1);
                end

            if export
                % I just want to add a new population to AllCells, Tree,
                % and GateTags and the Phenotype table

                gate_name = ['Rand_' gate_name];
                
                valid_gate_name = Helper.valid_gate(gate_name);
                if iscell(valid_gate_name)
                    valid_gate_name = valid_gate_name{1};
                end
                for smpl_idx = 1:numel(smpls)
                    smpl = smpls{smpl_idx};

                    vrnmstmp = app.data.(smpl).AllCells.Properties.VariableNames;
                    % Initialize the random points table
% %                     typ = cell(numel(vrnmstmp),1);
% %                     typ(:) = {'double'};
                    % Get the variable types
                    typ = varfun(@class,app.data.(smpl).AllCells,'OutputFormat','cell');
                    % Initialize the table
                    datTMP = table('Size',[num_rand_pts, numel(vrnmstmp)],'VariableTypes',typ,'VariableNames', vrnmstmp);
                    
                    % X axis
                    datTMP.(vrnms{app.figOpts.xaxIF.(num).Value}) = dat.(vrnms{app.figOpts.xaxIF.(num).Value});
                    % Y axis
                    datTMP.(vrnms{app.figOpts.yaxIF.(num).Value}) = dat.(vrnms{app.figOpts.yaxIF.(num).Value});
                    if app.figOpts.D3.(num) % 3D plot
                        % Z axis
                        datTMP.(vrnms{app.figOpts.zaxIF.(num).Value}) = dat.(vrnms{app.figOpts.zaxIF.(num).Value});
                    end
                    
                    % Find the Current Number of Gates
                    Gates = app.data.(smpl).GateTags.Properties.VariableNames;
                    NGates = max(str2double(strrep(Gates, Constants.gate_tag, '')));

                    % Add the new population to the tree, right under
                    % the top of the tree.
                    rand_tree = tree(gate_name, 'gate_type', 'logic');
                    app.data.(smpl).tree = app.data.(smpl).tree.add_kid(rand_tree);

                    % Add the new population to GateTags
                    valid_name = [app.data.(smpl).tree.name '/' valid_gate_name];
                    short_name = [app.data.(smpl).tree.full_name '/' char(gate_name)];
                    if iscell(short_name)
                        short_name = short_name{1};
                    end
                    tag = Helper.get_tag(app, short_name, smpl);
                    if iscell(tag)
                        tag = tag{1};
                    end
                    app.data.(smpl).GateTags.(tag) = {valid_name; short_name};
                    % Deal with Cells containing strings
                    INDCell = contains(typ, 'cell');
                    datTMP(:, INDCell) = {'NA'};
                    
                    % Add the new population to AllCells
                    IND = zeros(size(app.data.(smpl).AllCells, 1) + size(datTMP, 1), 1);
                    IND((size(app.data.(smpl).AllCells, 1)+1):end) = ones(size(datTMP, 1), 1);
                    app.data.(smpl).AllCells = [app.data.(smpl).AllCells; datTMP];
% % %                     app.data.(smpl).AllCells.([ Constants.gate_tag, num2str(NGates+1)]) = IND;
                    app.data.(smpl).AllCells.(tag) = IND;
                end
            else
                X = dat.(vrnms{app.figOpts.xaxIF.(num).Value});
                Y = dat.(vrnms{app.figOpts.yaxIF.(num).Value});
                Z = dat.(vrnms{app.figOpts.zaxIF.(num).Value});
                if app.figOpts.D3.(num) % 3D plot
                    hold on
                    scatter3(X, Y, Z, 5, 'filled', 'MarkerFaceColor', 'r')
                    hold off
                else
                    hold on
                    scatter(X, Y, 5, 'filled', 'MarkerFaceColor', 'r')
                    hold off
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Surfaces
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function frontend_surf(app, num, smpls, name0, edit, Type, web)
            % FRONTEND_SURF It is the front-end to surface edition parts of
            % program.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - num - Figure identifier (typically 'fig1' or similar)
            %   - smpls - cell, char, string - samples on which to change
            %       surfaces.
            %   - name0 - char, string - name of the surface to be initally
            %       filled in the uitable generated for user. If it's empty
            %       then defaults to 'New_Surface'.
            %   - edit - single logical - if user wants to define a surface
            %       whether to define it as a new surface or overwrite an
            %       existing one.
            %   - Type - char - Type of surface. For example:
            %       UDS for User Defined Surfaces.
            %       All of different types of surfaces can be found under
            %       app.Surfaces.
            %
            % Modifies:
            %   - app - More specifically field
            %       app.data.(smpl).Surfaces.(Type) for each smpl in smpls
            %       either by adding, removing or editing a surface with
            %       name chosen by user.
            
            if nargin<7
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            %% Build the options menus
            options = cell(9,4);

            if ~iscell(smpls)
                smpls = {smpls};
            end

            if isempty(name0)
                options(:,1) = {100, 50, 50, 'r', 0.1, 'New_Surface', false, 'New_Population', false};
            else
                for s_idx = 1:numel(smpls)
                    s = smpls{s_idx};
                    VolT = app.data.(s).Surfaces.(Type).(name0).Surf;
                end
                options(:,1) = {VolT.Alpha, VolT.RegionThreshold, VolT.HoleThreshold,...
                                'r', 0.1, name0, false, 'New_Population', false};
                if size(VolT.Points, 2)==3
                    options(1:4,3) = {'', '', VolT.surfaceArea, VolT.volume};
                else
                    options(1:4,3) = {VolT.perimeter, VolT.area, '', ''};
                end
            end
            options(1:4, 4) = {'Perimeter', 'Area', 'Surface Area', 'Volume'};
            options(:, 2) = { ...
                'Surface Creation Radius', ...
                'Volume Exclusion Threshold', ...
                'Hole Fill Threshold', ...
                'Color', ...
                'Opacity', ...
                'Surface Name', ...
                'Redefine Surface?', ...
                'Export Name', ...
                'External Cells?' ...
            };

            UIfig = uifigure('Name', 'Surface Options');
            if web==1
                UIfig.Visible='OFF';
            end

            UIfig.Position(3:4) = alpha*[620 250];
            % Create the table of options
            t = uitable(UIfig);
            t.Data = options;
            t.Position = alpha*[0 40 UIfig.Position(3) 210];
            t.ColumnName = {'Value','Option', 'Statistics'};
            t.ColumnEditable = [true, false, false, false];
            t.ColumnWidth = {alpha*100, alpha*210, alpha*100, alpha*210};

            % Plot/Edit Surface
            btn = uibutton(UIfig, 'push', ...
                'ButtonPushedFcn', @(btn,event) Plt_Helper.edit_surf( ...
                    app, ...
                    edit, ...
                    t.Data{7, 1}, ... Redefine a surface boolean
                    num, ...
                    smpls, ...
                    Type, ...
                    t.Data{6, 1}, ... Name of Surface
                    t.Data{1, 1}, ... Alpha of Surface
                    t.Data{2, 1}, ... Region Threshold
                    t.Data{3, 1}, ... Hole Threshold
                    t.Data{4, 1}, ... Color of Face of Surface
                    t.Data{5, 1} ... Alpha of Face of Surface
                ) ...
            );
            btn.Position = alpha*[UIfig.Position(3)-75, 5, 75, 30];
            btn.Text = 'Plot';
            Helper.func_SetCLR(app, btn, 'button')

            % Plot delete surface
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Plt_Helper.del_surf( ...
                app, ...
                smpls, ...
                Type, ...
                t.Data{6, 1} ... Name of Surface
            ));
            btn.Position = alpha*[UIfig.Position(3)-150, 5, 75, 30];
            btn.Text = 'Delete';
            Helper.func_SetCLR(app, btn, 'button')

            % Gate on cells within surface
            btn = uibutton(UIfig, 'push', ...
                'ButtonPushedFcn', @(btn,event) Plt_Gating.save_surf_gate( ...
                    app, ...
                    num, ...
                    smpls, ...
                    Type, ...
                    t.Data{6, 1}, ... Surface Name
                    t.Data{8, 1}, ... Population Name
                    t.Data{9, 1} ... Include things outside and not inside
                ) ...
            );
            btn.Position = alpha*[10, 5, 115, 30];
            btn.Text = 'Export Population';
            Helper.func_SetCLR(app, btn, 'button')
            
            % Dropdown with Parent Names
            parents = uidropdown(UIfig, ...
                'Items', Helper.get_common_surfaces(app, smpls, Type) ...
            );
            parents.Position = alpha*[125, 5, 115, 30];
            
            % Make subsurface button
            btn = uibutton(UIfig, 'push', ...
                'ButtonPushedFcn', @(btn, event) Plt_Helper.add_sub_surf( ...
                    app, ...
                    num, ...
                    smpls, ...
                    Type, ...
                    t.Data{6, 1}, ... Surface Name
                    t.Data{8, 1}, ... Population Name
                    parents.Value, ... Parent Surface Name
                    t.Data{4, 1}, ... Color of Face of Surface
                    t.Data{5, 1}, ... Alpha of Face of Surface
                    t.Data{9, 1} ... Include things outside and not inside
                ) ...
            );
            btn.Position = alpha*[240, 5, 115, 30];
            btn.Text = 'Make Subsurface';

        end % end draw surfaces function

        function edit_surf(app, edit, redefine, num, smpls, Type, name, alpha, reg_thres, hole_thres, face_color, face_alpha)
            % EDIT_SURF Edits or creates a surface in given samples.
            % Partially a wrap around Plt_Helper.surf_create.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - edit - Whether to edit currently existing surface or not.
            %   - redefine - Whether to redefine, and make a completely new
            %       surface, even if one with the same name exists (true),
            %       or not (false).
            %   - num - Figure identifier on which surface will be plotted,
            %       and from which data will be pulled.
            %   - smpls - Samples on which to define/modify surface
            %   - Type - Type of surface. This determines place to which
            %       surface will be saved. Example: UDS - User Defined
            %       Surface
            %   - name - Name of this particular surface.
            %   - alpha - Transparency of the surface. Value between 0-1.
            %   - reg_thers - Maximum threshold of supressing regions in
            %       the surface. For more details, better go to alphashape
            %       in matlab docs.
            %   - hole_thres - Maximum threshold of filling in holes in the
            %       surface. For more details, better go to alphashape in
            %       matlab docs.
            %   - face_color - Color of the face of the surface.
            %   - face_alpha - Transparency of the face of the surface
            %
            % Modifies:
            %   - app - Adds in a surface struct, under a given sample,
            %       type and name. Fills it in with information needed to
            %       reference, and reconstruct the surface (except for 
            %       region/hole threshold information).
            
            if isempty(name)
                errordlg('Surface name empty. Aborting.');
                return;
            end
            name = Helper.valid_var(name);

% %             smpls = app.figOpts.smpls.(num);
            
            if isempty(smpls)
                smpls = app.figOpts.smpls.(num);
            end
            if ~iscell(smpls)
                smpls = {smpls};
            end

            edit_smpl_check = true;
            for smpl_idx=1:numel(smpls)
                smpl = smpls{smpl_idx};
                if ~isfield(app.data.(smpl).Surfaces.(Type), name) || strcmp(name, 'New_Surface')
                    edit_smpl_check = false;
                end
            end

            if edit_smpl_check
                edit = true;
            end

            % re-define the surface
            if redefine
                edit = false;
                % if you are redefining the surface, check to see if the
                % plotted sample changed
                smpls = app.figOpts.smpls.(num);
                if ~iscell(smpls)
                    smpls = {smpls};
                end
            end

            if edit
                for smpl_idx = 1:numel(smpls)
                    smpl = smpls{smpl_idx};
                    Vol = app.data.(smpl).Surfaces.(Type).(name).Surf;
                    plt = app.data.(smpl).Surfaces.(Type).(name).Plot;
                    if ~isnumeric(app.data.(smpl).Surfaces.(Type).(name).Plot) && isvalid(app.data.(smpl).Surfaces.(Type).(name).Plot)
                        delete(plt)
                    end
                    Vol.Alpha = alpha;
                    Vol.RegionThreshold = reg_thres;
                    Vol.HoleThreshold = hole_thres;

                    LimX = xlim;
                    LimY = ylim;
                    if app.figOpts.D3.(num) % 3D plot
                        LimZ = zlim;
                    end

                    hold on
                    plt = plot(Vol, 'DisplayName', name);
                    hold off
                    xlim(LimX)
                    ylim(LimY)
                    if app.figOpts.D3.(num) % 3D plot
                        zlim(LimZ)
                    end

                    plt.FaceColor = face_color;
                    plt.FaceAlpha = face_alpha;
                    plt.EdgeColor = 'none';

                    app.data.(smpl).Surfaces.(Type).(name).Surf = Vol;
                    app.data.(smpl).Surfaces.(Type).(name).Plot = plt;
                end
            else % This is a new surface
                %% Pull the selected data
                %% Get plotted data
                dat = app.figOpts.dat.(num);
                
                % Do not use UpdateFigMenu unless you really need to
                if ~istable(dat)
                    if numel(dat) == 1
                        dat = dat{1};
                        app.figOpts.dat.(num) = dat;
                    else
                        [dat, vrnms, ~] = Plt_Helper.UpdateFigMenu(app, num, 'Combine', true);
                        app.figOpts.dat.(num) = dat;
                        app.figOpts.vrnms.(num) = vrnms;
                    end
                end
                Xax = app.figOpts.xaxIF.(num).String{app.figOpts.xaxIF.(num).Value};
                Yax = app.figOpts.yaxIF.(num).String{app.figOpts.yaxIF.(num).Value};
                Zax = app.figOpts.zaxIF.(num).String{app.figOpts.zaxIF.(num).Value};

                LimX = xlim;
                LimY = ylim;
                if app.figOpts.D3.(num) % 3D plot
                    LimZ = zlim;
                end
                for smpl_idx = 1:numel(smpls)
                    smpl = smpls{smpl_idx};
                    if ~isfield(app.data.(smpl).Surfaces, Type)
                        app.data.(smpl).Surfaces.(Type) = struct;
                    end

                    Plt_Helper.surf_create(app, smpl, app.figOpts.phenos.(num), dat, Xax, Yax, Zax, alpha, Type, name, reg_thres, hole_thres, app.figOpts.D3.(num));

                    hold on
                    plt = plot(app.data.(smpl).Surfaces.(Type).(name).Surf, 'DisplayName', name);

                    xlim(LimX)
                    ylim(LimY)
                    if app.figOpts.D3.(num) % 3D plot
                        zlim(LimZ)
                    end

                    plt.FaceColor = face_color;
                    plt.FaceAlpha = face_alpha;
                    plt.EdgeColor = 'k';
                    plt.EdgeAlpha = '0.1';
                    hold off

                    % Add the surface to the options menu
                    n = find(strcmp(app.DataN.Items, smpl));

                    NewSurf = uimenu(app.PlotMenu.(num).Surfaces.(Type).(['SMPL' (num2str(n))]));
                    NewSurf.MenuSelectedFcn = @(~,~) Plt_Helper.frontend_surf(app, num, smpl, name, true, Type);
                    NewSurf.Text = name;
                end

            end % end if edit
        end % end Edit Surface

        function del_surf(app, smpls, Type, name)
            % DEL_SURF Deletes a surface in all samples with a given name
            % and a given type
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smpls - Samples in which to remove a surface
            %   - Type - Type of surface. This determines place from which
            %       surface will be removed. Example: UDS - User Defined
            %       Surface
            %   - name - Name of surface to be removed.
            %
            % Modifies:
            %   - app - More specifically app.data.(smpl).Surfaces.(Type)
            %       has a field name removed.
            
            if ~iscell(smpls)
                smpls = {smpls};
            end
            for smpl_idx=1:numel(smpls)
                smpl = smpls{smpl_idx};
                % Delete the surface
                plt = app.data.(smpl).Surfaces.(Type).(name).Plot;
                if ~isnumeric(app.data.(smpl).Surfaces.(Type).(name).Plot) && isvalid(app.data.(smpl).Surfaces.(Type).(name).Plot)
                    delete(plt)
                end
                app.data.(smpl).Surfaces.(Type) = rmfield(app.data.(smpl).Surfaces.(Type), name);
            end
        end

        function add_sub_surf(app, num, smpls, Type, parent1_name, parent2_name, pop_name, face_color, face_alpha, outside)
            % ADD_SUB_SURF Creates a new surface which is an overlap
            % between 2 given surfaces. It will take points of second
            % parent, and gate them of first parent.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - num - Figure identifier on which surface will be plotted,
            %       and from which data will be pulled.
            %   - smpls - Samples on which to define/modify surface
            %   - Type - Type of surface. This determines place to which
            %       new surface will be saved. Example: UDS - User Defined
            %       Surface
            %   - parent1_name - Name of first parent surface.
            %   - parent2_name - Name of second parent surface.
            %   - pop_name - Name of new surface which will be created from
            %       overlap of name and parent_surf_name.
            %   - face_color - Color of the face of the new surface.
            %   - face_alpha - Transparency of the face of the new surface.
            %   - outside - Whether to choose points which are in parent_2
            %       inside of parent_1 (false), or outside of parent_1
            %       (true).
            %
            % Modifies:
            %   - app - Adds in a surface struct, under a given sample,
            %       type and name. Fills it in with information needed to
            %       reference, and reconstruct the surface (except for 
            %       region/hole threshold information).
            
            %% Pull the currently plotted data
            vrnms = app.figOpts.vrnms.(num);

            % Get the currently plotted surface
            if ~iscell(smpls)
                smpls = {smpls};
            end
            for smpl_idx=1:numel(smpls)
                smpl = smpls{smpl_idx};
                Vol = app.data.(smpl).Surfaces.(Type).(parent1_name).Surf;

                % Define the new populations name
                pop_name = Helper.valid_gate(pop_name); % Make it compatible with other gates.
                if iscell(pop_name)
                    pop_name = pop_name{1};
                end

                % Get the Parent surface
                VolParent = app.data.(smpl).Surfaces.(Type).(parent2_name).Surf;
                IND = inShape(Vol, VolParent.Points);
                if outside
                    IND = ~IND;
                end
                VolParent.Points = VolParent.Points(IND, :);
                % Plot the gatted child surface
                hold on
                plt = plot(VolParent, 'DisplayName', pop_name);
                hold off
                plt.FaceColor = face_color;
                plt.FaceAlpha = face_alpha;
                plt.EdgeColor = 'k';
                plt.EdgeAlpha = '0.1';

                % add the surface to the list of surfaces
                app.data.(smpl).Surfaces.(Type).(pop_name).Surf = VolParent;
                app.data.(smpl).Surfaces.(Type).(pop_name).Plot = plt;
                app.data.(smpl).Surfaces.(Type).(pop_name).Phenotype = app.figOpts.phenos.(num);
                app.data.(smpl).Surfaces.(Type).(pop_name).Axes = ...
                    {vrnms{app.figOpts.xaxIF.(num).Value}, ...
                     vrnms{app.figOpts.yaxIF.(num).Value}, ...
                     vrnms{app.figOpts.zaxIF.(num).Value}};

                % Put the new surface in the menu bar
                n = find(strcmp(app.DataN.Items, smpl));

                NewSurf = uimenu(app.PlotMenu.(num).Surfaces.(Type).(['SMPL' (num2str(n))]));
                NewSurf.MenuSelectedFcn = @(~,~) Plt_Helper.frontend_surf(app, num, smpl, pop_name, true, Type);
                NewSurf.Text = pop_name;
            end
        end

        function Vol = surf_create(app, smpl, phn, dat, Xax, Yax, Zax, alpha, Type, name, RegionThreshold, HoleThreshold, is3D)
            % SURF_CREATE Creates a surface, and saves it to a given
            % sample.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smpl - Sample to on which to define surface, and to which
            %       save it
            %   - phn - Phenotypes on which to define surface. (Should
            %       agree with dat given)
            %   - dat - table - Currently plotted data. Output of Plt_Helper.UpdateFigMenu
            %   - Xax - Name of first axis. Can be in full version (i.e.
            %       one shown to the user).
            %   - Yax - Name of second axis. Can be in full version (i.e.
            %       one shown to the user).
            %   - Zax - Name of third axis. Can be in full version (i.e.
            %       one shown to the user).
            %   - alpha - Transparency of the surface. Value between 0-1.
            %   - Type - Type of surface. This determines place to which
            %       surface will be saved. Example: UDS - User Defined
            %       Surface
            %   - name - Name of this particular surface.
            %   - RegionThreshold - Maximum threshold of supressing
            %       regions in the surface. For more details, better go to
            %       alphashape in matlab docs.
            %   - HoleThreshold - Maximum threshold of filling in holes
            %       in the surface. For more details, better go to
            %       alphashape in matlab docs.
            %
            % Output:
            %   - Vol - alphashape - Shape of surface created which can be
            %       plotted on the 
            %
            % Modifies:
            %   - app - Adds in a surface struct, under a given sample,
            %       type and name. Fills it in with information needed to
            %       reference, and reconstruct the surface (except for 
            %       region/hole threshold information).
            if iscell(smpl)
                if numel(smpl) ~= 1
                    error('Func_SurfCreate accepts only single sample. Please loop through them, and call it for each one');
                else
                    smpl = smpl{1};
                end
            end
            Xax = Helper.valid(app, Xax, smpl);
            if iscell(Xax)
                Xax = Xax{1};
            end
            Yax = Helper.valid(app, Yax, smpl);
            if iscell(Yax)
                Yax = Yax{1};
            end
            Zax = Helper.valid(app, Zax, smpl);
            if iscell(Zax)
                Zax = Zax{1};
            end
            %% Pull the X, Y, Z data
            X = dat.(Xax);
            Y = dat.(Yax);
            if isempty(Zax) || is3D~=1
                Z = 0. * X;
            elseif ~isempty(Zax)
                Z = dat.(Zax);
            end
            %% Create the surface
            if numel(unique(Z))~=1 % 3D plot
                Vol = alphaShape(X, Y, Z, alpha);
            else % 2D plot
                Vol = alphaShape(X, Y, alpha);
            end
            Vol.RegionThreshold = RegionThreshold;
            Vol.HoleThreshold = HoleThreshold;

            % If no surfaces have been define on this data
            if ~isfield(app.data.(smpl).Surfaces, Type)
                app.data.(smpl).Surfaces.(Type) = struct;
            end

            % add the surface to the list of surfaces
            app.data.(smpl).Surfaces.(Type).(name).Surf = Vol;
            app.data.(smpl).Surfaces.(Type).(name).Plot = [];
            app.data.(smpl).Surfaces.(Type).(name).PlotNum = [];
            app.data.(smpl).Surfaces.(Type).(name).Phenotype = phn;
            app.data.(smpl).Surfaces.(Type).(name).Axes = {Xax, Yax, Zax};

            %% Add lighting
%             shading interp
%             h = light;
%             lightangle(h, 10, 35);
%             h.FaceLighting = 'gouraud';
%             plt.AmbientStrength = 0.3;
%             plt.DiffuseStrength = 0.8;
%             plt.SpecularStrength = 0.9;
%             plt.SpecularExponent = 25;
%             plt.BackFaceLighting = 'unlit';

        end

        function rightclick(~,~)
           'Congrats on clicking the do nothing option';
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Export Plots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function export_plot_data(app, num)
            % EXPORT_PLOT_DATA It is front-end to function which allows user
            % to export the data which is currently plotted in a given
            % figure
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            %   - num - char, string - Identifier of figure from which data
            %       will be saved
            
            % If the plot is a correlation heatmap
            if strcmp(num, 'Corr_Heatmap')
                ax = gca;
                if ~iscell(ax.Title)
                    ax.Title = {ax.Title};
                end
                dnm = strrep(strjoin(ax.Title, '_'), ' ', '');
% % %                 HDR = ax.XDisplayLabels;
% % %                 datMAT = ax.ColorData;
            elseif strcmp(num, 'Line_Point_Plot')
                ax = gca;
                if ~iscell(ax.Title.String)
                    ax.Title.String = {ax.Title.String};
                end
                dnm = strrep(strjoin(ax.Title.String, '_'), ' ', '');
            elseif strcmp(num, 'Distance_Box_Plot')
                ax = gca;
                if ~iscell(ax.Title.String)
                    ax.Title.String = {ax.Title.String};
                end
                dnm = strrep(strjoin(ax.Title.String, '_'), ' ', '');
            else
                PltDat = get(gca,'Children');
                for Child_i=1:numel(PltDat)
                    dnm = PltDat(Child_i).DisplayName;
                    dnm = strrep(dnm, '/', '_');
                    [file, path] = uiputfile('*.csv', ['Save file name for ' dnm], ...
                        dnm);
                    if isempty(file)
                        return;
                    end
                    Plt_Helper.export_plot_data_backend(app, num, file, path, Child_i);
                end
                return;
            end
            [file, path] = uiputfile('*.csv', ['Save file name for ' dnm], ...
                    dnm);
            if isempty(file)
                return;
            end
            if file==0
                return;
            end
            Plt_Helper.export_plot_data_backend(app, num, file, path, 0);
        end

        function export_plot_data_backend(app, num, file, path, Child_i)
            % EXPORT_PLOT_DATA_BACKEND It is back-end to function which allows
            % user to export the data which is currently plotted in a given
            % figure
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            %   - num - char, string - Identifier of figure from which data
            %       will be saved
            %   - file - char, string - Name of file to which data should
            %       be saved
            %   - path - char, string - path of where this file should be
            %       located
            %   - Child_i - Number, selecting which child is currently being
            %       saved fromt he plot
            
            curr_path = pwd;
            cd(path)

            %% If the plot is a correlation heatmap
            switch num
                case 'Corr_Heatmap'
                    ax = gca;
                    HDRcol = [{'.'}; ax.XDisplayLabels];
                    HDRrow = ax.YDisplayLabels;
    %                 datMAT = ax.ColorData;
                    datMAT = ax.ColorDisplayData;
                    fid = fopen(file,'w');
                    headerformat = '%s\n';
                    dataformat = '%f\n';
                    % If there were more than one columns...
                    for numcolumns=1:(size(datMAT, 2))
                        headerformat = ['%s, ' headerformat];
                        dataformat = ['%f, ' dataformat];
                    end
                    fprintf(fid,headerformat,HDRcol{:});
                    fclose(fid);

                    % COnvert to a cell and add row names
                    datTBL = [HDRrow,num2cell(datMAT)];
                    writecell(datTBL, file,'Delimiter',',','WriteMode','append')

                    return
                
                case 'Line_Point_Plot'
                    %% If the plot is a line plot
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
                            if numel(Y)==1
                                datMAT = [datMAT, [X; Y]];
                                HDR = [HDR, {YLBL}];
                                inclX = 3;
                            elseif size(datMAT, 2) < numel(Y)
                                % deal with non-unique X-Y pairs
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

                    if numel(Y)==1
                        HDR = ['.'; HDR(2:end)'];
                        YLBL = ax.YLabel.String;
                        YLBL = strrep(YLBL, ',', '');
                        HDRrow = {XLBL; YLBL};
                        datTBL = [HDRrow, num2cell(datMAT)];
                        datMAT = datMAT';
                    end

                    fid = fopen(file,'w');
                    headerformat = '%s\n';
                    dataformat = '%f\n';
                    % If there were more than one columns...
    % % %                 for numcolumns=1:(size(datMAT, 1)-1)
                    for numcolumns=1:(numel(HDR)-1)
                        headerformat = ['%s, ' headerformat];
                        dataformat = ['%f, ' dataformat];
                    end
                    fprintf(fid,headerformat,HDR{:});
                    fclose(fid);

                    if numel(Y)==1
                        writecell(datTBL, file,'Delimiter',',','WriteMode','append')
                    else
                        dlmwrite(file, datMAT', '-append', 'precision', 5)
                    end
                case 'Distance_Box_Plot'
                    %% If the plot is a distance box plot
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

                    % Reshape the data table
                    datMATNew = [datMAT(3:2:end, 1), datMAT(4:2:end, 1)];
                    [~, INDs] = sort(datMAT(3:2:end, 1));
                    datMATNew = datMATNew(INDs, :);

                    NSamples = floor(size(datMATNew, 1)/numel(unique(datMAT(1:2:end, 1))))+1;
                    NCells = numel(unique(datMAT(1:2:end, 1)))-1;
                    try
                        datMATNew = reshape(datMATNew(:,2), [NSamples, NCells]);
                    catch
                        NSamples = NSamples+1;
                        datMATNew = reshape(datMATNew(:,2), [NSamples, NCells]);
                    end

                    % Find the cell names
                    HDRNew = HDR(4:2:end);
                    HDRNew = HDRNew(INDs);
                    HDRNew = HDRNew(1:NSamples:end);

                    if numel(HDRNew) == 1
                        % Don't rename if one .csv is loaded in.
                        same_start_str = 1;
                    else
                        same_start_str = HDRNew{1};
                        while ~all(startsWith(HDRNew, same_start_str, 'IgnoreCase', true))
                            if isempty(same_start_str)
                                break;
                            else
                                same_start_str = same_start_str(1:end-1);
                            end
                        end
                        same_start_str = numel(same_start_str) + 1;  % MatLab has inclusive indexing.
                    end
                    cellnames = HDRNew;
                    for cellnmi = 1:numel(cellnames)
                        cellnames{cellnmi} = cellnames{cellnmi}(same_start_str:end);
                    end
                    cellnames = [{''}, cellnames];

                    % Find the sample names
                    HDRNew = HDR(4:2:end);
                    HDRNew = HDRNew(INDs);
                    HDRNew = HDRNew(1:(NSamples));
                    if numel(HDRNew) == 1
                        % Don't rename if one .csv is loaded in.
                        same_end_str = 1;
                    else
                        same_end_str = HDRNew{1};
                        while ~all(endsWith(HDRNew, same_end_str, 'IgnoreCase', true))
                            if isempty(same_end_str)
                                break;
                            else
                                same_end_str = same_end_str(2:end);
                            end
                        end
                        same_end_str = numel(same_end_str) + 1;  % MatLab has inclusive indexing.
                    end
                    smplnames = HDRNew;
                    for cellnmi = 1:numel(smplnames)
                        smplnames{cellnmi} = smplnames{cellnmi}(1:same_end_str);
                    end
                    smplnames = strrep(smplnames, ' ', '');
                    smplnames = [{''}, smplnames];

                    fid = fopen(file,'w');
                    headerformat = '%s\n';
                    dataformat = '%f\n';
                    % If there were more than one columns...
                    for numcolumns=1:(size(datMATNew, 2)+1)
                        headerformat = ['%s, ' headerformat];
                        dataformat = ['%f, ' dataformat];
                    end
                    fprintf(fid,headerformat,cellnames{:});
                    fclose(fid);

                    fid = fopen(file,'a');
                    headerformat = '%s\n';
                    dataformat = '%f\n';
                    % If there were more than one columns...
                    for numcolumns=1:(size(datMATNew, 1))
                        headerformat = ['%s\n ' headerformat];
                        dataformat = ['%f\n ' dataformat];
                    end
                    fprintf(fid,headerformat,smplnames{:});
                    fclose(fid);

                    dlmwrite(file, datMATNew, '-append',  'precision',5,'coffset', 1)
                    
                otherwise
                    %% For all other cases 
                    % Pull the data in the figure
                    PltDat = get(gca,'Children');
                    ax = gca;
                    if strcmp(app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, 'None') % If there is no color axis
                        if strcmp(app.figOpts.yaxIF.(num).String{app.figOpts.yaxIF.(num).Value}, 'HistogramNum') % If user selected histogram
                            X = PltDat(Child_i).XData;
                            Y = PltDat(Child_i).YData;
                            YLBL = ax.YLabel.String;
                            XLBL = ax.XLabel.String;
                            HDR = {XLBL, YLBL};
                            datMAT = [X; Y];
                        else
                            X = PltDat(Child_i).XData;
                            Y = PltDat(Child_i).YData;
                            YLBL = ax.YLabel.String;
                            XLBL = ax.XLabel.String;
                            if app.figOpts.D3.(num)
                                Z = PltDat(Child_i).ZData;
                                ZLBL = ax.ZLabel.String;
                                HDR = {XLBL, YLBL, ZLBL};
                                datMAT = [X; Y; Z];
                            else
                                HDR = {XLBL, YLBL};
                                datMAT = [X; Y];
                            end
                        end
                    elseif strcmp(app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, 'Number of cells / Neighborhood')
                        %%% This ploting section isn't working now
                        X = PltDat(Child_i).XData;
                        Y = PltDat(Child_i).YData;
                        Z = PltDat(Child_i).ZData;
                        C = PltDat(Child_i).CData;
                        CLBL = 'Cells / mm^3 / Neighborhood';
                        HDR = {'X', 'Y', 'Z', CLBL};
                        datMAT = [X; Y; Z; C];
                    else % If there is some color axis
                        X = PltDat(Child_i).XData;
                        Y = PltDat(Child_i).YData;
                        C = PltDat(Child_i).CData'; % why does MATLAB switch formats...

                        CLBL = app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value};
                        YLBL = ax.YLabel.String;
                        XLBL = ax.XLabel.String;
                        if app.figOpts.D3.(num)
                            Z = PltDat(Child_i).ZData;
                            ZLBL = ax.ZLabel.String;
                            HDR = {XLBL, YLBL, ZLBL, CLBL};
                            datMAT = [X; Y; Z; C];
                        else
                            HDR = {XLBL, YLBL, CLBL};
                            datMAT = [X; Y; C];
                        end
                    end

                    fid = fopen(file,'w');
                    headerformat = '%s\n';
                    dataformat = '%f\n';

                    % If there were more than one columns...
                    for numcolumns=1:(size(datMAT, 1)-1)
                        headerformat = ['%s, ' headerformat];
                        dataformat = ['%f, ' dataformat];
                    end

                    fprintf(fid,headerformat,HDR{:});
                    fclose(fid);
                    dlmwrite(file,datMAT', '-append', 'precision', 5)                
            end
            cd(curr_path);
        end
        
        function export_figure(app, num)
            % EXPORT_FIGURE It is front-end of function that allows user to
            % save current figure as an image.
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            %   - num - char, string - Identifier of figure which will be
            %       saved.
            
            % Pull the selected Phenotype
            phnm = app.figOpts.phenos.(num);
            phnm = strjoin(phnm, "_");

            % Pull the selected sample name
            smpl = app.figOpts.smpls.(num);
            smpl = strjoin(smpl, "_");

            [file, path] = uiputfile({'*.jpg'; '*.tif'; '*.png'; '*.fig'},['Save figure ' smpl '_' phnm], ...
                [smpl '_' phnm]);
            if isempty(file)
                return;
            end
            
            Plt_Helper.export_figure_backend(app, num, file, path);
        end

        function export_figure_backend(app, num, file, path)
            % EXPORT_FIGURE_BACKEND It is back-end of function that allows
            % user to save current figure as an image.
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            %   - num - char, string - Identifier of figure which will be
            %       saved.
            %   - file - char, string - Name of file to which figure will
            %       be saved.
            %   - path - char, string - Path to that file.
            curr_path = pwd;
            cd(path);

            figSV = figure;
            % add an export data option
            % Create Menu bar
            ExportMenu = uimenu(figSV);
            ExportMenu.Text = 'Export';
            % Create an export data option
            ExportPltDat = uimenu(ExportMenu);
            ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
            ExportPltDat.Text = 'Export Plot Data to .csv';
            figSV.Color = app.figOpts.bgclr.(num);
            figSV.InvertHardcopy = 'off';

            Plotting.func_plot(app, num)

            saveas(figSV,file)
%             close(figSV)

            cd(curr_path);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Edit Plots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plot_converter(app, num)
            % PLOT_CONVERTER It is a function which allows the user
            % to convert a plot from a line plot to a heatmap
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            %   - num - char, string - Identifier of figure 
            
            % If the plot is a correlation heatmap
            switch num
                case 'Corr_Heatmap'
% % %                 ax = gca;
% % %                 dnm = strrep(strjoin(ax.Title, '_'), ' ', '');
% % %                 HDR = ax.XDisplayLabels;
% % %                 datMAT = ax.ColorData;
                case 'Line_Point_Plot'
                    ax = gca;
                    if ~iscell(ax.Title.String)
                        ax.Title.String = {ax.Title.String};
                    end
                    PltDat = get(gca,'Children');
                    for i=1:numel(PltDat)
                        if strcmp(PltDat(i).Tag, 'boxplot')
                            PltDat(i) =[];
                        end
                    end
                    for i=1:numel(PltDat)
                        X = PltDat(i).XData;
                        Y = PltDat(i).YData;
                        YLBL = PltDat(i).DisplayName;
                        XLBL = ax.XLabel.String;
                        if i==1
                            HDR = {XLBL, YLBL};
                            datMAT = [X; Y];
                        else
                            HDR = [HDR, {YLBL}];
                            datMAT = [datMAT; Y];
                        end
                    end
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
                    
                    hm = heatmap(fliplr(datMAT(2:end,:)'));
                    hm.Colormap = app.GUIOPTS.redbluecmap;
                    hm.XLabel = 'Region';
                    hm.XDisplayLabels = fliplr(HDR(2:end));
                    hm.YDisplayLabels = ax.XTickLabels;
                    ttl = ax.Title.String;
                    if iscell(ttl)
                        ttl = ttl{1};
                    end
                    title(ttl)
                    axs = struct(gca); %ignore warning that this should be avoided
                    cb = axs.Colorbar;
                    % Set the color axes in a logical way
                    switch ttl
                        case 'Fold Change'
                            limits = [min(min(datMAT(2:end,:))), max(max(datMAT(2:end,:)))]-1;
                            [~,MAX] = max(abs(limits));
                            [~,MIN] = min(abs(limits));
                            limits(MIN) = -1*limits(MAX); % Make the color limits symetrical
                            hm.ColorLimits = limits+1;
                            % put a title on the colorbar
                            cb.Label.String = 'Fold Change';
                        case 'Difference from mean'
                            limits = [min(min(datMAT(2:end,:))), max(max(datMAT(2:end,:)))]-1;
                            [~,MAX] = max(abs(limits));
                            [~,MIN] = min(abs(limits));
                            limits(MIN) = -1*limits(MAX); % Make the color limits symetrical
                            hm.ColorLimits = limits;
                            cb.Label.String = 'Difference from mean';
                    end
                
                otherwise
% %                 PltDat = get(gca,'Children');
% %                 for Child_i=1:numel(PltDat)
% %                     dnm = PltDat(Child_i).DisplayName;
% %                     dnm = strrep(dnm, '/', '_');
% %                     [file, path] = uiputfile('*.csv', ['Save file name for ' dnm], ...
% %                         dnm);
% %                     if isempty(file)
% %                         return;
% %                     end
% %                     Plt_Helper.export_plot_data_backend(app, num, file, path, Child_i);
% %                 end
% %                 return;
            end % end of conditional statement
        end % End of plot convert
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display and other functionality
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function highight_reg_func(app, num)
            % Pull the available Model Names
            modelnames = fieldnames(app.net);
            
            % pull what is currently on the color-axis
            caxIF = app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value};
            caxIF = Helper.valid_var(caxIF);
            if strcmp(caxIF, 'None')
                % 'There is no color axis to highlight'
            elseif contains(caxIF, modelnames)
                % 'A model was chosen, highlight cluster or region'
                ind = strcmp(caxIF, modelnames);
                % Get regions
                if ismember('RegNames', fieldnames(app.net.(modelnames{ind})))
                   list = app.net.(modelnames{ind}).RegNames; 
                else
                    list = 1:app.net.(modelnames{ind}).NR;
                    list = num2str(list');
                end
                output = listdlg('PromptString', 'Select Regions to Highlight',...
                        'ListString', list, ...
                        'OKString', 'Highlight');
                
                %% Now pull stats about the selected regions
                
                figh = get(gca,'children');
                Cdat = get(figh,'CData');
                if iscell(Cdat)
                    Sdat = cell(size(Cdat));
                    for cell_i = 1:size(Cdat, 1)
                        Sdat{cell_i}= ones(numel(Cdat{cell_i}), 1);
                    end
                else
                    Sdat = ones(numel(Cdat), 1);
                end
                
                if numel(output) > 1
                    if iscell(Cdat)
                        IND = cell(size(Cdat));
                        for cell_i = 1:size(Cdat, 1)
                            IND{cell_i} = sum(Cdat{cell_i}==output, 2)~=0;
                        end
                    else
                        IND = sum(Cdat==output, 2)~=0;
                    end
                else
                    if iscell(Cdat)
                        IND = cell(size(Cdat));
                        for cell_i = 1:size(Cdat, 1)
                            IND{cell_i} = Cdat{cell_i}==output;
                        end
                    else
                        IND = Cdat==output;
                    end
                end
                if iscell(Cdat)
                    for cell_i = 1:size(Cdat, 1)
                        Sdat{cell_i}(IND{cell_i}) = 20;
                        set(figh(cell_i),'SizeData', Sdat{cell_i});
                    end
                    
                else
                    Sdat(IND) = 20;
                    set(figh,'SizeData', Sdat);
                end
                
                
                
            else
                % 'Somme other option was chosen, do something with that'
            end
            
           
        end
        
        function highight_Neigh_func(app, num)
            modelnames = fieldnames(app.net);
            
            % pull what is currently on the color-axis
            caxIF = app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value};
            caxIF = Helper.valid_var(caxIF);
            
            % Pull the current sample
            smpls = app.figOpts.smpls.(num);
            
            % Check to make sure there is a color axis to highlight
            if strcmp(caxIF, 'None')
                return
            elseif ~contains(caxIF, modelnames)
                return                
            end
            
            %% 
            
            % 'A model was chosen, highlight cluster or region'
            ind = strcmp(caxIF, modelnames);
            % Get regions
            if ismember('RegNames', fieldnames(app.net.(modelnames{ind})))
               list = app.net.(modelnames{ind}).RegNames; 
            else
                list = 1:app.net.(modelnames{ind}).NR;
                list = num2str(list');
            end
            Ri = listdlg('PromptString', 'Select Region to Highlight',...
                    'ListString', list, ...
                    'OKString', 'Highlight', 'SelectionMode', 'single');
                
           if isempty(Ri)
               return
           end
           
           switch app.net.(modelnames{ind}).userdata.DataType
               case 'Individual Cells'
                   type = 'AllCells';
               case 'Raster Scanned Neighborhoods'
                   type = 'MFIRSN';
               case 'Cell Centered Neighborhoods'
                   type = 'MFICCN';
           end
%%%%%%%%%%%%%%%%%%%%  
           for smpli = 1:numel(smpls)
                    
                radius = app.rwindowRSN;
                NNPlts = 27;
                % Make the figure
                fig = figure;
                fig.Color = 'w';
                cla

                subplot(1,2,1)
                cla
                hold on
                box off
                axis equal

                subplot(1,2,2)
                cla
                hold on
                box off
                axis equal

                % Select model
                mdlVAR = [Constants.other_tag caxIF];
                
                % Pull indeces for random clusters from type Ri
                Neigh = app.data.(smpls{smpli}).(type)( ...
                    app.data.(smpls{smpli}).(type).(mdlVAR)==Ri, ...
                    {'X', 'Y', 'Z', 'NCells'});

                NCellsU = unique(Neigh.NCells);
                NCellsi = 1:round(numel(NCellsU)/NNPlts):numel(NCellsU);
                NCellsi = NCellsU(NCellsi);

                %%
                Npos = zeros(numel(NCellsi),3);
                cellsI = app.data.(smpls{smpli}).AllCells{:, {'X', 'Y', 'Z'}};
                INDKeep = false(size(cellsI,1), 1);
                PosSub = zeros(size(cellsI,1), 2); 
                
                
                subplot(1,2,2)
                plot(cellsI(:,1), cellsI(:,2), '.', 'color', [0.8,0.8,0.8], 'MarkerSize', 1, 'DisplayName', 'All Cells')

                xi = 0;
                yi = 0;
                for i = 1:numel(NCellsi)
                    NCellsIND = find(Neigh.NCells == NCellsi(i));
                    Npos(i, :) = Neigh{NCellsIND(1),{'X', 'Y', 'Z'}};
                                        
                    % Find cells within neighborhood
                    dist = sqrt(sum((Npos(i,:)-cellsI).^2, 2));
                    IND = dist <= radius;
                    INDKeep(IND) = true;
                    
                    if i==1
                        cellsj = app.data.(smpls{smpli}).AllCells(IND, :);
                        cellsj.X = cellsj.X + xi+xi*radius*2.2-Npos(i,1);
                        cellsj.Y = cellsj.Y + yi+yi*radius*2.2-Npos(i,2);
                    else
                        cellsjtmp = app.data.(smpls{smpli}).AllCells(IND, :);
                        cellsjtmp.X = cellsjtmp.X + xi+xi*radius*2.2-Npos(i,1);
                        cellsjtmp.Y = cellsjtmp.Y + yi+yi*radius*2.2-Npos(i,2);
                        cellsj = [cellsj; cellsjtmp];
                    end
                                        
                    PosSub(IND, 1) = xi+xi*radius*2.2-Npos(i,1);
                    PosSub(IND, 2) = yi+yi*radius*2.2-Npos(i,2);
                    
                    x = xi+xi*radius*2.2;
                    y = yi+yi*radius*2.2;
                    subplot(1,2,1)
% %                     plot(x,y, 'xk')
%                     plot(cellsI(IND,1)-Npos(i,1)+x, cellsI(IND,2)-Npos(i,2)+y, '.')

                    subplot(1,2,2)
                    plot(cellsI(IND,1), cellsI(IND,2), '.')

                    if rem(i,round(sqrt(numel(NCellsi))))==0
                        yi = yi+1;
                        xi = 0;
                    else
                        xi = xi+1;
                    end
                end
                
               %%
                subplot(1,2,1)
                plot(cellsj.X, cellsj.Y, '.', 'color', [0.8,0.8,0.8], 'MarkerSize', 1, 'DisplayName', 'All Cells')
% % %                 plot((cellsI(INDKeep,1)+PosSub(INDKeep, 1)), (cellsI(INDKeep,2)+PosSub(INDKeep, 2)), '.')

                phenos = app.figOpts.phenos.(num);
                for phn_i = 1:numel(phenos)
                    
                    gatei = Helper.gate_full2tag(app, phenos{phn_i}, smpls{smpli});
% %                     INDCelli = cellsj.(gatei).*INDKeep==1;
                    INDCelli = cellsj.(gatei)==1;
                    
                    if isempty(cellsj{INDCelli,{'X'}})
                        plot( ...
                            [0],...
                            [0], ...
                            '.', 'DisplayName', phenos{phn_i})
                    else
                        plot( ...
                            (cellsj{INDCelli,{'X'}}),...
                            (cellsj{INDCelli,{'Y'}}), ...
                            '.', 'DisplayName', phenos{phn_i})
                    end

                end
                axis off
                ax = gca;
                ax.Title.String = ['Sample: ' smpls{smpli} ' Neighborhoods from: ' list(Ri)];
% % %                subplot(1,2,1)
% % %                plot(cellsI(INDKeep,1), cellsI(INDKeep,2), '.')
           end
%%%%%%%%%%%%%%%%%%%%%            
        end
        
        function txt = PlotClickFunc(obj, event_obj, app, num)
            
            datcurs = datacursormode(gcf);
            if strcmp(datcurs.Enable, 'off')
                txt = 'Data Tips Disabled';
                datcurs.Enable = 'on';
                return
            end
            
% %             datcurs.Enable = 'on';

            % pull the fieldnames
            fnms = fieldnames(event_obj.Target);

            % Pull the plotted phenotype names
            phenos = app.figOpts.phenos.(num);

            % Pull the sample names
            smpl = app.figOpts.smpls.(num);

            % Pull all available phenotypes for the plotted samples
            [~, Phenotypes] = Helper.get_gates(app, smpl);

            % get the index of the selected point
            IND = get(event_obj, 'DataIndex');

            % get the position of the selected point
            switch numel(event_obj.Position)
                case 2
                    % case for if the plot is 2D
                    x = event_obj.Position(1);
                    y = event_obj.Position(2);
                    z = 'none';
                case 3
                    %case for if the plot is 3 dimensional
                    x = event_obj.Position(1);
                    y = event_obj.Position(2);
                    z = event_obj.Position(3);
            end

            % if there is a color axis pull the value
            if ~strcmp(app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, 'None')
                % Pull the name of the color axis
                cname = strrep(app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, '_', ' ');
                % Pull the value of the color axis
                colors = event_obj.Target.CData;
                clr = colors(IND);
            else
                % Pull the name of the color axis
                cname = 'Color';
                clr = 'none';   
            end

            % Build the default display box
                txt = {strrep(event_obj.Target.DisplayName, '_', ' '), ...
                       ['X: ' num2str(x)], ...
                       ['Y: ' num2str(y)], ...
                       ['Z: ' num2str(z)], ...
                       [cname ': ' num2str(clr)], ...
                       };

            % this will only be true if the current color axis is a model
            %  contains(cname, strrep(fieldnames(app.net), '_', ' '))

            % only run if the plotted objects are neighborhoods
            if contains(cname, strrep(fieldnames(app.net), '_', ' '))

                % make the currently selected points larger
                SizeData = ones(size(colors)) * 2;
                SizeData(colors==clr) = 20;
                event_obj.Target.SizeData = SizeData;
                % find the number of points in the selected group
                NNeighbor = sum(colors==clr);
                % find the percentage of total points in the current group
                percent = round(NNeighbor/size(colors,1)*1000)/10;

                try
                    cellTypeOfRegions = Plt_Helper.modelstats_func(app, num, smpl, phenos);
                    cellType = cellTypeOfRegions(clr);
                catch
                    'Something went wrong';
                    txt = {['* ', 'Region : ', num2str(clr)],...
                        ['* ', 'N : ', num2str(NNeighbor)],...
                        ['* ', num2str(percent), '% of currently plotted points are R', num2str(clr)],...
                        };
                    return
                end

                txt = {['* ', 'Region : ', num2str(clr)],...
                       ['* ', 'N : ', num2str(NNeighbor)],...
                       ['* ', num2str(percent), '% of currently plotted points are R', num2str(clr)],...
                       ['* ', 'Composed of Mostly ', char(cellType), ' cells']};
            end
        end %end of plot click function
        
        function cellNames = modelstats_func(app, num, smplnms, Phenotypes)

            cname = strrep(app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, '_', ' ');
            % just in case this is run with color not from model
            if ~contains(cname, strrep(fieldnames(app.net), '_', ' '))
                return
            end

            % deal with edge cases
            if strcmp(Phenotypes, {'All Cells'}) ...
                    || strcmp(Phenotypes, {'Density/MFI RSN'}) ...
                    || all(startsWith(Phenotypes, 'Nei_'))
                % Pull all available phenotypes for the plotted samples
                [~, Phenotypes] = Helper.get_gates(app, smplnms);
                Phenotypes = Phenotypes';
            end

            % Convert the color axis name to a valid model name
            DataType = strrep(cname,' ','_');

            % Convert the model name to a valid channel name
            MFIList = {strcat(Constants.other_tag, DataType)};

            % make sure the selected model is in all plotted samples
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

            % Initialize data type
            if strcmp(app.net.(DataType).userdata.DataType, 'Raster Scanned Neighborhoods')
                type = 'MFIRSN';
            elseif strcmp(app.net.(DataType).userdata.DataType, 'Cell Centered Neighborhoods')
                type = 'MFICCN';
            elseif strcmp(app.net.(DataType).userdata.DataType, 'Individual Cells')
                type = 'AllCells';
            end

            CellNms = [Phenotypes; MFIList];

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
                app.net.(DataType).userdata.DataType, ...
                cellprep, ...
                NormOpt, ...
                rmvzeros, ...
                mfiprep);
            
            if isempty(Dat_Pre)
                return
            end
            
            RowMAP = app.net.(DataType).cmap;
            NGroups = app.net.(DataType).NR;

            % Initialize the matrices
            MeanCellsM = zeros(numel(smplnms),numel(CellNms),NGroups);
            Percentage = zeros(numel(smplnms),NGroups);

            % For each region
            for j = 1:NGroups
                % For each sample
                for i=1:numel(smplnms)
                    % Load in the data
                    DatSub = Dat_Pre(INDDatPre{i}{1}(1):INDDatPre{i}{1}(2),:);
                    % Load correct Row value (ROW is reduduntant and
                    % can be substituted later, for str value)
                    ROW  = app.data.(smplnms{i}).(type).([Constants.other_tag, DataType]);
                    % sort all of the scan data in the same order
                    if i==1
                        if isempty(DatSub(ROW==j, :))
                            MeanCells = DatSub(1, :);
                        else
                            MeanCells = DatSub(ROW==j, :);
                        end
                    else  % Order all data in the same way
                        if isempty(DatSub(ROW==j, :)) % If there are no neighborhoods classified as region j
                            MeanCells = DatSub(1, :);
                        else
                            MeanCells = DatSub(ROW==j, :);
                        end
                    end
                    MeanCells = mean(MeanCells);
                    if isempty(DatSub(ROW==j, :)) % If there are no neighborhoods classified as region j
                        MeanCellsM(i, :, j) = NaN;
                    else
                        MeanCellsM(i, :, j) = MeanCells;
                    end
                    %Calculate the percentage of sample making up that region
                    Percentage(i, j) =  numel(ROW(ROW==j))/numel(ROW(ROW~=0));
                end
            end

            MeanCellsM2 = mean(mean(MeanCellsM, 3, 'omitnan'), 1, 'omitnan'); %Find the average cellularity across the tissues, and regions
            MeanCellsM3 = reshape(mean(MeanCellsM, 1, 'omitnan'), numel(CellNms),NGroups); % Find the average cellularity for each region (across all tissues)

            cellTypes = zeros(1, NGroups);
            cellNames = cell(1, NGroups);
            for j = 1:NGroups
                % Calculate the fold change
                [~, cellTypes(j)] = max(MeanCellsM3(:,j)'./MeanCellsM2);
                cellNames(j) = Phenotypes(cellTypes(j));
            end
        end
        
        function func_rightlick(app, ClickFuncHandle)
            
            if endsWith(ClickFuncHandle.Label, '''on''')
                % 'Switch off'
                ClickFuncHandle.Label = strrep(ClickFuncHandle.Label, 'on', 'off');
            elseif endsWith(ClickFuncHandle.Label, '''off''')
                % 'Switch on'
                ClickFuncHandle.Label = strrep(ClickFuncHandle.Label, 'off', 'on');
            end

        end
        
    end % End of Methods
end % End of Class Definition