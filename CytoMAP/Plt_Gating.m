classdef Plt_Gating

% Plt_Gating defines a suite of functions which plot and gate on
% populations

    methods (Static)
        function Import_Definitions_Func
% % %             import Helper.*;
% % %             import Constants.*;
% % %             import Plt_Helper.*;
            set(groot, 'DefaultTextInterpreter', 'none')
            set(groot, 'DefaultLegendInterpreter', 'none')
        end

        %% Plot Gating Functions
        function new_gate(app, num, GateType, varargin)
            % NEW_GATE Adds a new gate to the CytoMAP, however does not yet
            % save it into samples.
            %
            % Inputs:
            %   - app - Instance of CytoMAP
            %   - num - Identifier of figure from which to pull data on
            %       which gating logical will be calculated.
            %   - GateType - Type of gate, either 'PolyGate' for a polygon
            %       gate, or 'RectGate' for rectangular gate.
            % 
            % Key-word Inputs (optional):
            %   - name - default: '' - Name of the gate. If empty user will
            %       be prompted to input name of the gate.
            %   - position - default: [] - Position of the new gate. If
            %       empty user will have to choose position himself by
            %       clicking on the plot
            %
            % Modifies:
            %   - app - Specifically app.NewGate by adding to it the gate
            %       currently shown.
            
            defaults = struct;
            defaults.name = '';
            defaults.position = [];

            if mod(length(varargin), 2) ~= 0
                error('Uneven number of input arguments! populate_table needs propertyName/propertyValue pairs, after app argument')
            end

            for pair = reshape(varargin, 2, [])
                if isfield(defaults, lower(pair{1}))
                    defaults.(lower(pair{1})) = pair{2}; % Change default axis value to the given value
                else  % Key is not in the defaults
                    error('%s is not a recognized parameter name: ', pair{1})
                end
            end

            name = defaults.name;
            pos = defaults.position;

            %% Initialize some parameters and access the plotted data
            % Pull the selected sample name
            smpls = app.figOpts.smpls.(num);

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

            %% Define the parent population
            Parent = app.figOpts.phenos.(num);
            if ~iscell(Parent)  % Sanity check.
                Parent = {Parent};
            end

            %% Get name and valid name of the gate - Ask user to name gates
            if isempty(name)
                name = inputdlg({'Gate Name:'}, 'New Gate Name', 1, {''});
            end
            if isempty(name)  % If the user hits cancel, end function
                return
            end

            valid_name = Helper.valid_gate(name);
            if iscell(valid_name)
                valid_name = valid_name{1};
            end
            if iscell(name)
                name = name{1};
            end
            X = dat.(vrnms{app.figOpts.xaxIF.(num).Value});
            Y = dat.(vrnms{app.figOpts.yaxIF.(num).Value});

            %% Update currently chosen gate
            if ~isfield(app.NewGate, 'internal__curr') || ~isstruct(app.NewGate.internal__curr)
                app.NewGate.internal__curr = struct;
            end
            app.NewGate.internal__curr.(num) = name;

            %% Pull out an ellipse
            % % %                 imellipse
            %% Gate with a polygon
            if strcmp(GateType, 'PolyGate')
                if isempty(pos)
                    POLY = impoly;
                else
                    POLY = impoly(gca, pos);
                end
                setColor(POLY, 'r')
                POLY.Deletable = false;
                
                app.NewGate.(valid_name).obj = POLY;
                app.NewGate.(valid_name).type = 'Polygon';
                pos = getPosition(app.NewGate.(valid_name).obj);
                Logical = inpolygon(X, Y, pos(:,1), pos(:,2));
            %% Gate with a rectangle
            elseif strcmp(GateType, 'RectGate')
                %% Ask user to make the gate
                if isempty(pos)
                    RECT = imrect;
                else
                    RECT = imrect(gca, pos);
                end
                setColor(RECT, 'r')
                RECT.Deletable = false;

                app.NewGate.(valid_name).obj = RECT;
                app.NewGate.(valid_name).type = 'Rectangle';
                pos = getPosition(app.NewGate.(valid_name).obj);
                minX = pos(1); maxX = pos(1) + pos(3);
                minY = pos(2); maxY = pos(2) + pos(4);
                Logical = inpolygon(X,Y,[minX, maxX],[minY, maxY]);
            end
            
            %%%%%% Try adding a context menu for gates

                % Add context menu to poly object
                cmenu_obj = findobj(app.NewGate.(valid_name).obj,'Type','line','-or','Type','patch');
                cmenu_old = get(cmenu_obj,'uicontextmenu');
                m1 = uimenu(cmenu_old{1}, 'Label', 'Plot Gated Points ''off''',...
                    'Callback', @(m1, ~) Plt_Helper.func_rightlick(app , m1));

            %%%%%%

            % Add common fields to new gate
            app.NewGate.(valid_name).name = name;
            app.NewGate.(valid_name).smpls = smpls;
            app.NewGate.(valid_name).phenos = app.figOpts.phenos.(num);
            app.NewGate.(valid_name).vrnms = vrnms;
            app.NewGate.(valid_name).parent = Parent;
            app.NewGate.(valid_name).GateAxis1 = vrnms{app.figOpts.xaxIF.(num).Value};
            app.NewGate.(valid_name).GateAxis2 = vrnms{app.figOpts.yaxIF.(num).Value};

            update_callback = @(pos) Plt_Gating.update_gate(app, name, num);

            app.NewGate.(valid_name).pos = getPosition(app.NewGate.(valid_name).obj);
            app.NewGate.(valid_name).clr = getColor(app.NewGate.(valid_name).obj);
            app.NewGate.(valid_name).obj.addNewPositionCallback(update_callback);

% %             % Plot gate.
% %             hold on
% %             app.NewGate.(valid_name).plt = plot(X(Logical), ...
% %                 Y(Logical), '.', 'Color', app.NewGate.(valid_name).clr, ...
% %                 'DisplayName', name);
% %             uistack(app.NewGate.(valid_name).plt, 'top');
% %             uistack(app.NewGate.(valid_name).plt, 'down');

            % Add pointers to this gate in samples.
            for smpl_idx=1:numel(smpls)
                smpl = smpls{smpl_idx};
                if ~isfield(app.data.(smpl), 'Gates')
                    app.data.(smpl).Gates = {};
                elseif any(strcmp(app.data.(smpl).Gates, name))
                    % Do nothing. Already pointer to it.
                    continue;
                end
                app.data.(smpl).Gates(end + 1) = {name};
            end
            %% Put the new gate in the menu tree, if it's new.
            already_exists = false;
            for i=1:numel(app.PlotMenu.(num).Gate.UDG.(['SelectSMPL' (num2str(1))]).Children)
                c_i = app.PlotMenu.(num).Gate.UDG.(['SelectSMPL' (num2str(1))]).Children(i);
                if strcmp(name, c_i.Text)
                    already_exists = true;
                    break;
                end
            end
            if ~already_exists
                NewGateMenu = uimenu(app.PlotMenu.(num).Gate.UDG.(['SelectSMPL' (num2str(1))]));
                NewGateMenu.MenuSelectedFcn = @(~,~) Plt_Gating.update_gate(app, name, num);
                NewGateMenu.Text = name;
            end
        end % end new gate function

        function update_gate(app, name, num)
            % UPDATE_GATE It is a function that updates a given gate on a
            % current plot.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - name - Name of the gate to be updated
            %   - num - Identifier of figure from which to pull data on
            %       which gating logical will be calculated.
            %
            % Modifies:
            %   - app - Specifically app.NewGate.(name) by modifying fields
            %       to reflect current state of gate on the plot.
            
            %% Get names
            valid_name = Helper.valid_gate(name);
            if iscell(valid_name)
                valid_name = valid_name{1};
            end
            if iscell(name)
                name = name{1};
            end

            new_smpls = app.figOpts.smpls.(num);
            new_phenos = app.figOpts.phenos.(num);
            old_smpls = app.NewGate.(valid_name).smpls;

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

            %% Update currently chosen gate
            if ~isfield(app.NewGate, 'internal__curr') || ~isstruct(app.NewGate.internal__curr)
                app.NewGate.internal__curr = struct;
            end
            app.NewGate.internal__curr.(num) = name;

            %% Update Gate Params
            % Fields dependend on inputs/axes
            app.NewGate.(valid_name).parent = app.figOpts.phenos.(num);
            app.NewGate.(valid_name).GateAxis1 = vrnms{app.figOpts.xaxIF.(num).Value};
            app.NewGate.(valid_name).GateAxis2 = vrnms{app.figOpts.yaxIF.(num).Value};
            app.NewGate.(valid_name).smpls = new_smpls;
            app.NewGate.(valid_name).phenos = new_phenos;

            if isvalid(app.NewGate.(valid_name).obj)
                pos = getPosition(app.NewGate.(valid_name).obj);
            else
                pos = app.NewGate.(valid_name).pos;
                hold on;
                if strcmp(app.NewGate.(valid_name).type, 'Polygon')
                    app.NewGate.(valid_name).obj = impoly(gca, pos);
                elseif strcmp(app.NewGate.(valid_name).type, 'Rectangle')
                    app.NewGate.(valid_name).obj = imrect(gca, pos);
                end
                setColor(app.NewGate.(valid_name).obj,app.NewGate.(valid_name).clr);
                app.NewGate.(valid_name).obj.addNewPositionCallback(@(pos) Plt_Gating.update_gate(app, name, num));
                hold off;
            end

            % Fields dependend on actual axes
            app.NewGate.(valid_name).pos = pos;
            app.NewGate.(valid_name).clr = getColor(app.NewGate.(valid_name).obj);

            %% Plot gate
            X = dat.(app.NewGate.(valid_name).GateAxis1);
            Y = dat.(app.NewGate.(valid_name).GateAxis2);

            if strcmp(app.NewGate.(valid_name).type, 'Polygon')
                Logical = inpolygon(X, Y,pos(:,1), pos(:,2));
            elseif strcmp(app.NewGate.(valid_name).type, 'Rectangle')
                minX = pos(1); maxX = pos(1) + pos(3);
                minY = pos(2); maxY = pos(2) + pos(4);
                Logical = inpolygon(X,Y,[minX, maxX],[minY, maxY]);
            end

            % Delete previous red points.
            if isfield(app.NewGate.(valid_name), 'plt') % Sanity check.
                delete(app.NewGate.(valid_name).plt)
            end
            
            % Pull poly object context menu setting
            cmenu_obj = findobj(app.NewGate.(valid_name).obj,'Type','line','-or','Type','patch');
            cmenu_old = get(cmenu_obj,'uicontextmenu');
            txt = cmenu_old{1}.Children(1).Text;

            if endsWith(txt, '''on''')
                % 'Plotting is on'
                
                % Plot new points
                hold on
                app.NewGate.(valid_name).plt = plot(X(Logical), ...
                    Y(Logical), '.', 'Color', app.NewGate.(valid_name).clr, ...
                    'DisplayName', name);
                uistack(app.NewGate.(valid_name).plt, 'top');
                uistack(app.NewGate.(valid_name).plt, 'down');
                
            elseif endsWith(txt, '''off''')
                % 'Plotting is off'

            else
                % Context menu element does not exist
                %%%%%% Try adding a context menu for gates

                    % Add context menu to poly object
                    cmenu_obj = findobj(app.NewGate.(valid_name).obj,'Type','line','-or','Type','patch');
                    cmenu_old = get(cmenu_obj,'uicontextmenu');
                    m1 = uimenu(cmenu_old{1}, 'Label', 'Plot Gated Points ''off''',...
                        'Callback', @(m1, ~) Plt_Helper.func_rightlick(app , m1));

                %%%%%%
            end

            %% Update samples cell arrays with gates names.
            % Both remove ones which no longer are in the gate.
            % And add the ones which previously were not in the samples.
            to_remove = setdiff(old_smpls, new_smpls);
            to_add = setdiff(new_smpls, old_smpls);
            for smpl_rmv_idx = 1:numel(to_remove)
                smpl_rmv = to_remove{smpl_rmv_idx};
                % Bunch of sanity checks, which might allow for removing samples, without actually thinking about gates etc.
                if isfield(app.data, smpl_rmv) && isfield(app.data.(smpl_rmv), name) && any(strcmp(app.data.(smpl_rmv).Gates, name))
                    idx = strcmp(app.data.(smpl_rmv).Gates, name);
                    app.data.(smpl_rmv).Gates(idx) = [];
                end
            end

            % Add pointers to this gate in samples.
            for smpl_add_idx=1:numel(to_add)
                smpl_add = to_add{smpl_add_idx};
                if ~isfield(app.data.(smpl_add), 'Gates')
                    app.data.(smpl_add).Gates = {};
                elseif any(strcmp(app.data.(smpl_add).Gates, name))
                    % Do nothing. Already pointer to it. Sanity Check.
                    continue;
                end
                app.data.(smpl_add).Gates(end + 1) = {name};
            end
        end

        function save_gate(app, num)
            % SAVE_GATE It is a function that saves last gate from given plot.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - num - Identifier of figure from which to pull gate and
            %       data which will be used in saving process.
            %
            % Modifies:
            %   - app - Specifically fields in app.data.(smpl), where smpl
            %       are all samples on which this gate was defined.\
            
            if  ~isfield(app.NewGate, 'internal__curr') || ~isfield(app.NewGate.internal__curr, num) || isempty(app.NewGate.internal__curr)
                errordlg('Could not find gates on current figure. Aborting.');
                return;
            end
            name = app.NewGate.internal__curr.(num);

            valid_gate_name = Helper.valid_gate(name);
            if iscell(valid_gate_name)
                valid_gate_name = valid_gate_name{1};
            end
            if ~isfield(app.NewGate, valid_gate_name)  % Sanity check
                errordlg("Something went wrong. Gate with that name no longer exists. Aborting");
                return;
            end
            if ~isvalid(app.NewGate.(valid_gate_name).obj)
                errordlg("Make sure your gate is currently plotted before saving it.");
                return;
            end
            wb = waitbar(0, 'Please wait. Saving the Gates.');

            % Currently plotted data
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

            Plt_Gating.update_gate(app, name, num);  % Update the gate to be currently plotted, on correct smpls and phenos.
            smpls = app.NewGate.(valid_gate_name).smpls;
            phenos = app.NewGate.(valid_gate_name).parent;
            type = app.NewGate.(valid_gate_name).type;

            % Precompute stuff for non-global neighborhoods (ones, starting with a tag).
            in_rsn = true(numel(phenos), 1);
            if all(startsWith(phenos, Constants.neigh_tag))
                % Determine if it's MFI RSN or CCN. Consistently for all Samples.
                for phn_idx = 1:numel(phenos)
                    phn = phenos{phn_idx};
                    for smpl_idx = 1:numel(smpls)
                        smpl = smpls{smpl_idx};
                        if ~isfield(app.data.(smpl), 'MFIRSN') || ...
                                ~ismember({phn}, app.data.(smpl).MFIRSN.Properties.VariableNames)
                            in_rsn(phn_idx) = false;
                        end
                    end
                end
            end

            waitbar(.05, wb, 'Saving to the samples');
            for smpl_idx=1:numel(smpls)
                waitbar(.05 + (smpl_idx * .9 / numel(smpls)), wb, ...
                    "Saving to the samples. Sample " + ...
                    num2str(smpl_idx), " out of ", num2str(numel(smpls))...
                );
                smpl = smpls{smpl_idx};
                if numel(phenos) == 1 && strcmp(phenos, 'Density/MFI RSN')
                    % Ignore channel names that are wonky
                    INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                    INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).MFIRSN.Properties.VariableNames);
                    smpl_dat = app.data.(smpl).MFIRSN(:, dat.Properties.VariableNames(INDIgnore));
                elseif numel(phenos) == 1 && strcmp(phenos, 'Density/MFI CCN')
                    % Ignore channel names that are wonky
                    INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                    INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).MFICCN.Properties.VariableNames);
                    smpl_dat = app.data.(smpl).MFICCN(:, dat.Properties.VariableNames(INDIgnore));
                elseif numel(phenos) == 1 && strcmp(phenos, 'All Cells')
                    % Ignore channel names that are wonky
                    INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                    INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).AllCells.Properties.VariableNames);
                    smpl_dat = app.data.(smpl).AllCells(:, dat.Properties.VariableNames(INDIgnore));
                elseif all(startsWith(phenos, Constants.neigh_tag))
                    smpl_dat = {};
                    for phn_idx = 1:numel(phenos)
                        phn = phenos{phn_idx};
                        % Load all the neighborhoods into dat.
                        if in_rsn(phn_idx)
                            % Ignore channel names that are wonky
                            INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                            INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).MFIRSN.Properties.VariableNames);
                            smpl_dat{phn_idx} = app.data.(smpl).MFIRSN(app.data.(smpl).MFIRSN.(phn), dat.Properties.VariableNames(INDIgnore));
                        else
                            % Ignore channel names that are wonky
                            INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                            INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).MFICCN.Properties.VariableNames);
                            smpl_dat{phn_idx} = app.data.(smpl).MFICCN(app.data.(smpl).MFICCN.(phn), dat.Properties.VariableNames(INDIgnore));
                        end
                    end
                else
                    smpl_dat = {};
                    % Pull the data that is to be plotted
                    for phn_idx = 1:numel(phenos)
                        phn = phenos{phn_idx};
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
                        % Ignore channel names that are wonky
                        INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                        INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).AllCells.Properties.VariableNames);
                        
                        smpl_dat{phn_idx} = app.data.(smpl).AllCells(app.data.(smpl).AllCells.(tag)==1, dat.Properties.VariableNames(INDIgnore));
                    end
                end
                X = dat.(app.NewGate.(valid_gate_name).GateAxis1);
                Y = dat.(app.NewGate.(valid_gate_name).GateAxis2);

                if strcmp(app.NewGate.(valid_gate_name).type, 'Polygon')
                    Logical = inpolygon(X, Y,app.NewGate.(valid_gate_name).pos(:,1), app.NewGate.(valid_gate_name).pos(:,2));
                    points = app.NewGate.(valid_gate_name).pos;
                elseif strcmp(app.NewGate.(valid_gate_name).type, 'Rectangle')
                    pos = getPosition(app.NewGate.(valid_gate_name).obj);
                    minX = pos(1); maxX = pos(1) + pos(3);
                    minY = pos(2); maxY = pos(2) + pos(4);
                    Logical = inpolygon(X,Y,[minX, maxX],[minY, maxY]);
                    points = [minX, minY;
                              minX, maxY;
                              maxX, maxY;
                              maxX, minY
                    ];
                end
                if ~iscell(smpl_dat)
                    %deal with some issues with mismatch between channel
                    %indeces between samples
                    INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                    INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, smpl_dat.Properties.VariableNames);
                    taginclude = dat.Properties.VariableNames(INDIgnore);
                    for ragi = 1:(numel(taginclude)-3)
                        if ~any(ismember(smpl_dat(:,taginclude), dat(Logical,taginclude)))
                            taginclude = taginclude(1:end-1);
                        elseif any(ismember(smpl_dat(:,taginclude), dat(Logical,taginclude)))
                            break 
                        end
                    end
                    % General Table Case
                    logical_smpl = ismember(smpl_dat(:,taginclude), dat(Logical, taginclude));
                    if numel(phenos) == 1 && strcmp(phenos, 'Density/MFI RSN')
                        % Add the gated data to the sample
                        valid_name = Helper.valid_neigh(name);
                        if iscell(valid_name)
                            valid_name = valid_name{1};
                        end
                        app.data.(smpl).MFIRSN.(valid_name) = logical_smpl;

                        continue;
                    elseif numel(phenos) == 1 && strcmp(phenos, 'Density/MFI CCN')
                        % Add the gated data to the sample
                        valid_name = Helper.valid_neigh(name);
                        if iscell(valid_name)
                            valid_name = valid_name{1};
                        end
                        app.data.(smpl).MFICCN.(valid_name) = logical_smpl;
                        continue;
                    elseif numel(phenos) == 1 && strcmp(phenos, 'All Cells')
                        short_name = strcat('All/', name);
                        if iscell(short_name)
                            short_name = short_name{1};
                        end
                        tag = Helper.get_tag(app, short_name, smpl);
                        if iscell(tag)
                            tag = tag{1};
                        end
                        path = strcat('Gate_All/', Helper.valid_gate(name));
                        if iscell(path)
                            path = path{1};
                        end

                        new_tree = tree(name, 'gate_points', points, ...
                            'gate_axes', {app.NewGate.(valid_gate_name).GateAxis1; ...
                                          app.NewGate.(valid_gate_name).GateAxis2}, ...
                            'gate_type', lower(type), ...
                            'tag', tag);

                        app.data.(smpl).AllCells.(tag) = logical_smpl;
                        app.data.(smpl).GateTags.(tag) = {path; short_name};
                        app.data.(smpl).tree = app.data.(smpl).tree.add_kid(new_tree);
                    end
                else
                    %% Phenotype Case or Neighborhood
                    if all(startsWith(phenos, Constants.neigh_tag))
                        %% Neighborhoods
                        valid_name = Helper.valid_neigh(name);
                        if iscell(valid_name)
                            valid_name = valid_name{1};
                        end

                        if isfield(app.data.(smpl), 'MFIRSN')
                            logical_rsn = zeros(size(app.data.(smpl).MFIRSN, 1), 1);
                        end
                        
                        if isfield(app.data.(smpl), 'MFICCN')
                            logical_ccn = zeros(size(app.data.(smpl).MFICCN, 1), 1);
                        end
                        
                        for phn_idx=1:numel(phenos)
                            phn = phenos{phn_idx};
                            
                            if in_rsn(phn_idx)
                                logical_rsn = logical_rsn | app.data.(smpl).MFIRSN.(phn);
                            else
                                logical_ccn = logical_ccn | app.data.(smpl).MFICCN.(phn);
                            end
                            
                            if any(logical_rsn)
                                INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                                INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).MFIRSN.Properties.VariableNames);
                                INDIgnore = dat.Properties.VariableNames(INDIgnore);
                                logical_rsn(logical_rsn == 1) = ismember(smpl_dat{phn_idx}, dat(Logical, INDIgnore));
                                app.data.(smpl).MFIRSN.(valid_name) = logical_rsn;
                            elseif any(logical_ccn)
                                INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                                INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).MFICCN.Properties.VariableNames);
                                INDIgnore = dat.Properties.VariableNames(INDIgnore);
                                logical_ccn(logical_ccn == 1) = ismember(smpl_dat{phn_idx}, dat(Logical, INDIgnore));
                                app.data.(smpl).MFICCN.(valid_name) = logical_ccn;
                            end
                            
                        end

                        continue;
                    else
                        %% Phenotype
                        valid_name = Helper.valid_gate(name);
                        for phn_idx=1:numel(phenos)
                            %deal with some issues with mismatch between channel
                            %indeces between samples
                            INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                            INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, smpl_dat{phn_idx}.Properties.VariableNames);
                            taginclude = dat.Properties.VariableNames(INDIgnore);
                            for taginclude_i = 1:(numel(taginclude)-3)
                                if ~any(ismember(smpl_dat{phn_idx}(:,taginclude), dat(Logical,taginclude)))
                                    taginclude = taginclude(1:end-1);
                                elseif any(ismember(smpl_dat{phn_idx}(:,taginclude), dat(Logical,taginclude)))
                                    break 
                                end
                            end
                            logical_smpl = ismember(smpl_dat{phn_idx}(:,taginclude), dat(Logical, taginclude));

                            phn = phenos{phn_idx};
                            short_name = split(phn, '/');
                            if strcmp(short_name(2), name)
                                continue;
                            end
                            short_name = strcat(short_name(2), '/', name);
                            if iscell(short_name)
                                short_name = short_name{1};
                            end
                            tag = Helper.get_tag(app, short_name, smpl);
                            if iscell(tag)
                                tag = tag{1};
                            end

                            new_path = Helper.gate_full2path(app, phn, smpl);
                            new_path = strcat(new_path, '/', valid_name);
                            if iscell(new_path)
                                new_path = new_path{1};
                            end

                            new_tree = tree(name, 'gate_points', points, ...
                                'gate_axes', {app.NewGate.(valid_gate_name).GateAxis1; ...
                                              app.NewGate.(valid_gate_name).GateAxis2}, ...
                                'gate_type', lower(app.NewGate.(valid_gate_name).type), ...
                                'tag', tag ...
                            );

                            parent_tag = Helper.gate_full2tag(app, phn, smpl);
                            if iscell(parent_tag)
                                parent_tag = parent_tag{1};
                            end
                            new_logical = app.data.(smpl).AllCells.(parent_tag);
                            new_logical(new_logical == 1) = logical_smpl;
                            app.data.(smpl).AllCells.(tag) = new_logical;
                            app.data.(smpl).GateTags.(tag) = {new_path; short_name};
                            try
                                app.data.(smpl).tree = app.data.(smpl).tree.add_kid(new_tree, new_path);
                            catch
                                print('Adding population to tree broke')
                            end
                        end
                    end
                end
            end
            waitbar(.95, wb, 'Updating Table for Current Plot.');
            if isfield(app.figOpts.phn_smpl_tab, num)
                app.figOpts.phn_smpl_tab.(num).Data = Plotting.update_phn_smpl_table(app, num, app.figOpts.phn_smpl_tab.(num).Data);
            end
            close(wb);
            Plt_Gating.update_gate(app, name, num);
        end

        function save_surf_gate(app, num, smpls, Type, name, pop_name, outside)
            %% SAVE_SURF_GATE Pull the currently plotted data
            
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
            X = dat.(vrnms{app.figOpts.xaxIF.(num).Value});
            Y = dat.(vrnms{app.figOpts.yaxIF.(num).Value});
            Z = dat.(vrnms{app.figOpts.zaxIF.(num).Value});

            % Get the currently plotted surface
            if ~iscell(smpls)
                smpls = {smpls};
            end
            for smpl_idx=1:numel(smpls)
                smpl = smpls{smpl_idx};
                Vol = app.data.(smpl).Surfaces.(Type).(name).Surf;

                % Define the new populations name
                pop_name = Helper.valid_gate(pop_name); % Make it compatible with other gates.
                if iscell(pop_name)
                    pop_name = pop_name{1};
                end

                % Either pull the surface or the plotted points within
                % the gating surface
                % Find the points that are in the surface
                if size(Vol.Points, 2) == 3  % 3D plot
                    IND = inShape(Vol, X, Y, Z);
                    if outside
                        IND = ~IND;
                    end
                    hold on
                    plot3(X(IND), Y(IND), Z(IND), '.r');
                    hold off
                else % 2D plot
                    IND = inShape(Vol, X, Y);
                    if outside
                        IND = ~IND;
                    end
                    hold on
                    plot(X(IND),Y(IND), '.r');
                    hold off
                end
                % Find the binary matrix for AllCells
                % Use this notation to Gate Cells
                GatedDat = ismember(app.data.(smpl).AllCells(:, dat.Properties.VariableNames), dat(IND, :), 'rows');
                for parent_idx=1:numel(app.figOpts.phenos.(num))
                    % Pull the name of the parent population
                    parent = app.figOpts.phenos.(num){parent_idx};

                    if strcmp('All Cells', parent)
                        parent_name = Helper.full_gate(app.data.(smpl).tree.name);
                        parent_path = strcat(app.data.(smpl).tree.name);
                        parent_logic = ones(size(GatedDat, 1), 1);
                    else
                        idx = find(strcmp({parent}, table2cell(app.data.(smpl).GateTags(2, :))));

                        parent_path = app.data.(smpl).GateTags{1, idx};
                        parent_name = app.data.(smpl).GateTags{2, idx};
                        parent_name = split(parent_name);
                        if numel(parent_name)>1
                            % This might break if the parent name isn't
                            % a cell array
                            parent_name = parent_name(2);
                        end
                        parent_logic = app.data.(smpl).AllCells.(app.data.(smpl).GateTags.Properties.VariableNames{idx});
                    end

                    % Define the full path and short name for the new gated
                    % population
                    path_to_gate = strcat(parent_path, '/', pop_name);
                    if iscell(path_to_gate)
                        path_to_gate = path_to_gate{1};
                    end

                    short_name = strcat(parent_name, '/', Helper.full_gate(pop_name));
                    if iscell(short_name)
                        short_name = short_name{1};
                    end
                    %% Add a new population to AllCells, Tree, GateTags and the Phenotype table

                    % Find the Current Number of Gates
                    Gates = app.data.(smpl).GateTags.Properties.VariableNames;

                    idx = find(strcmp({short_name}, table2cell(app.data.(smpl).GateTags(2, :))));
                    if isempty(idx)
                        tag = Helper.get_tag(app, short_name, smpl);
                        if iscell(tag)
                            tag = tag{1};
                        end
                        % Add the new population to GateTags
                        app.data.(smpl).GateTags.(tag) = {path_to_gate; short_name};
                        new_tree = tree(Helper.full_gate(pop_name), 'gate_type', 'logic', 'tag', tag);
                        app.data.(smpl).tree = app.data.(smpl).tree.add_kid(new_tree, parent_path);
                    else
                        tag = Gates{idx};
                    end

                    % Add the new population to Phenotype Table
                    app.data.(smpl).AllCells.(tag) = GatedDat & parent_logic;
                end
            end
        end
    end % End of Methods
end % End of Class Definition