classdef Plotting

% Plotting defines a suite of functions which plot and visualize data
%
% Last edited by CS 2019/01/24

% Warning supressions.
%#ok<*IMPOLY>
%#ok<*IMRECT>

    methods (Static)
        function Import_Definitions_Func
% % %             import Helper.*;
% % %             import Constants.*;
% % %             import Plt_Helper.*;
% % %             import Plt_Gating.*;
            set(groot, 'DefaultTextInterpreter', 'none')
            set(groot, 'DefaultLegendInterpreter', 'none')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main Plot/Figure Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function func_newfig(app, varargin)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % func_newfig creates a new figure, including all controls
            % which CytoMAP gives to user on a plot
            %
            % Input:
            %   - app - Instance of CytoMAP
            %
            % Key-word Input (optional):
            %   - X - char - default: 'X' - what to plot on the X axis, and
            %       what to fill into X-axis dropdown
            %   - Y - char - default: 'Y' - what to plot on the Y axis, and
            %       what to fill into Y-axis dropdown
            %   - Z - char - default: 'Z' - what to plot on the Z axis, and
            %       what to fill into Z-axis dropdown
            %   - C - char - default: 'None' - what to plot on the Color
            %       axis, and what to fill into Color-axis dropdown
            %   - P - cell, string, char - default: 'All Cells' - which
            %       phenotypes to include on a plot, and which ones to mark
            %       in the Sample/Phenotype Table.
            %   - S - cell, string, char - default: app.DataN.Items{1} -
            %       which samples to include on a plot, and which ones to
            %       mark in the Sample/Phenotype Table.
            %   - Sldr - char - default: 'Z' - what axis to have a slider
            %       which filters the plotted data on.
            %
            % Modifies:
            %   - app - Specifically app.figOpts.___.(num), where ___ is
            %       bunch of different things CytoMAP keeps track of on the
            %       plot, and num is identifier of the figure.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~Helper.any_sample(app)
                return;
            end
            
            %% Create the figure
            fig = figure;

            fig.Name = 'RSZ';
            Helper.func_SetCLR(app, fig, 'figure')
            fig.Position(3:4) = [725 450];
            num = [ 'fig' num2str(fig.Number)]; %(Use this to assign handles)
            fig.SizeChangedFcn = @(dd, ~) resize_figure_wrap(dd, app);
            fig.CloseRequestFcn = @(dd, ~) Plotting.close_figure(dd, app, num);

            % Make fields in figOpts which are recently added.
            new_fields = {'opts',...
                'xaxTXT', 'yaxTXT', 'zaxTXT', 'caxTXT', 'phnTXT', 'sampleTXT', 'sldrIFTXT', ...
                'phenos', 'smpls', ...
                'plt_ax', ...
                'phn_smpl_tab', 'close_phn_smpl_tab_btn', ...
                'order_list_btn', 'order_list', 'order_stack', 'order_items', ...
                'dat', 'vrnms'
            };

            for f_i=1:numel(new_fields)
                if ~isfield(app.figOpts, new_fields{f_i})
                    app.figOpts.(new_fields{f_i}) = struct;
                end
            end

            if ~isfield(app.figOpts.smpls, num)
                app.figOpts.smpls.(num) = {app.DataN.Value};
            end

            % Get current phenotypes in selected samples if they do not exist
            if ~isfield(app.figOpts.phenos, num)
                app.figOpts.phenos.(num) = 'All Cells';
            end

            app.figOpts.bgclr.(num) = app.GUIOPTS.bgclr;
            app.figOpts.txtclr.(num) = app.GUIOPTS.txtclr;

            % Define defualt arguments for varargin
            app.figOpts.defaults.(num) = struct( ...
                'X', 'X', ...
                'Y', 'Y', ...
                'Z', 'Z', ...
                'C', 'None', ...
                'P', 'All Cells', ...
                'Sldr', 'Z', ...
                'S', app.DataN.Items{1} ...
            );

            % Check that varargin is of even length
            if mod(length(varargin), 2) ~= 0
                error('Uneven number of input arguments! func_newfig needs propertyName/propertyValue pairs, after app argument')
            end

            % Process key, i.e. change default values.
            for pair = reshape(varargin, 2, [])
                if isfield(app.figOpts.defaults.(num), pair{1})
                    app.figOpts.defaults.(num).(pair{1}) = pair{2}; % Change default axis value to the given value
                else  % Key is not in the defaults
                    error('%s is not a recognized parameter name', pair{1})
                end
            end
            
            if ~iscell(app.figOpts.defaults.(num).S)
                app.figOpts.defaults.(num).S = {app.figOpts.defaults.(num).S};
            end
            
            %% Create Menu bar options
            % Create Menu bar
            ExportMenu = uimenu(fig);
            ExportMenu.Text = 'Export';

            % Create an export data option
            ExportPltDat = uimenu(ExportMenu);
            ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, num);
            ExportPltDat.Text = 'Export Plot Data to .csv';

            % Create an export figure option
            ExportPltDat = uimenu(ExportMenu);
            ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_figure(app, num);
            ExportPltDat.Text = 'Export Plot';

            %% Create General plot menu options
            app.PlotMenu.(num).Main = uimenu(fig);
            app.PlotMenu.(num).Main.Text = 'Plot';

            % Create asmooth data options
            app.PlotMenu.(num).Smooth = uimenu(app.PlotMenu.(num).Main);
            app.PlotMenu.(num).Smooth.Text = 'Plot Smoothed Data';
            
            app.PlotMenu.(num).Smooth_mean = uimenu(app.PlotMenu.(num).Smooth);
            app.PlotMenu.(num).Smooth_mean.Text = 'Smooth with running mean';
            app.PlotMenu.(num).Smooth_mean.Checked = 'off';
            app.PlotMenu.(num).Smooth_mean.MenuSelectedFcn = ...
                @(~, ~) func_Checked_PLT(app, num,app.PlotMenu.(num).Smooth_mean);

            app.PlotMenu.(num).Smooth_Gauss = uimenu(app.PlotMenu.(num).Smooth);
            app.PlotMenu.(num).Smooth_Gauss.Text = 'Smooth with gaussian';
            app.PlotMenu.(num).Smooth_Gauss.Checked = 'off';
            app.PlotMenu.(num).Smooth_Gauss.MenuSelectedFcn = ...
                @(~, ~) func_Checked_PLT(app, num,app.PlotMenu.(num).Smooth_Gauss);

            app.PlotMenu.(num).Smooth_CumSum = uimenu(app.PlotMenu.(num).Smooth);
            app.PlotMenu.(num).Smooth_CumSum.Text = 'Smooth with normalized cumulative sum';
            app.PlotMenu.(num).Smooth_CumSum.Checked = 'off';
            app.PlotMenu.(num).Smooth_CumSum.MenuSelectedFcn = ...
                @(~, ~) func_Checked_PLT(app, num,app.PlotMenu.(num).Smooth_CumSum);
            
            %% Create Plot Gate menu options
            app.PlotMenu.(num).GateMain = uimenu(app.PlotMenu.(num).Main);
            app.PlotMenu.(num).GateMain.Text = 'Plot Gates';
            app.PlotMenu.(num).Gate.UDG = struct;
            for i=1:numel(app.DataN.Items)
                app.PlotMenu.(num).Gate.UDG.(['SelectSMPL' (num2str(i))]) = uimenu(app.PlotMenu.(num).GateMain);
                app.PlotMenu.(num).Gate.UDG.(['SelectSMPL' (num2str(i))]).Text = app.DataN.Items{i};
                % If there are user defined gates put them in the drop down menu
                if isfield(app.data.(app.DataN.Items{i}), 'Gates')
                    gatenms = app.data.(app.DataN.Items{i}).Gates;
                    for j = 1:numel(gatenms)
                        SelectGate = uimenu(app.PlotMenu.(num).Gate.UDG.(['SelectSMPL' (num2str(i))]));
                        SelectGate.MenuSelectedFcn = @(~,~) Plt_Gating.update_gate(app, gatenms{j}, num);
                        SelectGate.Text = gatenms{j};
                    end
                end
            end

            %% Create plot Surface menu options
            app.PlotMenu.(num).SurfacesMain = uimenu(app.PlotMenu.(num).Main);
            app.PlotMenu.(num).SurfacesMain.Text = 'Plot Surface';
            app.PlotMenu.(num).Surfaces.UDS = struct;
            for i=1:numel(app.DataN.Items)
                app.PlotMenu.(num).Surfaces.UDS.(['SMPL' (num2str(i))]) = uimenu(app.PlotMenu.(num).SurfacesMain);
                app.PlotMenu.(num).Surfaces.UDS.(['SMPL' (num2str(i))]).Text = app.DataN.Items{i};
                if sum(contains(fieldnames( app.data.(app.DataN.Items{i})), 'Surfaces'))==0
                    app.data.(app.DataN.Items{i}).Surfaces = struct;
                end
                % Check for user defined surfaces Surfaces.UDS
                if sum(contains(fieldnames( app.data.(app.DataN.Items{i}).Surfaces), 'UDS'))==0
                    app.data.(app.DataN.Items{i}).Surfaces.UDS = struct;
                else % If there are user defined surfaces put them in the drop down menu
                    surfnms = fieldnames(app.data.(app.DataN.Items{i}).Surfaces.UDS);
                    for j = 1:numel(surfnms)
                        SelectSurf = uimenu(app.PlotMenu.(num).Surfaces.UDS.(['SMPL' (num2str(i))]));
                        SelectSurf.MenuSelectedFcn = @(~,~) Plt_Helper.frontend_surf(app, num, app.DataN.Items{i}, surfnms{j}, true, 'UDS');
                        SelectSurf.Text = surfnms{j};
                    end
                end
            end

            %% Create plot Regions menu options
            app.PlotMenu.(num).RegionSurfacesMain = uimenu(app.PlotMenu.(num).Main);
            app.PlotMenu.(num).RegionSurfacesMain.Text = 'Plot Regions';
            app.PlotMenu.(num).Surfaces.RSNMain = uimenu(app.PlotMenu.(num).RegionSurfacesMain);
            app.PlotMenu.(num).Surfaces.RSNMain.Text = 'Raster Scanned Neighborhoods';
            app.PlotMenu.(num).Surfaces.CCNMain = uimenu(app.PlotMenu.(num).RegionSurfacesMain);
            app.PlotMenu.(num).Surfaces.CCNMain.Text = 'Cell Centered Neighborhoods';
            app.PlotMenu.(num).Surfaces.RSN = struct;
            app.PlotMenu.(num).Surfaces.CCN = struct;
            for i=1:numel(app.DataN.Items)
                app.PlotMenu.(num).Surfaces.RSN.(['SMPL' (num2str(i))]) = uimenu(app.PlotMenu.(num).Surfaces.RSNMain);
                app.PlotMenu.(num).Surfaces.RSN.(['SMPL' (num2str(i))]).Text = app.DataN.Items{i};

                app.PlotMenu.(num).Surfaces.CCN.(['SMPL' (num2str(i))]) = uimenu(app.PlotMenu.(num).Surfaces.CCNMain);
                app.PlotMenu.(num).Surfaces.CCN.(['SMPL' (num2str(i))]).Text = app.DataN.Items{i};

                if sum(contains(fieldnames( app.data.(app.DataN.Items{i}).Surfaces), 'RSN'))==0
                    app.data.(app.DataN.Items{i}).Surfaces.RSN = struct;
                else % If there are surfaces put them in the drop down menu
                    RSNsurfnms = fieldnames(app.data.(app.DataN.Items{i}).Surfaces.RSN);
                    for j = 1:numel(RSNsurfnms)
                        SelectSurf = uimenu(app.PlotMenu.(num).Surfaces.RSN.(['SMPL' (num2str(i))]));
                        SelectSurf.MenuSelectedFcn = @(~,~) Plt_Helper.frontend_surf(app, num, app.DataN.Items{i}, RSNsurfnms{j}, true, 'RSN');
                        SelectSurf.Text = RSNsurfnms{j};
                    end
                end

                if sum(contains(fieldnames( app.data.(app.DataN.Items{i}).Surfaces), 'CCN'))==0
                    app.data.(app.DataN.Items{i}).Surfaces.CCN = struct;
                else % If there are surfaces put them in the drop down menu
                    CCNsurfnms = fieldnames(app.data.(app.DataN.Items{i}).Surfaces.CCN);
                    for j = 1:numel(CCNsurfnms)
                        SelectSurf = uimenu(app.PlotMenu.(num).Surfaces.CCN.(['SMPL' (num2str(i))]));
                        SelectSurf.MenuSelectedFcn = @(~,~) Plt_Helper.frontend_surf(app, num, app.DataN.Items{i}, CCNsurfnms{j}, true, 'CCN');
                        SelectSurf.Text = CCNsurfnms{j};
                    end
                end
            end

            %% Create Plot Polygons menu options
            app.PlotMenu.(num).PolygonMain = uimenu(app.PlotMenu.(num).Main);
            app.PlotMenu.(num).PolygonMain.Text = 'Plot Polygons';
            app.PlotMenu.(num).Polygon.UDP = struct;
            polynms = fieldnames(app.polygons);
            for i=1:numel(app.DataN.Items)
                app.PlotMenu.(num).Polygon.UDP.(['SelectSMPL' (num2str(i))]) = uimenu(app.PlotMenu.(num).PolygonMain);
                app.PlotMenu.(num).Polygon.UDP.(['SelectSMPL' (num2str(i))]).Text = app.DataN.Items{i};
                % If there are user defined polygons put them in the drop down menu
                if ~isempty(polynms)
                    for j = 1:numel(polynms)
                        SelectGate = uimenu(app.PlotMenu.(num).Polygon.UDP.(['SelectSMPL' (num2str(i))]));
                        SelectGate.MenuSelectedFcn = @(~,~) Plt_Helper.func_poly(app, polynms{j});
                        SelectGate.Text = polynms{j};
                    end
                end
            end
            
            %% Create Highlight Region Option
            
            app.PlotMenu.(num).HighlightReg = uimenu(app.PlotMenu.(num).Main);
            app.PlotMenu.(num).HighlightReg.Text = 'Highlight Regions';
            app.PlotMenu.(num).HighlightReg.MenuSelectedFcn = @(~,~) Plt_Helper.highight_reg_func(app, num);
%             app.PlotMenu.(num).HRegs = struct;       

            app.PlotMenu.(num).HighlightReg = uimenu(app.PlotMenu.(num).Main);
            app.PlotMenu.(num).HighlightReg.Text = 'Highlight Neighborhoods';
            app.PlotMenu.(num).HighlightReg.MenuSelectedFcn = @(~,~) Plt_Helper.highight_Neigh_func(app, num);
            
            %% Create toolbar options

            % Make sure the standard toolbar is displayed
            fig.ToolBar = 'figure';

            % Find the standard toolbar handle
% % %             toolbar = struct;
            Main = findall(fig,'Type','uitoolbar');

            % Clear fugure
            ClearFig  =  uipushtool('Parent',Main);
            ClearFig.TooltipString = 'Clear Plot';
            ClearFig.CData = app.figOpts.icons.ClearFig;
            ClearFig.Separator = 'on';
            ClearFig.ClickedCallback = 'cla';

            % Update figure
            UpdateFig  =  uipushtool('Parent',Main);
            UpdateFig.TooltipString = 'Update Plot';
            UpdateFig.CData = app.figOpts.icons.RefreshPlot;
            UpdateFig.ClickedCallback =  @(~,~) Plt_Helper.func_refresh(app, num);

            % Make a new rectangular gate
            NewRectGate  =  uipushtool('Parent',Main);
            NewRectGate.TooltipString = 'New Rectangular Gate';
            NewRectGate.CData = app.figOpts.icons.RectGate;
            NewRectGate.Separator = 'on';
            NewRectGate.ClickedCallback = @(~,~) Plt_Gating.new_gate(app, num, 'RectGate');

            % Make a new Polygon gate
            NewRectGatePoly  =  uipushtool('Parent',Main);
            NewRectGatePoly.TooltipString = 'New Polygon Gate';
            NewRectGatePoly.CData = app.figOpts.icons.PolyGate;
            NewRectGatePoly.Separator = 'off';
            NewRectGatePoly.ClickedCallback = @(~,~) Plt_Gating.new_gate(app, num, 'PolyGate');

            % Export Gated population
            SaveGate = uipushtool('Parent',Main);
            SaveGate.TooltipString = 'Save Last Chosen Gate';
            SaveGate.CData = app.figOpts.icons.Save;
            SaveGate.Separator = 'off';
            SaveGate.ClickedCallback = @(~,~) Plt_Gating.save_gate(app, num);

            % Axes even
            eqax = uitoggletool('Parent',Main);
            eqax.TooltipString = 'Square Axes';
            eqax.CData = app.figOpts.icons.EvenAxes;
            app.figOpts.eqax.(num) = true;
            eqax.ClickedCallback = @(~,~) Plt_Helper.pltax(eqax, app, num);
            eqax.State= 'on';

            % Plot 3D
            D3 = uitoggletool('Parent',Main);
            D3.TooltipString = 'Make plot 3D';
            D3.CData = app.figOpts.icons.Plot3D;
            app.figOpts.D3.(num) = false;
            D3.ClickedCallback = @(~,~) Plt_Helper.plt3D(D3, app, num);
            D3.State= 'off';

            % Invert Colors
            IC = uitoggletool('Parent',Main);
            IC.TooltipString = 'Invert Colors';
            IC.CData = app.figOpts.icons.INVERT;
            if strcmp(app.figOpts.bgclr.(num), 'w')
                app.figOpts.IC.(num)=false;
                IC.State= 'off';
            else
                app.figOpts.IC.(num)=true;
                IC.State= 'on';
            end
            IC.ClickedCallback = @(~,~) Plt_Helper.pltIC(IC, app, num);

            % Rotate figure TODO: make a unique icon for this
            % (This rotates the data in the figure for display reasons it doesn't change the original data)
            UpdateFig  =  uipushtool('Parent',Main);
            UpdateFig.TooltipString = 'Rotate about Z axis';
            UpdateFig.CData = app.figOpts.icons.Rotate;
            UpdateFig.ClickedCallback =  @(~,~) Plt_Helper.frontend_rotate(app);
            % New Spots
            PltSpts  =  uipushtool('Parent',Main);
            PltSpts.TooltipString = 'Draw Spots';
            PltSpts.CData = app.figOpts.icons.NewPoints;
            PltSpts.Separator = 'on';
            PltSpts.ClickedCallback =  @(~, ~) Plt_Helper.func_Spots(app);
            % New Polygons
            PltPoly  =  uipushtool('Parent',Main);
            PltPoly.TooltipString = 'Draw Polygons';
            PltPoly.CData = app.figOpts.icons.NewPoly;
            PltPoly.ClickedCallback =  @(~, ~) Plt_Helper.func_poly(app, '');
            % Plot Points
            PltSpts  =  uipushtool('Parent',Main);
            PltSpts.TooltipString = 'Plot Spots';
            PltSpts.CData = app.figOpts.icons.PlotPoints;
            PltSpts.ClickedCallback =  @(~, ~) Plt_Helper.plotpts(app, num);
            % Make a surface around a group of points
            VolWrap  =  uipushtool('Parent',Main);
            VolWrap.TooltipString = 'Make Surface';
            VolWrap.CData = app.figOpts.icons.Wrap;
            VolWrap.Separator = 'on';
            VolWrap.ClickedCallback =  @(~, ~) Plt_Helper.frontend_surf(app, num, app.figOpts.smpls.(num), [], false, 'UDS');
            % Generate Random Points
            VolWrap  =  uipushtool('Parent',Main);
            VolWrap.TooltipString = 'Generate Random Points';
            VolWrap.CData = app.figOpts.icons.RandPoints;
            VolWrap.Separator = 'off';
            VolWrap.ClickedCallback =  @(~, ~) Plt_Helper.frontend_random_pts(app, num);
            
            %% Create some buttons
            % Create an options button
            app.figOpts.opts.(num) = uicontrol('Style', 'pushbutton', ...
                'String', 'Options',...
                'Position', [Constants.plt_menu_x0 20 50 30],...
                'Callback', @(~, ~) Plt_Helper.plot_options(app, num));
            Helper.func_SetCLR(app,  app.figOpts.opts.(num), 'UICpushbutton')
            % Create a show table button
            app.figOpts.smplPhenoBtn.(num) = uicontrol('Style', 'pushbutton', ...
                'String', 'Show Table', ...
                'Callback', @(~, ~) Plotting.show_phn_smpl_table(app, num, fig));
            Helper.func_SetCLR(app,  app.figOpts.smplPhenoBtn.(num), 'UICpushbutton')
            app.figOpts.smplPhenoBtn.(num).Position = [Constants.plt_menu_x0, fig.Position(4) - Constants.plt_menu_y0 - 30, 150, 40];
            app.figOpts.smplPhenoBtn.(num).Position(3:4) = app.figOpts.smplPhenoBtn.(num).Extent(3:4);
            
            %% Put axis options in the figure
            app.figOpts.xaxIF.(num) = uicontrol('Style', 'popup','Callback', @(~,~) Plotting.func_plot(app, num));
            app.figOpts.xaxIF.(num).Position= [210 10 100 15];
            app.figOpts.xaxTXT.(num) = uicontrol('Style','text');
            app.figOpts.xaxTXT.(num).Position = [160 10 50 15];
            app.figOpts.xaxTXT.(num).String = 'X-Axis: ';
            app.figOpts.xaxTXT.(num).Position(3:4) = app.figOpts.xaxTXT.(num).Extent(3:4);
            app.figOpts.xaxTXT.(num).HorizontalAlignment = 'right';
            Helper.func_SetCLR(app, app.figOpts.xaxTXT.(num), 'UICpopup')

            app.figOpts.yaxIF.(num) = uicontrol('Style', 'popup','Callback', @(~,~) Plotting.func_plot(app, num));
            app.figOpts.yaxIF.(num).Position= [210 30 100 15];
            app.figOpts.yaxTXT.(num) = uicontrol('Style','text');
            app.figOpts.yaxTXT.(num).Position = [160 30 50 15];
            app.figOpts.yaxTXT.(num).String = 'Y-Axis: ';
            app.figOpts.yaxTXT.(num).Position(3:4) = app.figOpts.yaxTXT.(num).Extent(3:4);
            app.figOpts.yaxTXT.(num).HorizontalAlignment = 'right';
            Helper.func_SetCLR(app, app.figOpts.yaxTXT.(num), 'UICpopup')

            app.figOpts.zaxIF.(num) = uicontrol('Style', 'popup','Callback', @(~,~) Plotting.func_plot(app, num));
            app.figOpts.zaxIF.(num).Position= [370 30 100 15];
            app.figOpts.zaxTXT.(num) = uicontrol('Style','text');
            app.figOpts.zaxTXT.(num).Position = [320 30 50 15];
            app.figOpts.zaxTXT.(num).String = 'Z-Axis: ';
            app.figOpts.zaxTXT.(num).Position(3:4) = app.figOpts.zaxTXT.(num).Extent(3:4);
            app.figOpts.zaxTXT.(num).HorizontalAlignment = 'right';
            Helper.func_SetCLR(app, app.figOpts.zaxTXT.(num), 'UICpopup')

            app.figOpts.caxIF.(num) = uicontrol('Style', 'popup','Callback', @(~,~) Plotting.func_plot(app, num));
            app.figOpts.caxIF.(num).Position= [370 10 100 15];
            app.figOpts.caxTXT.(num) = uicontrol('Style','text');
            app.figOpts.caxTXT.(num).Position = [320 10 50 15];
            app.figOpts.caxTXT.(num).String = 'C-Axis: ';
            app.figOpts.caxTXT.(num).Position(3:4) = app.figOpts.caxTXT.(num).Extent(3:4);
            app.figOpts.caxTXT.(num).HorizontalAlignment = 'right';
            Helper.func_SetCLR(app, app.figOpts.caxTXT.(num), 'UICpopup')

            % Add a slidebar with its own axis
            app.figOpts.SldrZ.(num) = uicontrol('style','slider');
            app.figOpts.SldrZ.(num).Min = 0;
            app.figOpts.SldrZ.(num).Value = 0.5;
            app.figOpts.SldrZ.(num).Max = 1;
            app.figOpts.SldrZ.(num).Position = [720 65 15 365];
            app.figOpts.SldrZ.(num).Visible = 'on';
            addlistener(app.figOpts.SldrZ.(num), 'Value', 'PostSet', @(~, ~) Plotting.func_plot(app, num));
            Helper.func_SetCLR(app, app.figOpts.SldrZ.(num), 'UICpopup')
            
            % add slidebar label
            app.figOpts.SldrTXT.(num) = uicontrol('Style','text');
            app.figOpts.SldrTXT.(num).Position = [650 430 200 15];
            app.figOpts.SldrTXT.(num).Position(3:4) = app.figOpts.SldrTXT.(num).Extent(3:4);
            app.figOpts.SldrTXT.(num).String = 'Slidebar: ';
            app.figOpts.SldrTXT.(num).HorizontalAlignment = 'left';
            app.figOpts.SldrTXT.(num).Visible = 'on';
            Helper.func_SetCLR(app,  app.figOpts.SldrTXT.(num), 'UICpopup')

            % Add a Slider Axis Selection
            app.figOpts.sldrIF.(num) = uicontrol('Style', 'popup','Callback', @(~,~) Plotting.func_plot(app, num));
            app.figOpts.sldrIF.(num).Position= [555 30 100 15];
            app.figOpts.sldrIFTXT.(num) = uicontrol('Style','text');
            app.figOpts.sldrIFTXT.(num).Position = [475 22 50 15];
            app.figOpts.sldrIFTXT.(num).HorizontalAlignment = 'right';
            app.figOpts.sldrIFTXT.(num).String = 'Slider-Axis: ';
            app.figOpts.sldrIFTXT.(num).Position(3:4) = app.figOpts.sldrIFTXT.(num).Extent(3:4);
            Helper.func_SetCLR(app,  app.figOpts.sldrIFTXT.(num), 'UICpopup')

            % Add a Z-range selection
            app.figOpts.sldrRNG.(num) = uicontrol('Style','edit');
            app.figOpts.sldrRNG.(num).String = '1';
            app.figOpts.sldrRNG.(num).Position = [600 1 85 15];
            Helper.func_SetCLR(app,  app.figOpts.sldrRNG.(num), 'UICpopup')

            % Add label for Z range
            app.figOpts.sldrRNGTXT.(num) = uicontrol('Style','text');
            app.figOpts.sldrRNGTXT.(num).Position = [475 1 70 15];
            app.figOpts.sldrRNGTXT.(num).String = 'Slider-Axis Range: ';
            app.figOpts.sldrRNGTXT.(num).Position(3:4) = app.figOpts.sldrRNGTXT.(num).Extent(3:4);
            app.figOpts.sldrRNGTXT.(num).HorizontalAlignment = 'right';
            app.figOpts.sldrRNGTXT.(num).Visible = 'on';
            Helper.func_SetCLR(app,  app.figOpts.sldrRNGTXT.(num), 'UICpopup')

% % %             % create a menu that responds to right clicks
% % %             rgtclk = uicontextmenu(fig);
% % %             % Assign the uicontextmenu to the figure
% % %             fig.UIContextMenu = rgtclk;
% % %             % Create child menu items for the uicontextmenu
% % %             uimenu(rgtclk,'Label','This does nothing yet','Callback',@rightclick);

            ax = gca;
            ax.Position(4) = ax.Position(4)-0.1*ax.Position(4);
            ax.Position(2) = ax.Position(2)+0.1;
            
            % Update the figure dropdown menus (This is takes too long)
            app.figOpts.phenos.(num) = app.figOpts.defaults.(num).P; 
            Plt_Helper.UpdateFigMenu(app, num);  

            % Change Plot Axes to default axes
            app.figOpts.xaxIF.(num).Value  = find(strcmp(app.figOpts.xaxIF.(num).String, app.figOpts.defaults.(num).X));
            app.figOpts.yaxIF.(num).Value  = find(strcmp(app.figOpts.yaxIF.(num).String, app.figOpts.defaults.(num).Y));
            app.figOpts.zaxIF.(num).Value  = find(strcmp(app.figOpts.zaxIF.(num).String, app.figOpts.defaults.(num).Z));
            app.figOpts.caxIF.(num).Value  = find(strcmp(app.figOpts.caxIF.(num).String, app.figOpts.defaults.(num).C));
            app.figOpts.sldrIF.(num).Value = find(strcmp(app.figOpts.sldrIF.(num).String, app.figOpts.defaults.(num).Sldr));

            app.figOpts.smpls.(num) = app.figOpts.defaults.(num).S;
            if isStringScalar(app.figOpts.defaults.(num).P) && strcmp(app.figOpts.defaults.(num).P, "AllCells")
                [~, app.figOpts.phenos.(num)] = Helper.get_gates(app, app.figOpts.defaults.(num).S);
            else
                app.figOpts.phenos.(num) = app.figOpts.defaults.(num).P;
            end

            app.figOpts.xaxIFTMP.(num) = 'NewFig';
            app.figOpts.yaxIFTMP.(num) = 'NewFig';
            app.figOpts.zaxIFTMP.(num) = 'NewFig';
            app.figOpts.caxIFTMP.(num) = 'NewFig';
            app.figOpts.sldrIFTMP.(num) = 'NewFig';
            app.figOpts.phnTMP.(num) = 'NewFig';
            app.figOpts.sampleTMP.(num) = 'NewFig';
            app.figOpts.View.(num) = [0 90];
            
            % resize the plot window
            Plotting.resize_figure(fig, app);

            % Generate text if the user clicks on a point
            Cursor_obj = datacursormode(fig);
            Cursor_obj.Enable = 'on';
            Cursor_obj.UpdateFcn = @(obj, event_obj) Plt_Helper.PlotClickFunc(obj, event_obj, app, num);
%             Cursor_obj.UpdateFcn = {@Plt_Helper.PlotClickFunc, app, num};

            % Update the figure
            Plotting.func_plot(app, num);

            % Fix the figure size
            Plotting.resize_figure(fig, app)

            %% UX function wrap-arounds
            function resize_figure_wrap(dd, app)
                persistent last_pos;
                if ~isempty(last_pos) && all(abs(last_pos(3:4)-dd.Position(3:4))<1) && ispc
                    return;
                else
                    Plotting.resize_figure(dd, app);
                end
                last_pos = dd.Position;
            end

            function func_Checked_PLT(app, num, object)
                switch object.Checked
                    case 'on'
                        object.Checked = 'off';
                        Plt_Helper.func_refresh(app, num);
                    case 'off'
                        object.Checked = 'on';
                        Plt_Helper.func_refresh(app, num);
                end
            end
        end

        function func_plot(app, num, overlay)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % func_plot actually plots stuff on the given figure according
            % to user choices.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - num - Figure Identifier (For example 'fig1')
            %
            % Modifies:
            %   - app - Specifically all fields in app.figOpts.___.(num)
            %       related to what is currently plotted, and plot choices.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % This is to prevent func_plot from calling itself when it
            % changes the z-axis slider value
            if strcmp(app.figOpts.SldrTXT.(num).String , 'StopListner')
                return
            end
            
            if nargin<3
                overlay = 0;
            else
                % If there is an overlayed image keep it
                if overlay == 1
                    ax = gca;
                    for ch_i = 1:numel(ax.Children)
                        if strcmp(ax.Children(ch_i).Type, 'image')
                            imgax = struct;
                            imgax.XData = ax.Children(ch_i).XData;
                            imgax.YData = ax.Children(ch_i).YData;
                            imgax.CData = ax.Children(ch_i).CData;
                            imgax.ing = ch_i;
                        end
                    end
                end
            end

            %% Initialize defaults and check axes
            % Make sure the selected phenotype and sample are cells
            if ~iscell(app.figOpts.phenos.(num))
                app.figOpts.phenos.(num) = {app.figOpts.phenos.(num)};
            end
            if ~iscell(app.figOpts.smpls.(num))
                app.figOpts.smpls.(num) = {app.figOpts.smpls.(num)};
            end
            
            %check if any axes or phenotypes changed
            app.figOpts.defaults.(num).xchnd = ~strcmp(app.figOpts.xaxIFTMP.(num), app.figOpts.xaxIF.(num).String{app.figOpts.xaxIF.(num).Value});
            app.figOpts.defaults.(num).ychnd = ~strcmp(app.figOpts.yaxIFTMP.(num), app.figOpts.yaxIF.(num).String{app.figOpts.yaxIF.(num).Value});
            app.figOpts.defaults.(num).zchnd = ~strcmp(app.figOpts.zaxIFTMP.(num), app.figOpts.zaxIF.(num).String{app.figOpts.zaxIF.(num).Value});
            % Why would this ever have more than one value?
            if numel(app.figOpts.caxIF.(num).Value)>1
                app.figOpts.caxIF.(num).Value = app.figOpts.caxIF.(num).Value(1);
            elseif isempty(app.figOpts.caxIF.(num).Value)
                app.figOpts.caxIF.(num).Value = find(strcmp(app.figOpts.caxIF.(num).String, "None"));
            end
            app.figOpts.defaults.(num).cchnd = ~strcmp(app.figOpts.caxIFTMP.(num), app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value});
            app.figOpts.defaults.(num).sldrchnd = ~strcmp(app.figOpts.sldrIFTMP.(num), app.figOpts.sldrIF.(num).String{app.figOpts.sldrIF.(num).Value});

% % % % Check if phenos/smpls sets are equal
% % % [{'phnTMP'}; app.figOpts.phnTMP.(num); {'phenos'}; app.figOpts.phenos.(num)]
            
            app.figOpts.defaults.(num).pchnd = ~Helper.setequal(app.figOpts.phnTMP.(num), app.figOpts.phenos.(num));
            app.figOpts.defaults.(num).schnd = ~Helper.setequal(app.figOpts.sampleTMP.(num), app.figOpts.smpls.(num));
            
            % Preserve the current viewing angle
            ax = gca;
            app.figOpts.View.(num) = ax.View;
            
            % Pull the sorted model names
            ModelNames = fieldnames(app.net);
            
            if numel(app.figOpts.phenos.(num)) > 1 || numel(app.figOpts.smpls.(num)) > 1
                % Plot multiple samples at a time.
                cla
                
                hold on
                % if you change phenotypes
                if app.figOpts.defaults.(num).pchnd
                    if strcmp(app.figOpts.phnTMP.(num),'NewFig')
                        app.figOpts.phenos.(num) = app.figOpts.defaults.(num).P;
                    else
                        app.figOpts.defaults.(num).P = app.figOpts.phenos.(num);
                    end
                end
                % if you changed samples
                if app.figOpts.defaults.(num).schnd
                    if strcmp(app.figOpts.sampleTMP.(num),'NewFig')
                        app.figOpts.smpls.(num) = app.figOpts.defaults.(num).S;
                    else
                        app.figOpts.defaults.(num).S = app.figOpts.smpls.(num);
                    end
                    
                end
                % Find the number of samples that have been overlayed
                if iscell(app.figOpts.defaults.(num).S)
                    N_smp = numel(app.figOpts.defaults.(num).S);
                else
                    N_smp = 1;
                end
                % Find the number of phenotypes that have been overlayed
                if iscell(app.figOpts.defaults.(num).P)
                    N_phn = numel(app.figOpts.defaults.(num).P);
                else
                    app.figOpts.defaults.(num).P = {app.figOpts.defaults.(num).P};
                    N_phn = 1;
                end

                [datALL, vrnms, ~] = Plt_Helper.UpdateFigMenu(app, num, 'Combine', false);
                app.figOpts.dat.(num) = datALL;
                app.figOpts.vrnms.(num) = vrnms;
                if iscell(app.figOpts.defaults.(num).P)
                    % Pull the selected sample name
                    app.figOpts.phenos.(num) = app.figOpts.defaults.(num).P;
                end
                if app.figOpts.defaults.(num).pchnd || app.figOpts.defaults.(num).schnd || app.figOpts.defaults.(num).sldrchnd

                    % Update Slider ranges if the axis changed.
                    maxZ = 0;
                    minZ = 0;
                    for smp_i=1:N_smp
                        for phn_i=1:N_phn
                            Slider = datALL{smp_i, phn_i}.(vrnms{app.figOpts.sldrIF.(num).Value});
                            if isempty(Slider)
                                Slider = [-10, 10];
                            end
                            if smp_i == 1 && phn_i == 1
                                maxZ = max(Slider);
                                minZ = min(Slider);
                            else
                                maxZ = max([maxZ max(Slider)]);
                                minZ = min([minZ min(Slider)]);
                            end
                            if maxZ - minZ == 0
                                minZ = -1.;
                                maxZ = 1.;
                            end
                        end
                    end
                    app.figOpts.SldrZ.(num).Min = minZ-0.1*minZ;
                    app.figOpts.SldrZ.(num).Max = maxZ+0.1*maxZ;
                    app.figOpts.sldrRNG.(num).String = num2str(app.figOpts.SldrZ.(num).Max-app.figOpts.SldrZ.(num).Min);
                    app.figOpts.SldrZ.(num).Visible = 'on';
                    RNG = str2double(app.figOpts.sldrRNG.(num).String);
                    app.figOpts.SldrTXT.(num).String = 'StopListner';
                    app.figOpts.SldrZ.(num).Value = app.figOpts.SldrZ.(num).Min+RNG/2;
                end

                RNG = str2double(app.figOpts.sldrRNG.(num).String);
                Slider_VAL = app.figOpts.SldrZ.(num).Value;
                app.figOpts.SldrTXT.(num).Position(1:2) = app.figOpts.SldrTXT.(num).Position(1:2) + ...
                    app.figOpts.SldrTXT.(num).Position(3:4) + ...
                    [Constants.plt_menu_x0, Constants.plt_menu_y0];
                full_name = Helper.full(app, vrnms{app.figOpts.sldrIF.(num).Value}, app.figOpts.defaults.(num).S(1));
                full_name = full_name{1};
                app.figOpts.SldrTXT.(num).String = ['Slidebar: ' full_name ' ' num2str(Slider_VAL-RNG/2, 3) ' to ' num2str(Slider_VAL+RNG/2, 3)];
                app.figOpts.SldrTXT.(num).Position(3:4) = app.figOpts.SldrTXT.(num).Extent(3:4) + [30 0];
                app.figOpts.SldrTXT.(num).Position(1:2) = app.figOpts.SldrTXT.(num).Position(1:2) - ...
                    app.figOpts.SldrTXT.(num).Position(3:4) - ...
                    [Constants.plt_menu_x0, Constants.plt_menu_y0];
                
            else % if you are only plotting one sample/phenotype

                cla
                N_smp = 1;
                N_phn = 1;
                app.figOpts.defaults.(num).P = app.figOpts.phenos.(num);
                app.figOpts.defaults.(num).S =  app.figOpts.smpls.(num);

                % Update the drop down menus with the current phenotypes
                [datALL, vrnms, ~] = Plt_Helper.UpdateFigMenu(app, num, 'Combine', true);
                app.figOpts.dat.(num) = datALL;
                app.figOpts.vrnms.(num) = vrnms;

                % If you changed the slider-axis
                if app.figOpts.defaults.(num).sldrchnd 

                    Slider = datALL.(vrnms{app.figOpts.sldrIF.(num).Value});
                    % Only really want to define the range if the
                    % slider axis has been changed
                    if max(Slider) - min(Slider) == 0
                        minZ = min(Slider) - 1.;
                        maxZ = max(Slider) + 1.;
                    else
                        minZ = min(Slider);
                        maxZ = max(Slider);
                    end
                    app.figOpts.SldrZ.(num).Min = minZ-0.1*minZ;
                    app.figOpts.SldrZ.(num).Max = maxZ+0.1*maxZ;

                    app.figOpts.sldrRNG.(num).String = num2str(app.figOpts.SldrZ.(num).Max-app.figOpts.SldrZ.(num).Min);

                    app.figOpts.SldrZ.(num).Visible = 'on';
                    RNG = str2double(app.figOpts.sldrRNG.(num).String);
                    % Make it so setting this doesn't trigger the listner function
                    app.figOpts.SldrTXT.(num).String = 'StopListner';
                    app.figOpts.SldrZ.(num).Value = app.figOpts.SldrZ.(num).Min+RNG/2;
                    Slider_VAL = app.figOpts.SldrZ.(num).Value;

                    % Update the label of the Slider, and re-center it.
                    app.figOpts.SldrTXT.(num).Position(1:2) = app.figOpts.SldrTXT.(num).Position(1:2) + ...
                        app.figOpts.SldrTXT.(num).Position(3:4) + ...
                        [Constants.plt_menu_x0, Constants.plt_menu_y0];
                    full_name = Helper.full(app, vrnms{app.figOpts.sldrIF.(num).Value}, app.figOpts.defaults.(num).S);
                    full_name = full_name{1};
                    app.figOpts.SldrTXT.(num).String = ['Slidebar: ' full_name ' ' num2str(Slider_VAL-RNG/2, 3) ' to ' num2str(Slider_VAL+RNG/2, 3)];
                    app.figOpts.SldrTXT.(num).Position(3:4) = app.figOpts.SldrTXT.(num).Extent(3:4) + [30 0];
                    app.figOpts.SldrTXT.(num).Position(1:2) = app.figOpts.SldrTXT.(num).Position(1:2) - ...
                        app.figOpts.SldrTXT.(num).Position(3:4) - ...
                        [Constants.plt_menu_x0, Constants.plt_menu_y0];
                end
            end

            %% Pull the selected phenotype, if you change phenotype keep axes
            if app.figOpts.defaults.(num).pchnd || app.figOpts.defaults.(num).schnd
                % Make default x, y axes
                ind_i = find(strcmp(app.figOpts.xaxIF.(num).String, app.figOpts.xaxIFTMP.(num)));
                if ~isempty(ind_i)
                    app.figOpts.xaxIF.(num).Value = ind_i;
                else
                    app.figOpts.xaxIF.(num).Value = find(strcmp(app.figOpts.xaxIF.(num).String,...
                        app.figOpts.defaults.(num).X));
                end
                
                ind_i = find(strcmp(app.figOpts.yaxIF.(num).String, app.figOpts.yaxIFTMP.(num)));
                if ~isempty(ind_i)
                    app.figOpts.yaxIF.(num).Value = ind_i;
                else
                    app.figOpts.yaxIF.(num).Value = find(strcmp(app.figOpts.yaxIF.(num).String, ...
                        app.figOpts.defaults.(num).Y));
                end
                
                ind_i = find(strcmp(app.figOpts.zaxIF.(num).String, app.figOpts.zaxIFTMP.(num)));
                if ~isempty(ind_i)
                    app.figOpts.zaxIF.(num).Value = ind_i;
                else
                    app.figOpts.zaxIF.(num).Value = find(strcmp(app.figOpts.zaxIF.(num).String, ...
                        app.figOpts.defaults.(num).Z));
                end

                ind_i = find(strcmp(app.figOpts.caxIF.(num).String, app.figOpts.caxIFTMP.(num)));
                if ~isempty(ind_i)
                    % Sometimes ind_i has multiple elements
                    if numel(ind_i) > 1
                        ind_i = ind_i(1);
                    end
                    app.figOpts.caxIF.(num).Value = ind_i;
                else
                    app.figOpts.caxIF.(num).Value = find(strcmp(app.figOpts.caxIF.(num).String, ...
                        app.figOpts.defaults.(num).C));
                end

            end

            app.figOpts.xaxIFTMP.(num) = app.figOpts.xaxIF.(num).String{app.figOpts.xaxIF.(num).Value};
            app.figOpts.yaxIFTMP.(num) = app.figOpts.yaxIF.(num).String{app.figOpts.yaxIF.(num).Value};
            app.figOpts.zaxIFTMP.(num) = app.figOpts.zaxIF.(num).String{app.figOpts.zaxIF.(num).Value};
            app.figOpts.caxIFTMP.(num) = app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value};
            app.figOpts.sldrIFTMP.(num)= app.figOpts.sldrIF.(num).String{app.figOpts.sldrIF.(num).Value};

            app.figOpts.phnTMP.(num) = app.figOpts.phenos.(num);
            app.figOpts.sampleTMP.(num) = app.figOpts.smpls.(num);

            % Find the index shift for the color axis data list
            IndOff = numel(app.figOpts.caxIF.(num).String) - numel(vrnms);

            %% Plot Things
            plttypeCAX = app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value};
            plttypeYAX = app.figOpts.yaxIF.(num).String{app.figOpts.yaxIF.(num).Value};
            
            if strcmp(plttypeCAX, 'Point Density Grid') ...
                    || strcmp(plttypeCAX, 'Point Density Contour') ...
                    || strcmp(plttypeCAX, 'Point Density')...
                    && ~strcmp(plttypeYAX, 'HistogramNum') && ~strcmp(plttypeYAX, 'HistogramFreq')
                %% Plot point density
                gran = 100;
                [dat, vrnms, ~] = Plt_Helper.UpdateFigMenu(app, num, 'Combine', true);
                app.figOpts.dat.(num) = dat;
                app.figOpts.vrnms.(num) = vrnms;
                X = dat.(vrnms{app.figOpts.xaxIF.(num).Value});
                Y = dat.(vrnms{app.figOpts.yaxIF.(num).Value});
                
                if strcmp(plttypeCAX, 'Point Density Grid')
                    C = hist3([Y,X], 'Nbins', [gran, gran]);
                    % Get plotting coordiates
                    if max(X) == max(X)
                        if max(X) == 0
                            X_plot = linspace(-0.1, 0.1, size(C, 2));
                        else
                            X_plot = linspace(min(X), max(X), size(C, 2));
                        end
                    else
                        X_plot = linspace(min(X), max(X), size(C, 2));
                    end
                    if max(Y) == max(Y)
                        if max(Y) == 0
                            Y_plot = linspace(-0.1, 0.1, size(C, 1));
                        else
                            Y_plot = linspace(min(Y), max(Y), size(C, 1));
                        end
                    else
                        Y_plot = linspace(min(Y), max(Y), size(C, 1));
                    end
    %                 pcolor(X_plot, Y_plot, C);
                    imagesc(X_plot, Y_plot, C);
                    set(gca, 'YDir', 'normal')
                    view(0, 90)
                    cbar = colorbar;
                    colormap(jet)
                    ylabel(cbar, 'Points per bin');
                    xlabel(app.figOpts.xaxIF.(num).String{app.figOpts.xaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                    ylabel(app.figOpts.yaxIF.(num).String{app.figOpts.yaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                    if max(C, [], 'all') ~= min(C, [], 'all')
                        ax.Colorbar.Limits = [min(C, [], 'all'), max(C, [], 'all')];
                        app.figOpts.limits.(num)(4, 1:2) = {min(C, [], 'all'), max(C, [], 'all')};
                    end
                    caxis(ax.Colorbar.Limits);
                elseif strcmp(plttypeCAX, 'Point Density Contour') || strcmp(plttypeCAX, 'Point Density')
                    %% Make a plot with a color axis that is point density with contours
                    Z = dat.(vrnms{app.figOpts.zaxIF.(num).Value});
                    Slider = dat.(vrnms{app.figOpts.sldrIF.(num).Value});
                    Slider_VAL = app.figOpts.SldrZ.(num).Value;
                    RNG = str2double(app.figOpts.sldrRNG.(num).String);
                    app.figOpts.SldrTXT.(num).String = ['Slidebar: ' vrnms{app.figOpts.sldrIF.(num).Value} ' ' num2str(Slider_VAL-RNG/2, 3) ' to ' num2str(Slider_VAL+RNG/2, 3)];
                    IND = (Slider_VAL-RNG/2)<=Slider & Slider<=(Slider_VAL+RNG/2);
                    if sum(IND) == 0
                        return;
                    end
                    Z = Z(IND);
                    X = X(IND);
                    Y = Y(IND);
                    if app.figOpts.D3.(num)
                        ax = gca;
                        plt_PointDensity.PlotDensity3D(ax, [X, Y, Z], 300);
% % %                                 plt = scatter3(X, Y, Z, 5, C, 'filled', 'DisplayName', dsnm);
% % %                                 zlabel(app.figOpts.zaxIF.(num).String{app.figOpts.zaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                    else
                        ax = gca;
                        plt_PointDensity.Draw(ax, [X, Y], strcmp(plttypeCAX, 'Point Density Contour'));
                    end
% %                     plt.MarkerEdgeColor = 'none';
                    cbar = colorbar;
% %                             cbar.Color = app.GUIOPTS.txtclr;
                    % Add  smoothed line
                    switch app.PlotMenu.(num).Smooth_Gauss.Checked
                        case 'on'
                        hold on
                        [X, IND_X] = sort(X);
                        plot(X, smoothdata(Y(IND_X),'gaussian', 100), '.', 'DisplayName', 'Smoothed' );
                    end
                    switch app.PlotMenu.(num).Smooth_mean.Checked
                        case 'on'
                        hold on
                        [X, IND_X] = sort(X);
                        plot(X, movmean(Y(IND_X), 100), '.', 'DisplayName', 'Smoothed' );
                    end
                    switch app.PlotMenu.(num).Smooth_CumSum.Checked
                        case 'on'
                        hold on
                        [X, IND_X] = sort(X);
                        plot(X, max(Y).*(cumsum(Y(IND_X))/max(cumsum(Y(IND_X)))), '.-', 'DisplayName', 'Smoothed' );
                    end
                    ylabel(cbar, 'Points per bin');
                    xlabel(app.figOpts.xaxIF.(num).String{app.figOpts.xaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                    ylabel(app.figOpts.yaxIF.(num).String{app.figOpts.yaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                end % end type of density plot conditional
            else
                min_c = 0;
                max_c = 1;
                num_ticks = 1;
                for smp_i=1:N_smp
                    for phn_i=1:N_phn
                        % Pull out the sample that is to be plotted
                        if iscell(datALL)
                            dat = datALL{smp_i, phn_i};
                        else
                            dat = datALL;
                        end
                        % pull out the current phenotype
                        if iscell(app.figOpts.defaults.(num).P)
                            phnm = app.figOpts.defaults.(num).P{phn_i};
                        else
                            phnm = app.figOpts.defaults.(num).P;
                        end
                        % Pull the current sample
                        if iscell(app.figOpts.defaults.(num).S)
                            smpl = app.figOpts.defaults.(num).S{smp_i};
                        else
                            smpl = app.figOpts.defaults.(num).S;
                        end
                        % define the data name for plot title etc
                        dsnm = strcat(smpl, '_', phnm);
                        if iscell(dsnm)
                            dsnm = dsnm{1};
                        end

                        %% Decide the type of plot and display it
                        bins = 100; % for histograms
                        if strcmp(plttypeCAX, 'None') || strcmp(plttypeYAX, 'HistogramNum') || strcmp(plttypeYAX, 'HistogramFreq') % If there is no color axis
                            if strcmp(plttypeYAX, 'HistogramNum') % If user selected histogram
                                X = dat.(vrnms{app.figOpts.xaxIF.(num).Value});

                                % Adjust what data is plotted with the slider
                                Slider = dat.(vrnms{app.figOpts.sldrIF.(num).Value});
                                Slider_VAL = app.figOpts.SldrZ.(num).Value;
                                RNG = str2double(app.figOpts.sldrRNG.(num).String);
                                full_name = Helper.full(app, vrnms{app.figOpts.sldrIF.(num).Value}, app.figOpts.defaults.(num).S);
                                full_name = full_name{1};
                                app.figOpts.SldrTXT.(num).String = ['Slidebar: ' full_name ' ' num2str(Slider_VAL-RNG/2, 3) ' to ' num2str(Slider_VAL+RNG/2, 3)];
                                IND = (Slider_VAL-RNG/2)<=Slider & Slider<=(Slider_VAL+RNG/2);
                                X = X(IND);
                                % if there are not points in the selected Z
                                % range, don't plot anything
                                if sum(IND) == 0
                                    return;
                                end
                                %% Custom Histogram
                                [Y, N] = histcounts(X, bins);
                                Centers = N+abs((N(2)-N(1))/2);
                                Centers = Centers(1:end-1);
                                hold on
    %                             plot(Centers, Y, '.k', 'DisplayName', dsnm)
                                % Make a splin interpolant function
                                % Find coefficients for spline interpolant
                                splinePara = spline(Centers,Y);
                                % Evaluate piecewise polynomial
                                xplot = min(Centers):((max(Centers)-min(Centers))/(4*bins)):max(Centers);
                                splineDat = ppval(splinePara,xplot);
                                % Plot the fit
                                plot(xplot,splineDat,'-','DisplayName', [dsnm ', Fit'], 'LineWidth', 1);
    %                             area(xplot,splineDat, 'FaceAlpha', 0.5, 'DisplayName', [dsnm ', Fit'])
                                ylabel('Number of cells', 'Color', app.GUIOPTS.txtclr, 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                                xlabel(app.figOpts.xaxIF.(num).String{app.figOpts.xaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                            elseif strcmp(plttypeYAX, 'HistogramFreq') % If user selected histogram
                                %% Custom Histogram
                                X = dat.(vrnms{app.figOpts.xaxIF.(num).Value});

                                % Adjust what data is plotted with the slider
                                Slider = dat.(vrnms{app.figOpts.sldrIF.(num).Value});
                                Slider_VAL = app.figOpts.SldrZ.(num).Value;
                                RNG = str2double(app.figOpts.sldrRNG.(num).String);
                                full_name = Helper.full(app, vrnms{app.figOpts.sldrIF.(num).Value}, app.figOpts.defaults.(num).S);
                                full_name = full_name{1};
                                app.figOpts.SldrTXT.(num).String = ['Slidebar: ' full_name ' ' num2str(Slider_VAL-RNG/2, 3) ' to ' num2str(Slider_VAL+RNG/2, 3)];
                                IND = (Slider_VAL-RNG/2)<=Slider & Slider<=(Slider_VAL+RNG/2);
                                if sum(IND) == 0
                                    return;
                                end
                                X = X(IND);
                                [Y, N] = histcounts(X, bins);
                                % Normalize the y axis
                                Y = Y./max(Y);
                                Centers = N+abs((N(2)-N(1))/2);
                                Centers = Centers(1:end-1);
                                hold on
                                % Make a splin interpolant function
                                % Find coefficients for spline interpolant
                                splinePara = spline(Centers,Y);
                                % Evaluate piecewise polynomial
                                xplot = min(Centers):((max(Centers)-min(Centers))/(4*bins)):max(Centers);
                                splineDat = ppval(splinePara,xplot);
                                % Plot the fit
                                plot(xplot,splineDat,'-','DisplayName', [dsnm ', Fit'], 'LineWidth', 1);
                                ylabel('Percentage of cells', 'Color', app.GUIOPTS.txtclr, 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                                xlabel(app.figOpts.xaxIF.(num).String{app.figOpts.xaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                            else
                                %% Normal no-color axis case
                                X = dat.(vrnms{app.figOpts.xaxIF.(num).Value});
                                Y = dat.(vrnms{app.figOpts.yaxIF.(num).Value});
                                % Pull the Z axis
                                Z = dat.(vrnms{app.figOpts.zaxIF.(num).Value});
                                Slider = dat.(vrnms{app.figOpts.sldrIF.(num).Value});
                                Slider_VAL = app.figOpts.SldrZ.(num).Value;
                                RNG = str2double(app.figOpts.sldrRNG.(num).String);
                                full_name = Helper.full(app, vrnms{app.figOpts.sldrIF.(num).Value}, app.figOpts.defaults.(num).S{1});
                                full_name = full_name{1};
                                app.figOpts.SldrTXT.(num).String = ['Slidebar: ' full_name ' ' num2str(Slider_VAL-RNG/2, 3) ' to ' num2str(Slider_VAL+RNG/2, 3)];
                                IND = (Slider_VAL-RNG/2) <= Slider & Slider <= (Slider_VAL+RNG/2);
                                % If scope of Z axis is no cells then plot nothing.
                                if sum(IND) == 0
                                    return;
                                end

                                X = X(IND);
                                Y = Y(IND);
                                Z = Z(IND);

                                % Display only the
                                if app.figOpts.D3.(num)
%                                     plt = plot3(X, Y, Z, '.', 'DisplayName', dsnm);
                                    plt = scatter3(X, Y, Z, 'o', 'filled', 'DisplayName', dsnm);
                                    zlabel(app.figOpts.zaxIF.(num).String{app.figOpts.zaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                                else
%                                     plt = plot(X, Y, '.', 'DisplayName', dsnm);
                                    plt = scatter(X, Y, 'o', 'filled', 'DisplayName', dsnm);
                                    hold on
                                    % view(0,90) % could also use view(2) = default 2D view
                                end
                                plt.SizeData = 5;
%                                 plt.MarkerSize = 5;
                                ylabel(app.figOpts.yaxIF.(num).String{app.figOpts.yaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                                xlabel(app.figOpts.xaxIF.(num).String{app.figOpts.xaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                                
                                switch app.PlotMenu.(num).Smooth_Gauss.Checked
                                    case 'on'
                                    hold on
                                    [X, IND_X] = sort(X);
%                                     plot(smoothdata(X,'gaussian', 100), smoothdata(Y(IND_X),'gaussian', 100), '.', 'DisplayName', ['Smoothed_' dsnm]);
                                    plot(X, smoothdata(Y(IND_X),'gaussian', 100), '.', 'DisplayName', ['Smoothed_' dsnm]);
                                end
                                
                                switch app.PlotMenu.(num).Smooth_mean.Checked
                                    case 'on'
                                    hold on
                                    [X, IND_X] = sort(X);
                                    plot(X, movmean(Y(IND_X), 100), '.', 'DisplayName', ['Smoothed_' dsnm]);
                                end
                                
                                switch app.PlotMenu.(num).Smooth_CumSum.Checked
                                    case 'on'
                                    hold on
                                    [X, IND_X] = sort(X);
                                    plot(X, max(Y).*(cumsum(Y(IND_X))/max(cumsum(Y(IND_X)))), '.-', 'DisplayName', ['Smoothed_' dsnm]);
                                end
                            end
                        % If the user selected some color axis and not a histogram
                        elseif strcmp(plttypeCAX, 'Number of cells / Neighborhood')
                            %% Plot cell densities from MFIRSN table
                            if isfield(app.data.(app.DataN.Value), 'MFIRSN')
                                X = app.data.(smpl).MFIRSN.X;
                                Y = app.data.(smpl).MFIRSN.Y;
                                Z = app.data.(smpl).MFIRSN.Z;
                                if numel(app.figOpts.phenos.(num)) == 1 && strcmp(app.figOpts.phenos.(num)(1), 'All Cells')
                                    C = app.data.(smpl).MFIRSN.('NCells');
                                elseif numel(app.figOpts.phenos.(num)) == 1 && strcmp(app.figOpts.phenos.(num)(1), 'Density/MFI RSN')
                                    C = app.data.(smpl).MFIRSN.('NCells');
                                else
                                    tags = Helper.get_tag(app, app.figOpts.phenos.(num), smpl);
                                    if all(~ismember(tags, app.data.(smpl).MFIRSN.Properties.VariableNames))
                                        return;
                                    else
                                        tags = tags(ismember(tags, app.data.(smpl).MFIRSN.Properties.VariableNames));
                                    end
                                    C = app.data.(smpl).MFIRSN(:, tags);
                                    C = sum(table2array(C), 2);  % Sum all the columns, to have all cells counts together.
                                end

                                % 3D implementation of heatmap
                                C = reshape(C, [numel(unique(Y)), numel(unique(X)), numel(unique(Z))]);
                                % Do sum along z
                                C = sum(C, 3);
                                % Select specific Z plane
                                imagesc('XData', unique(X), 'YData', unique(Y), 'CData', C, [0 max(max(C))])
                                cbar = colorbar;
                                colormap(app.map.jet)
                                ylabel(cbar, 'Cells / Neighborhood')
                                ylabel('Y', 'Color', app.GUIOPTS.txtclr, 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                                xlabel('X', 'Color', app.GUIOPTS.txtclr, 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                                title(dsnm, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                                legend off
                                colorbar;
                                break;  % We are doing all phenos in one go, so we can skip rest of for loop.
                            else
                                % If you haven't run cellularity yet
                                errordlg(sprintf('I"m afraid I can"t do that. Please run Raster Scan Neighborhoods first'),'Error');
                            end
                        elseif strcmp(plttypeCAX, 'Point Density Contour') || strcmp(plttypeCAX, 'Point Density')
                            % Do nothing here
                        else % If there is some other color axis
                                %% Make a plot with a color axis
                                X = dat.(vrnms{app.figOpts.xaxIF.(num).Value});
                                Y = dat.(vrnms{app.figOpts.yaxIF.(num).Value});
                                C = dat.(vrnms{app.figOpts.caxIF.(num).Value-IndOff});
                                Z = dat.(vrnms{app.figOpts.zaxIF.(num).Value});
                                Slider = dat.(vrnms{app.figOpts.sldrIF.(num).Value});
                                Slider_VAL = app.figOpts.SldrZ.(num).Value;
                                RNG = str2double(app.figOpts.sldrRNG.(num).String);
                                app.figOpts.SldrTXT.(num).String = ['Slidebar: ' vrnms{app.figOpts.sldrIF.(num).Value} ' ' num2str(Slider_VAL-RNG/2, 3) ' to ' num2str(Slider_VAL+RNG/2, 3)];
                                IND = (Slider_VAL-RNG/2)<=Slider & Slider<=(Slider_VAL+RNG/2);
                                if sum(IND) == 0
                                    return;
                                end
                                Z = Z(IND);
                                X = X(IND);
                                Y = Y(IND);
                                C = C(IND);

                                if app.figOpts.D3.(num)
                                    plt = scatter3(X, Y, Z, 5, C, 'filled', 'DisplayName', dsnm);
                                    zlabel(app.figOpts.zaxIF.(num).String{app.figOpts.zaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                                else
                                    plt = scatter(X, Y, 5, C, 'filled', 'DisplayName', dsnm);
                                end
                                plt.MarkerEdgeColor = 'none';
                                cbar = colorbar;
                                cbar.Color = app.GUIOPTS.txtclr;

                                % If you selected a clustering model colorscale
                                if ismember(...
                                        app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, ...
                                        Helper.full_var(ModelNames))
                                    % This is never run, since the
                                    % Well... Now it runs
                                    % IGNORE
                                else
                                    min_c = min(min_c, round(min(min(C))));
                                    max_c = max(max_c, round(max(max(C))));
                                    num_ticks = max(num_ticks, numel(unique(C)));
                                end
                                % Add  smoothed line
                                switch app.PlotMenu.(num).Smooth_Gauss.Checked
                                    case 'on'
                                    hold on
                                    [X, IND_X] = sort(X);
%                                     plot(smoothdata(X,'gaussian', 100), smoothdata(Y(IND_X),'gaussian', 100), '.', 'DisplayName', ['Smoothed_' dsnm]);
                                    plot(X, smoothdata(Y(IND_X),'gaussian', 100), '.', 'DisplayName', ['Smoothed_' dsnm]);
                                end
                                
                                switch app.PlotMenu.(num).Smooth_mean.Checked
                                    case 'on'
                                    hold on
                                    [X, IND_X] = sort(X);
                                    plot(X, movmean(Y(IND_X), 100), '.', 'DisplayName', ['Smoothed_' dsnm]);
                                end
                                
                                switch app.PlotMenu.(num).Smooth_CumSum.Checked
                                    case 'on'
                                    hold on
                                    [X, IND_X] = sort(X);
                                    plot(X, max(Y).*(cumsum(Y(IND_X))/max(cumsum(Y(IND_X)))), '.-', 'DisplayName', ['Smoothed_' dsnm]);
                                end
                        end
                    end % end of loop through all phenotypes
                end % end of loop through all samples
                app.figOpts.plt_ax.(num) = ax;
                if strcmp(plttypeCAX, 'None') && ~ismember(plttypeYAX, {'HistogramFreq', 'HistogramNum'})
                    if isfield(app.figOpts.order_stack, num)
                        Plt_Helper.update_plot_order(app, num);

                        ax = gca;
                        colors = zeros(numel(ax.Children), 3);
                        for k_idx=1:numel(ax.Children)
                            if ismember('Color',fieldnames(ax.Children(k_idx)))
                                colors(k_idx, :) = ax.Children(k_idx).Color;
                            elseif ismember('CData',fieldnames(ax.Children(k_idx)))
                                colors(k_idx, :) = ax.Children(k_idx).CData;
                            end
                        end
                        ax.Children = ax.Children(numel(app.figOpts.order_stack.(num)) + 1 - app.figOpts.order_stack.(num));
                        for k_idx=1:numel(ax.Children)
                            if ismember('Color',fieldnames(ax.Children(k_idx)))
                                ax.Children(k_idx).Color = colors(k_idx, :);
                            elseif ismember('CData',fieldnames(ax.Children(k_idx)))
                                ax.Children(k_idx).CData = colors(k_idx, :);
                            end
                        end
                    end
                end
            end % end of all plotting

            
            % Set Colorbar
            if ~ismember(plttypeCAX, {'None', 'Number of cells / Neighborhood', 'Point Density Grid', 'Point Density Contour', 'Point Density'}) && ...
                    ~ismember(plttypeYAX, {'HistogramNum', 'HistogramFreq'})
                if ismember(...
                        app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, ...
                        Helper.full_var(ModelNames))
                    INDModel = strcmp( ...
                        Helper.valid_var(app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}), ...
                        ModelNames ...
                    );
                    ModelName = ModelNames{INDModel};
                    cbar.Limits = [0 app.net.(ModelName).NR];
                    cbar.Ticks = app.net.(ModelName).NR;
                    % Change the color scale to match the limits
                    caxis([cbar.Limits(1), cbar.Limits(2)])
                    colormap(app.net.(ModelName).cmap)
                    
                else
                    if min_c == max_c && min_c == 0
                        min_c = -0.1;
                        max_c = 0.1;
                    elseif min_c == max_c
                        min_c = 0.9 * min_c;
                        max_c = 1.1 * max_c;
                    end
                    max_ticks = 100;
                    num_ticks = min(max_ticks, num_ticks);
                    cbar.Ticks = linspace(min_c, max_c, num_ticks);
                    cbar.Limits = [min_c, max_c];
                    caxis(cbar.Limits);
                    app.map.jet = colormap(jet(numel(cbar.Ticks)));
                    colormap(jet);
                    if size(jet, 1)==1
                        colormap(jet(2));
                    end
                end
                % clabel(cbar, app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value})
                ylabel(app.figOpts.yaxIF.(num).String{app.figOpts.yaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                xlabel(app.figOpts.xaxIF.(num).String{app.figOpts.xaxIF.(num).Value}, 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                title(strjoin(app.figOpts.phenos.(num), " "), 'Color', app.figOpts.txtclr.(num), 'FontSize', app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI})
                legend off
                colorbar;

            end
            %% Set some plot properties
            if app.figOpts.eqax.(num)
                axis equal
            else
                axis normal
            end
            ax = gca;
            axis tight
            ax.Color = app.figOpts.bgclr.(num);
            ax.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
            ax.YColor = app.figOpts.txtclr.(num);
            ax.XColor = app.figOpts.txtclr.(num);
            if app.figOpts.D3.(num)
                ax.ZColor = app.figOpts.txtclr.(num);
            end
            % Set the viewing angle
            ax.View = app.figOpts.View.(num);
            ax.Box = 'off';
            ax.XGrid = 'off';
            ax.YGrid = 'off';
            ax.ZGrid = 'off';

            % Set the plot axes
            % if the x axis changed
            if app.figOpts.defaults.(num).xchnd || app.figOpts.defaults.(num).pchnd || app.figOpts.defaults.(num).schnd
                app.figOpts.limits.(num)(1,1) = {ax.XLim(1)};
                app.figOpts.limits.(num)(1,2) = {ax.XLim(2)};
            else
                ax.XLim = [app.figOpts.limits.(num){1,1:2}];
            end
            % if the y axis changed
            if app.figOpts.defaults.(num).ychnd || app.figOpts.defaults.(num).pchnd || app.figOpts.defaults.(num).schnd
                app.figOpts.limits.(num)(2,1) = {ax.YLim(1)};
                app.figOpts.limits.(num)(2,2) = {ax.YLim(2)};
            else
                ax.YLim = [app.figOpts.limits.(num){2,1:2}];
            end
            % If the Z axis changed
            if size(app.figOpts.limits.(num), 1)>=3
                if app.figOpts.defaults.(num).zchnd || app.figOpts.defaults.(num).sldrchnd ||app.figOpts.defaults.(num).pchnd || app.figOpts.defaults.(num).schnd && app.figOpts.D3.(num)
                    app.figOpts.limits.(num)(3,1) = {ax.ZLim(1)};
                    app.figOpts.limits.(num)(3,2) = {ax.ZLim(2)};
                elseif ~any(cellfun(@isempty, app.figOpts.limits.(num)(3,:)))
                    ax.ZLim = [app.figOpts.limits.(num){3,1:2}];
                end
            end
            % if the c axis changed
            if app.figOpts.defaults.(num).cchnd ...
                    || app.figOpts.defaults.(num).pchnd ...
                    || app.figOpts.defaults.(num).schnd ...
                    && ~strcmp(app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, 'None') ...
                    && ~ismember(...
                            app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, ...
                            Helper.full_var(ModelNames))
                cbar = colorbar;
                limits = ax.Colorbar.Limits;
                caxis(limits);
                ax.Colorbar.Limits = limits;
                ax.Colorbar.Ticks = limits(1):((limits(2)-limits(1))/10):limits(2);

                app.figOpts.limits.(num)(4,1) = {ax.Colorbar.Limits(1)};
                app.figOpts.limits.(num)(4,2) = {ax.Colorbar.Limits(2)};
                colormap(app.map.jet)
                if size(app.map.jet, 1)==1
                    app.map.jet = jet(2);
                    colormap(app.map.jet);
                end
                ylabel(cbar, app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value})
                colorbar;
            elseif ~strcmp(app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, 'None') ...
                    && ~ismember(...
                            app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, ...
                            Helper.full_var(ModelNames))
                cbar = colorbar;
                limits = [app.figOpts.limits.(num){4,1:2}];
                caxis(limits);
                ax.Colorbar.Limits = limits;
                ax.Colorbar.Ticks = limits(1):((limits(2)-limits(1))/10):limits(2);
                colormap;
                ylabel(cbar, app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value})
                colorbar;
            end
            if strcmp(app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, 'None')
                colorbar off
            end

            if ismember(...
                    app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}, ...
                    Helper.full_var(ModelNames))
                INDModel = strcmp(Helper.valid_var(app.figOpts.caxIF.(num).String{app.figOpts.caxIF.(num).Value}), ModelNames);
                ModelName = ModelNames{INDModel};
                ax.Colorbar.Limits = [0 app.net.(ModelName).NR];
                ax.Colorbar.Ticks = round(ax.Colorbar.Limits(1):ax.Colorbar.Limits(2));
                % Change the color scale to match the limits
                caxis([ax.Colorbar.Limits(1), ax.Colorbar.Limits(2)])
                colormap(app.net.(ModelName).cmap)
                colorbar;
            end
            
            % if there was an image overlay, add it back in
            if overlay == 1
                if exist('imgax', 'var')
                    % replot the image
                    image(imgax.XData, imgax.YData, imgax.CData);
                    % Put the image on bottom
                    chH = get(gca,'Children');
                    set(gca,'Children',[chH(2:end);chH(1)]);
                end
            end
% %             %%%%%%%%%%%%%% I don't understand why this is needed
% %             % flip the plot order for some reason??
% %             revorder = numel(ax.Children):-1:1;
% %             ax.Children = ax.Children(revorder);
% %             %%%%%%%%%%%%%%
            
            app.figOpts.plt_ax.(num) = ax;
        end

        function resize_figure(src, app)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % resize_figure allows to resize and scale the whole plot and
            % figure alongside with all of the buttons and drop downs to
            % fit the new size of the figure.
            %
            % Input:
            %   - src - Handle to UIFigure
            %   - app - Instance of CytoMAP
            %
            % Modifies:
            %   - app - More specififcally fields in app.figOpts, by
            %       setting Positions to new, scaled stuff.
            %   - src - If size of src is smaller then something set in
            %       CytoMAP, then it's resize to at least that size.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fig = gcf;
            num = ['fig' num2str(fig.Number)];
            if ~strcmp(fig.Name, 'RSZ')
                return;
            end
            %% Settings for resizing
            delta_x = Constants.plt_menu_dx;  % Distance between 2 things horizontally
            delta_y = Constants.plt_menu_dy;  % Distance between 2 things vertically
            x_0 = Constants.plt_menu_x0;  % Initial displacement in x from corner
            y_0 = Constants.plt_menu_y0;  % Initial displacement in y from corner
            min_fig_size = Constants.min_fig_size;

            %% Set minimal size of figure
            if ~strcmp(src.WindowStyle, 'docked')  % Cannot resize docked stuff
                src.Position(3:4) = max(min_fig_size, src.Position(3:4));
            end

            %% 1st Column - Options and (Table Button?)
            app.figOpts.opts.(num).Position(1) = x_0;
            app.figOpts.opts.(num).Position(2) = y_0;
            app.figOpts.opts.(num).Position(3) = min([src.Position(3) * 0.1, app.figOpts.opts.(num).Extent(3)]);
            app.figOpts.opts.(num).Position(4) = app.figOpts.opts.(num).Extent(4);

            %% Width of each column in columns 2 through 5
            % It is width of both text and dropdown
            width = src.Position(3) - app.figOpts.opts.(num).Position(1) - app.figOpts.opts.(num).Position(3);
            width = width / 3;  % Length per column
            width = max(1, width);  % Make sure it ain't negative

            %% 2nd Column - X and Y
            % Do texts in front of the dropdowns first
            app.figOpts.yaxTXT.(num).Position(1) = app.figOpts.opts.(num).Position(3) + ...
                app.figOpts.opts.(num).Position(1) + delta_x;
            app.figOpts.yaxTXT.(num).Position(2) = y_0;
            app.figOpts.yaxTXT.(num).Position(3) = max([app.figOpts.yaxTXT.(num).Position(3), app.figOpts.xaxTXT.(num).Position(3)]);
            app.figOpts.yaxTXT.(num).Position(4) = max([app.figOpts.yaxTXT.(num).Position(4), app.figOpts.yaxIF.(num).Position(4)]);

            app.figOpts.xaxTXT.(num).Position(1) = app.figOpts.yaxTXT.(num).Position(1);
            app.figOpts.xaxTXT.(num).Position(2) = app.figOpts.yaxTXT.(num).Position(2) + ...
                app.figOpts.yaxTXT.(num).Position(4) + delta_y;
            app.figOpts.xaxTXT.(num).Position(3) = app.figOpts.yaxTXT.(num).Position(3);
            app.figOpts.xaxTXT.(num).Position(4) = max([app.figOpts.xaxTXT.(num).Position(4), app.figOpts.xaxIF.(num).Position(4)]);

            % Dropdowns
            app.figOpts.yaxIF.(num).Position(1) = app.figOpts.yaxTXT.(num).Position(1) + ...
                app.figOpts.yaxTXT.(num).Position(3);
            app.figOpts.yaxIF.(num).Position(2) = y_0;
            app.figOpts.yaxIF.(num).Position(3) = max(1, width - app.figOpts.yaxTXT.(num).Position(3) - delta_x);
            app.figOpts.yaxIF.(num).Position(4) = app.figOpts.yaxTXT.(num).Position(4);

            app.figOpts.xaxIF.(num).Position(1) = app.figOpts.yaxIF.(num).Position(1);
            app.figOpts.xaxIF.(num).Position(2) = app.figOpts.xaxTXT.(num).Position(2);
            app.figOpts.xaxIF.(num).Position(3) = app.figOpts.yaxIF.(num).Position(3);
            app.figOpts.xaxIF.(num).Position(4) = app.figOpts.xaxTXT.(num).Position(4);

            %% 3rd Column - Z and Color
            app.figOpts.caxTXT.(num).Position(1) = app.figOpts.yaxIF.(num).Position(3) + ...
                app.figOpts.yaxIF.(num).Position(1) + delta_x;
            app.figOpts.caxTXT.(num).Position(2) = y_0;
            app.figOpts.caxTXT.(num).Position(3) = max([app.figOpts.caxTXT.(num).Position(3), app.figOpts.zaxTXT.(num).Position(3)]);
            app.figOpts.caxTXT.(num).Position(4) = max([app.figOpts.caxTXT.(num).Position(4), app.figOpts.caxIF.(num).Position(4)]);

            app.figOpts.zaxTXT.(num).Position(1) = app.figOpts.caxTXT.(num).Position(1);
            app.figOpts.zaxTXT.(num).Position(2) = app.figOpts.caxTXT.(num).Position(2) + ...
                app.figOpts.caxTXT.(num).Position(4) + delta_y;
            app.figOpts.zaxTXT.(num).Position(3) = app.figOpts.caxTXT.(num).Position(3);
            app.figOpts.zaxTXT.(num).Position(4) = max([app.figOpts.zaxTXT.(num).Position(4), app.figOpts.zaxIF.(num).Position(4)]);

            app.figOpts.caxIF.(num).Position(1) = app.figOpts.caxTXT.(num).Position(1) + ...
                app.figOpts.caxTXT.(num).Position(3);
            app.figOpts.caxIF.(num).Position(2) = app.figOpts.caxTXT.(num).Position(2);
            app.figOpts.caxIF.(num).Position(3) = max(1, width - app.figOpts.caxTXT.(num).Position(3) - delta_x);
            app.figOpts.caxIF.(num).Position(4) = app.figOpts.caxTXT.(num).Position(4);

            app.figOpts.zaxIF.(num).Position(1) = app.figOpts.caxIF.(num).Position(1);
            app.figOpts.zaxIF.(num).Position(2) = app.figOpts.zaxTXT.(num).Position(2);
            app.figOpts.zaxIF.(num).Position(3) = app.figOpts.caxIF.(num).Position(3);
            app.figOpts.zaxIF.(num).Position(4) = app.figOpts.zaxTXT.(num).Position(4);

            %% 4th Column - Slider Axis and Range
            app.figOpts.sldrRNGTXT.(num).Position(1) = app.figOpts.caxIF.(num).Position(3) + ...
            app.figOpts.caxIF.(num).Position(1) + delta_x;
            app.figOpts.sldrRNGTXT.(num).Position(2) = y_0;
            app.figOpts.sldrRNGTXT.(num).Position(3) = max([app.figOpts.sldrIFTXT.(num).Position(3), app.figOpts.sldrRNGTXT.(num).Position(3)]);
            app.figOpts.sldrRNGTXT.(num).Position(4) = max([app.figOpts.sldrRNGTXT.(num).Position(4), app.figOpts.sldrRNG.(num).Position(4)]);

            app.figOpts.sldrIFTXT.(num).Position(1) = app.figOpts.sldrRNGTXT.(num).Position(1);
            app.figOpts.sldrIFTXT.(num).Position(2) = app.figOpts.sldrRNGTXT.(num).Position(2) + ...
                app.figOpts.sldrRNGTXT.(num).Position(4) + delta_y;
            app.figOpts.sldrIFTXT.(num).Position(3) = app.figOpts.sldrRNGTXT.(num).Position(3);
            app.figOpts.sldrIFTXT.(num).Position(4) = max([app.figOpts.sldrIFTXT.(num).Position(4), app.figOpts.sldrIF.(num).Position(4)]);

            app.figOpts.sldrRNG.(num).Position(1) = app.figOpts.sldrRNGTXT.(num).Position(1) + ...
                app.figOpts.sldrRNGTXT.(num).Position(3);
            app.figOpts.sldrRNG.(num).Position(2) = app.figOpts.sldrRNGTXT.(num).Position(2);
            app.figOpts.sldrRNG.(num).Position(3) = max(1, width - app.figOpts.sldrRNGTXT.(num).Position(3) - delta_x);
            app.figOpts.sldrRNG.(num).Position(4) = app.figOpts.sldrRNGTXT.(num).Position(4);

            app.figOpts.sldrIF.(num).Position(1) = app.figOpts.sldrRNG.(num).Position(1);
            app.figOpts.sldrIF.(num).Position(2) = app.figOpts.sldrRNG.(num).Position(2) + app.figOpts.sldrRNG.(num).Position(4) + delta_y;
            app.figOpts.sldrIF.(num).Position(3) = app.figOpts.sldrRNG.(num).Position(3);
            app.figOpts.sldrIF.(num).Position(4) = app.figOpts.sldrRNGTXT.(num).Position(4);

            %% Slider on right top
            app.figOpts.SldrTXT.(num).Position(1) = src.Position(3) - app.figOpts.SldrTXT.(num).Position(3) - x_0;
            app.figOpts.SldrTXT.(num).Position(2) = src.Position(4) - app.figOpts.SldrTXT.(num).Position(4) - y_0;
            app.figOpts.SldrZ.(num).Position(1) = src.Position(3) - app.figOpts.SldrZ.(num).Position(3) - x_0;
            app.figOpts.SldrZ.(num).Position(2) = app.figOpts.sldrIF.(num).Position(2) + app.figOpts.sldrIF.(num).Position(4) + delta_y;
            app.figOpts.SldrZ.(num).Position(4) = app.figOpts.SldrTXT.(num).Position(2) - app.figOpts.SldrZ.(num).Position(2);

            %% Axes and (possibly) Sample/Phenotype Table
            % Create Axes for the plot

            ax = gca;
            % Check what the view is right now.
            if strcmp(app.figOpts.smplPhenoBtn.(num).Visible, 'off')
                % Split screen between table and axes (and close button)
                width_tab_axis = (src.Position(3) - (2.0 * x_0 + app.figOpts.SldrZ.(num).Position(3))) / 2;
                height_tab_axis = src.Position(4) - (y_0 + app.figOpts.sldrIF.(num).Position(2) + app.figOpts.sldrIF.(num).Position(4));

                app.figOpts.phn_smpl_tab.(num).Position = [ ...
                    x_0, ...
                    app.figOpts.sldrIF.(num).Position(2) + app.figOpts.sldrIF.(num).Position(4) + delta_y, ...
                    width_tab_axis, ...
                    height_tab_axis
                ];

                % Set widths of columns inside of the table
                app.figOpts.phn_smpl_tab.(num).ColumnWidth{1} = 55;
                app.figOpts.phn_smpl_tab.(num).ColumnWidth{2} = (app.figOpts.phn_smpl_tab.(num).Position(3) - 110) * 0.6;
                app.figOpts.phn_smpl_tab.(num).ColumnWidth{3} = 55;
                app.figOpts.phn_smpl_tab.(num).ColumnWidth{4} = app.figOpts.phn_smpl_tab.(num).Position(3) - sum(cell2mat(app.figOpts.phn_smpl_tab.(num).ColumnWidth(1:3)));
                
                % Set axis to right half of screen. Offsets are there to keep names visible
                ax.Position(1) = (app.figOpts.phn_smpl_tab.(num).Position(1) + app.figOpts.phn_smpl_tab.(num).Position(3) + 70.0) / src.Position(3);
                ax.Position(2) = (app.figOpts.sldrIF.(num).Position(2) + app.figOpts.sldrIF.(num).Position(4) + delta_y + 40.0) / src.Position(4);
                ax.Position(3) = 1 - ax.Position(1) - ((x_0 + app.figOpts.SldrZ.(num).Position(3) + 90.0) / src.Position(3));
                ax.Position(4) = 1 - ax.Position(2) - ((y_0 + 55.0) / src.Position(4));

                app.figOpts.close_phn_smpl_tab_btn.(num).Position(1) = ...
                    app.figOpts.phn_smpl_tab.(num).Position(1) + ...
                    app.figOpts.phn_smpl_tab.(num).Position(3);
                app.figOpts.close_phn_smpl_tab_btn.(num).Position(2) = ...
                    app.figOpts.phn_smpl_tab.(num).Position(2) + ...
                    app.figOpts.phn_smpl_tab.(num).Position(4) - ...
                    app.figOpts.close_phn_smpl_tab_btn.(num).Position(4);

                app.figOpts.order_list_btn.(num).Position(1) = ...
                    app.figOpts.close_phn_smpl_tab_btn.(num).Position(1);
                app.figOpts.order_list_btn.(num).Position(2) = ...
                    app.figOpts.close_phn_smpl_tab_btn.(num).Position(2) - ...
                    app.figOpts.order_list_btn.(num).Position(4) - delta_y;
            else
                app.figOpts.smplPhenoBtn.(num).Position(1) = x_0;
                app.figOpts.smplPhenoBtn.(num).Position(2) = src.Position(4) - app.figOpts.smplPhenoBtn.(num).Position(4) - y_0;

                % Axes only and show table button
                ax.Position(1) = (x_0 + 70) * 1.0 / src.Position(3);
                ax.Position(2) = (app.figOpts.sldrIF.(num).Position(2) + app.figOpts.sldrIF.(num).Position(4) + delta_y + 50.0) / src.Position(4);
                ax.Position(3) = 1 - ((2.0 * x_0 + 150. + app.figOpts.SldrZ.(num).Position(3)) / src.Position(3));
                ax.Position(4) = 1 - ax.Position(2) - ((y_0 + 35.0) / src.Position(4));
            end

        end

        function close_figure(src, app, num)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % close_figure closes given plot, and removes any trace of it
            % in the CytoMAP
            %
            % Input:
            %   - src - Handle to UIfigure of the plot
            %   - app - Instance of CytoMAP
            %   - num - Figure Identifier (For example 'fig1')
            %
            % Modifies:
            %   - src - deletes it
            %   - app - Removes all app.figOpts.___.(num) fields, where
            %       ___ is every field in figOpts.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            delete(src);
            for f_i=1:numel(fieldnames(app.figOpts))
                f = fieldnames(app.figOpts);
                f = f{f_i};

                if isfield(app.figOpts.(f), num)
                    app.figOpts.(f) = rmfield(app.figOpts.(f), num);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sample/Phenotype Table Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function show_phn_smpl_table(app, num, main_fig, overlay)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % show_phn_smpl_table splits the view in given figure into 2,
            % with Phenotype/Sample table on left, and actual plot on
            % right. This allows user to choose what is plotted in an
            % efficient way
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - num - Figure Identifier (For example 'fig1')
            %   - main_fig - Handle to UIfigure of the plot
            %
            % Modifies:
            %   - app - Modifies, creates fields in app.figOpts
            %       corresponding to table, choices in the table. Also
            %       modifies positions and visibility of few on the plot
            %       which are tracked in the app.figOpts
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get current phenotypes in selected samples if they do not exist
            app.figOpts.smplPhenoBtn.(num).Visible = 'off';
            if ~isfield(app.figOpts.smpls, num) || isempty(app.figOpts.smpls.(num))
                app.figOpts.smpls.(num) = {app.DataN.Value};
            end
                        
            if nargin < 4
                overlay = 0;
            end
            
            % Get current phenotypes in selected samples if they do not exist
            if ~isfield(app.figOpts.phenos, num) || isempty(app.figOpts.phenos.(num))
                smpl = app.figOpts.smpls.(num);
                smpl = smpl{1};
                app.figOpts.phenos.(num) = app.data.(smpl).GateTags(2, 1);
            end

            % Load or make new ui table
            if ~isfield(app.figOpts.phn_smpl_tab, num)
                %% Get initial table and it's properties
                tmp_tab = Helper.populate_table(app, 'smpls', app.figOpts.smpls.(num));
                tab = Plotting.phn_from_populate(tmp_tab, app, num);

                % Make default first selection
                tab(:, 1) = {false};
                % Deal with empty elements
                indempty = ~cellfun(@isempty,tab(:, 2));
                
                INDmem = ismember(tab(indempty, 2), app.figOpts.phenos.(num));
                tab(~indempty, 1) = {''};
                tab(INDmem, 1) = {true};

                % UI
                app.figOpts.phn_smpl_tab.(num) = uitable(main_fig, ...
                    'Data', tab, ...
                    'ColumnName', {'Select', 'Phenotype', 'Select Sample', 'Sample'}, ...
                    'ColumnEditable', [true, false, true, false] ...
                );
            else
                app.figOpts.phn_smpl_tab.(num).Visible = 'on';
            end

            %% Set Positions
            % Set position based on current main figure.
            width_tab_axis = (main_fig.Position(3) - (2.0 * Constants.plt_menu_x0 + app.figOpts.SldrZ.(num).Position(3))) / 2;
            height_tab_axis = main_fig.Position(4) - (Constants.plt_menu_y0 + app.figOpts.sldrIF.(num).Position(2) + app.figOpts.sldrIF.(num).Position(4));

            app.figOpts.phn_smpl_tab.(num).Position = [...
                Constants.plt_menu_x0, ...
                app.figOpts.sldrIF.(num).Position(2) + app.figOpts.sldrIF.(num).Position(4) + Constants.plt_menu_dy, ...
                width_tab_axis, ...
                height_tab_axis
            ];

            ax = gca;
            % Set axis to current position/extent
            ax.Position(1) = (app.figOpts.phn_smpl_tab.(num).Position(1) + app.figOpts.phn_smpl_tab.(num).Position(3) + 50.0) / main_fig.Position(3);
            ax.Position(2) = (app.figOpts.phn_smpl_tab.(num).Position(2) + 40.0) / main_fig.Position(4);
            ax.Position(3) = 1 - ax.Position(1) - ((Constants.plt_menu_x0 + app.figOpts.SldrZ.(num).Position(3) + 0.0) / main_fig.Position(3));
            ax.Position(4) = 1 - ax.Position(2) - ((Constants.plt_menu_y0 + 0.0) / main_fig.Position(4));
            app.figOpts.phn_smpl_tab.(num).ColumnWidth = {55, 125, 55, 82};
            app.figOpts.phn_smpl_tab.(num).CellEditCallback = @(dd, p) cell_edit(dd, p, app, num, overlay);

            %% Make Exit button
            app.figOpts.close_phn_smpl_tab_btn.(num) = uicontrol(main_fig, 'Style', 'pushbutton', ...
                'String', 'Exit Table', ...
                'Callback', @(~, ~) btn_exit(app, num, main_fig));            
            app.figOpts.close_phn_smpl_tab_btn.(num).Position(3:4) = [91 22];
            app.figOpts.close_phn_smpl_tab_btn.(num).Position(1) = ...
                app.figOpts.phn_smpl_tab.(num).Position(1) + ...
                app.figOpts.phn_smpl_tab.(num).Position(3);
            app.figOpts.close_phn_smpl_tab_btn.(num).Position(2) = ...
                app.figOpts.phn_smpl_tab.(num).Position(2) + ...
                app.figOpts.phn_smpl_tab.(num).Position(4) - ...
                app.figOpts.close_phn_smpl_tab_btn.(num).Position(4);
            Helper.func_SetCLR(app,  app.figOpts.close_phn_smpl_tab_btn.(num), 'UICpopup')


            app.figOpts.order_list_btn.(num) = uicontrol(main_fig, 'Style', 'pushbutton', ...
                'String', 'Show Plot Order', ...
                'Callback', @(~, ~) btn_order_list(app, num));
            app.figOpts.order_list_btn.(num).Position(3:4) = [120 22];
            app.figOpts.order_list_btn.(num).Position(1) = app.figOpts.close_phn_smpl_tab_btn.(num).Position(1);
            app.figOpts.order_list_btn.(num).Position(2) = app.figOpts.close_phn_smpl_tab_btn.(num).Position(2) - ...
                app.figOpts.order_list_btn.(num).Position(4) - Constants.plt_menu_dy;
            Helper.func_SetCLR(app,  app.figOpts.order_list_btn.(num), 'UICpopup')

            % Set widths of columns inside of the table
            app.figOpts.phn_smpl_tab.(num).ColumnWidth{1} = 55;
            app.figOpts.phn_smpl_tab.(num).ColumnWidth{2} = (app.figOpts.phn_smpl_tab.(num).Position(3) - 110) * 0.6;
            app.figOpts.phn_smpl_tab.(num).ColumnWidth{3} = 55;
            app.figOpts.phn_smpl_tab.(num).ColumnWidth{4} = app.figOpts.phn_smpl_tab.(num).Position(3) - sum(cell2mat(app.figOpts.phn_smpl_tab.(num).ColumnWidth(1:3)));
            %% Helper Functions
            function cell_edit(dd, p, app, num, overlay)
                if p.Indices(2) == 3
                    %% Update samples if they are not empty. Otherwise do not change anything and abort
                    ind = cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 3)), 3));
                    if any(ind)
                        % Update selection
                        smpls = dd.Data(~cellfun('isempty', dd.Data(:, 4)), 4);
                        app.figOpts.smpls.(num) = smpls(ind);
                    else
                        % Samples are empty - Go back and don't change anything
                        dd.Data(p.Indices(1), p.Indices(2)) = {p.PreviousData};
                        return;
                    end
                    tmp_dat = Plotting.update_phn_smpl_table(app, num, dd.Data);
                    non_empty = tmp_dat(~cellfun('isempty', tmp_dat(:, 2)), 1);
                    if sum(cell2mat(non_empty)) == 0
                        % New sample not compatible with the chosen phenotype
                        % Fall back to All Cells
                        tmp_dat{1, 1} = true;
                        app.figOpts.phenos.(num) = dd.Data(1, 2);
                    end
                    app.figOpts.phn_smpl_tab.(num).Data = tmp_dat;

                elseif p.Indices(2) == 1
                    ind = cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 1)), 1));
                    if p.Indices(1) == 1
                        dd.Data(ind, 1) = {false};
                        dd.Data(1, 1) = {true};
                        app.figOpts.phenos.(num) = dd.Data(1, 2);
                    else
                        if ~any(ind)
                            % No samples chosen after this one (i.e. nothing to plot). Fall back.
                            dd.Data{p.Indices(1), p.Indices(2)} = p.PreviousData;
                            return;
                        end

                        % If AllCells where chosen, then uncheck them,
                        % since something else got checked now.
                        if dd.Data{1, 1}
                            dd.Data(1, 1) = {false};
                        end

                        if ismember(dd.Data{p.Indices(1), 2}, {'Density/MFI CCN', 'Density/MFI RSN'})
                            % Only MFI can be chosen.
                            % TODO make it so I can
                            % overlay cells and MFI
                            dd.Data(ind, 1) = {false};
                            dd.Data(p.Indices(1), 1) = {true};
                        elseif startsWith(dd.Data{p.Indices(1), 2}, Constants.neigh_tag) && p.NewData == 1
                            % Get current mask
                            ph_chosen = cell2mat(dd.Data(ind, 1));
                            % Update mask to include only ones starting with neighborhood tag
                            ph_chosen_tmp = dd.Data(ind, 2);
                            ph_chosen = ph_chosen_tmp(ph_chosen);
                            ph_chosen = ph_chosen(startsWith(ph_chosen, Constants.neigh_tag));
                            % Scale mask to whole data
                            ph_chosen = ismember(dd.Data(ind, 2), ph_chosen);
                            dd.Data(ind, 1) = Helper.logical2cell(ph_chosen);
                        elseif p.NewData == 1
                            ph_chosen = cell2mat(dd.Data(ind, 1));
                            % Update mask to include only normal phenotypes
                            ph_chosen_tmp = dd.Data(ind, 2);
                            ph_chosen = ph_chosen_tmp(ph_chosen);
                            ph_chosen = ph_chosen(~startsWith(ph_chosen, Constants.neigh_tag));
                            ph_chosen = ph_chosen(~ismember(ph_chosen, {'Density/MFI CCN', 'Density/MFI RSN'}));
                            % Scale mask to whole data
                            ph_chosen = ismember(dd.Data(ind, 2), ph_chosen);
                            dd.Data(ind, 1) = Helper.logical2cell(ph_chosen);
                        end
                        % Get new indices
                        ind = cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 1)), 1));

                        phenos = dd.Data(~cellfun('isempty', dd.Data(:, 2)), 2);
                        app.figOpts.phenos.(num) = phenos(ind);
                    end
                end
                if (p.Indices(2) == 3 || p.Indices(2) == 1) && isfield(app.figOpts.order_stack, num)
                    Plt_Helper.update_plot_order(app, num);
                end
                Plotting.func_plot(app, num , overlay);
            end

            function btn_exit(app, num, main_fig)
                app.figOpts.close_phn_smpl_tab_btn.(num).Visible = 'off';
                app.figOpts.order_list_btn.(num).Visible = 'off';
                app.figOpts.phn_smpl_tab.(num).Visible = 'off';
                app.figOpts.smplPhenoBtn.(num).Visible = 'on';
                Plotting.resize_figure(main_fig, app);
            end

            function btn_order_list(app, num, web)
                if nargin<3
                    web=0;
                end
                alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
                
                ax_child = app.figOpts.plt_ax.(num).Children;
                combined_list = {};
                for child_idx = 1:numel(ax_child)
                    combined_list{child_idx} = ax_child(child_idx).DisplayName;
                end
                % Make sure that order list is synchronized with phenos from table.
                if ~isfield(app.figOpts.order_list, num) || ~ishghandle(app.figOpts.order_list.(num)) % New or previously closed stuff
                    fig = uifigure('Position', alpha*[100 100 300 250]);
                    if web==1
                        fig.Visible='OFF';
                    end

                    app.figOpts.order_list.(num) = uilistbox(fig, ...
                        'Items', combined_list, ...
                        'Position', alpha*[0 50 fig.Position(3) fig.Position(4) - 50] ...
                    );
                    ax = gca;
                    uibutton(fig, 'Position', alpha*[0 25 fig.Position(3) / 2.0 25], 'Text', 'Move Up', 'ButtonPushedFcn', @(btn, p) move_up(app.figOpts.order_list.(num), ax));
                    uibutton(fig, 'Position', alpha*[0 0 fig.Position(3) / 2.0 25], 'Text', 'Move Top', 'ButtonPushedFcn', @(btn, p) move_top(app.figOpts.order_list.(num), ax));
                    uibutton(fig, 'Position', alpha*[fig.Position(3) / 2.0, 25, fig.Position(3) / 2.0, 25], 'Text', 'Move Down', 'ButtonPushedFcn', @(btn, p) move_down(app.figOpts.order_list.(num), ax));
                    uibutton(fig, 'Position', alpha*[fig.Position(3) / 2.0, 0, fig.Position(3) / 2.0, 25], 'Text', 'Move Bottom', 'ButtonPushedFcn', @(btn, p) move_bottom(app.figOpts.order_list.(num), ax));

                    if ~isfield(app.figOpts.order_stack, num)
                        app.figOpts.order_stack.(num) = 1:numel(combined_list);
                        app.figOpts.order_items.(num) = combined_list(app.figOpts.order_stack.(num));
                        app.figOpts.order_list.(num).Items = app.figOpts.order_items.(num);
                    else
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
                        app.figOpts.order_list.(num).Items = app.figOpts.order_items.(num);
                    end
                else
                    if ~Helper.setequal(app.figOpts.order_items.(num), combined_list)
                        new_idxs = ismember(app.figOpts.order_items.(num), combined_list);
                        new_idxs = app.figOpts.order_stack.(num)(new_idxs);
                        if numel(new_idxs) < numel(combined_list)
                            add_idxs = (numel(new_idxs) + 1):numel(combined_list);
                            new_idxs(end + 1:numel(combined_list)) = add_idxs;
                        else % numel(new_idxs) == numel(combined_list)
                            % Element in new_idxs can be higher than numel combined_list
                            % If so, re-index these elements, while keeping relative orderings.
                            over_idxs = new_idxs(new_idxs > numel(combined_list));
                            if ~isempty(over_idxs)
                                new_diff = setdiff(1:numel(new_idxs), new_idxs);
                                srt_old_diff = sort(over_idxs);
                                for so_idx = 1:numel(srt_old_diff)
                                    so_i = srt_old_diff(so_idx);
                                    new_idxs(so_i == new_idxs) = new_diff(so_idx);
                                end
                            end
                        end
                        app.figOpts.order_stack.(num) = new_idxs;
                        app.figOpts.order_items.(num) = combined_list(new_idxs);
                        app.figOpts.order_list.(num).Items = app.figOpts.order_items.(num);
                    end
                end
                
                function move_up(btn, ax)
                    idx = find(strcmp(btn.Items, btn.Value));
                    order = 1:numel(btn.Items);
                    if idx == 1
                        return;
                    end

                    order(idx) = idx - 1;
                    order(idx - 1) = idx;
                    tmp = btn.Items(idx - 1);
                    btn.Items(idx - 1) = {btn.Value};
                    btn.Items(idx) = tmp;
                    ax.Children = ax.Children(order);

% % %                     %order = cat(2, 1, order + 1);
% % %                     ax.Children(end).CData = ax.Children(end).CData(:, order);
% % %                     ax.XTickLabel = ax.XTickLabel(order);
                end

                function move_down(btn, ax)
                    idx = find(strcmp(btn.Items, btn.Value));
                    order = 1:numel(btn.Items);
                    if idx == numel(btn.Items)
                        return;
                    end

                    order(idx) = idx + 1;
                    order(idx + 1) = idx;
                    tmp = btn.Items(idx + 1);
                    btn.Items(idx + 1) = {btn.Value};
                    btn.Items(idx) = tmp;
                    ax.Children = ax.Children(order);
                end

                function move_top(btn, ax)
                    idx = find(strcmp(btn.Items, btn.Value));
                    order = 1:numel(btn.Items);
                    if idx == 1
                        return;
                    end

                    order(idx) = 1;
                    order(1) = idx;
                    tmp = btn.Items(1);
                    btn.Items(1) = {btn.Value};
                    btn.Items(idx) = tmp;
                    ax.Children = ax.Children(order);
% % %                     %order = cat(2, 1, order + 1);
% % %                     ax.Children(end).CData = ax.Children(end).CData(:, order);
% % %                     ax.XTickLabel = ax.XTickLabel(order);
                end

                function move_bottom(btn, ax)
                    idx = find(strcmp(btn.Items, btn.Value));
                    order = 1:numel(btn.Items);
                    if idx == numel(btn.Items)
                        return;
                    end

                    order(idx) = numel(btn.Items);
                    order(end) = idx;
                    tmp = btn.Items(end);
                    btn.Items(end) = {btn.Value};
                    btn.Items(idx) = tmp;
                    ax.Children = ax.Children(order);

% % % %                     %order = cat(2, 1, order + 1);
% % % %                     ax.Children(end).CData = ax.Children(end).CData(:, order);
% % % %                     ax.XTickLabel = ax.XTickLabel(order);
                end
            end
        end
        
        function new_table_data = update_phn_smpl_table(app, num, prev_table_data)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update_phn_smpl_table updates the contents of the
            % Sample/Phenotype table according to what user chooses, as
            % well as modifies app.figOpts to reflect the changes.
            % This allows user to choose what is plotted in an
            % efficient way.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - num - Figure Identifier (For example 'fig1')
            %   - prev_table_data - cell - Previous data from
            %       Sample/Phenotype table.
            %
            % Output:
            %   - new_table_data - Data which should be shown in the new
            %       version of Sample/Phenotype table.
            %
            % Modifies:
            %   - app - Modifies fields in app.figOpts corresponding
            %       choices in the table.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Save phenotypes from last table
            previous_phenos = {};
            for ph_idx=1:size(prev_table_data, 1)
                if ~isempty(prev_table_data{ph_idx, 1}) && prev_table_data{ph_idx, 1} == 1
                    previous_phenos{end + 1} = prev_table_data{ph_idx, 2}; %#ok<*AGROW>
                end
            end
            %% Process Data from Table
            tmp = Plotting.phn_to_populate(prev_table_data);

            tmp = Helper.populate_table(app, 'smpls', app.figOpts.smpls.(num), 'prev_table', tmp);
            tmp = Plotting.phn_from_populate(tmp, app, num);

            if numel(previous_phenos) == 1
                if strcmp(previous_phenos{1}, "All Cells")
                    tmp(1, 1) = {true};
                    tmp(2:sum(~cellfun('isempty', tmp(:, 2))), 1) = {false};
                elseif strcmp(previous_phenos{1}, "Density/MFI CCN")
                    idx_mfi = 0;
                    for ph_idx=1:size(tmp, 1)
                        if ~isempty(tmp{ph_idx, 2}) && strcmp(tmp{ph_idx, 2}, "Density/MFI CCN")
                            idx_mfi = ph_idx;
                        end
                    end
                    if idx_mfi ~= 0
                        tmp(1:sum(~cellfun('isempty', tmp(:, 2))), 1) = {false};
                        tmp(idx_mfi, 1) = {true};
                    end
                elseif strcmp(previous_phenos{1}, "Density/MFI RSN")
                    idx_mfi = 0;
                    for ph_idx=1:size(tmp, 1)
                        if ~isempty(tmp{ph_idx, 2}) && strcmp(tmp{ph_idx, 2}, "Density/MFI RSN")
                            idx_mfi = ph_idx;
                        end
                    end
                    if idx_mfi ~= 0
                        tmp(1:sum(~cellfun('isempty', tmp(:, 2))), 1) = {false};
                        tmp(idx_mfi, 1) = {true};
                    end
                elseif startsWith(previous_phenos{1}, Constants.neigh_tag)
                    nonempty = ~cellfun('isempty', tmp(:, 2));
                    idxs = ismember(tmp(nonempty, 2), previous_phenos);
                    if any(idxs)
                        % Some of neighborhood are still there
                        % Merge idxs and nonempty (they got different sizes)
                        idxs_idx = 1;
                        for ne_idx = 1:numel(nonempty)
                            if nonempty(ne_idx)
                                nonempty(ne_idx) = idxs(idxs_idx);
                                idxs_idx = idxs_idx + 1;
                            end
                        end
                        tmp(1:numel(nonempty), 1) = {false};
                        tmp(nonempty, 1) = {true};
                    end
                end
            elseif all(startsWith(previous_phenos, Constants.neigh_tag))
                nonempty = ~cellfun('isempty', tmp(:, 2));
                idxs = ismember(tmp(nonempty, 2), previous_phenos);
                if any(idxs)
                    % Some of neighborhood are still there
                    % Merge idxs and nonempty (they got different sizes)
                    idxs_idx = 1;
                    for ne_idx = 1:numel(nonempty)
                        if nonempty(ne_idx)
                            nonempty(ne_idx) = idxs(idxs_idx);
                            idxs_idx = idxs_idx + 1;
                        end
                    end
                    tmp(1:numel(nonempty), 1) = {false};
                    tmp(nonempty, 1) = {true};
                end
            end

            %% Update phenotypes in app
            app.figOpts.phenos.(num) = {};
            for ph_idx=1:size(tmp, 1)
                if ~isempty(tmp{ph_idx, 1}) && tmp{ph_idx, 1} == 1
                    app.figOpts.phenos.(num){end + 1} = tmp{ph_idx, 2};
                end
            end

            %% Update Data on plot and in app
            new_table_data = tmp;
        end

        function data = phn_to_populate(data)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % phn_to_populate - takes a data from smpl/phn table, and
            %       remakes it into format applicable for populate table.
            %       It's main purposes is to ensure that 'All Cells' and
            %       MFIs/Neighborhoods do not interfere with populate table
            %       function.
            %       Inverse to phn_to_populate.
            %
            % Input:
            %   data - cell - data currently in phenotype/sample table on
            %       given figure
            %
            % Output:
            %   data - data which can be passed into populate_table
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            valid_idxs = true(size(data, 1), 1);
            for idx_ph = 1:size(data, 1)
                if ~isempty(data{idx_ph, 2})
                    if startsWith(data{idx_ph, 2}, Constants.neigh_tag)
                        valid_idxs(idx_ph) = false;
                    elseif ismember(data(idx_ph, 2), {'All Cells', 'Density/MFI CCN', 'Density/MFI RSN'})
                        valid_idxs(idx_ph) = false;
                    end
                end
            end

            new_ph_choices = data(valid_idxs, 1);
            new_ph = data(valid_idxs, 2);

            count_choices = 0;
            for ph_c_idx=1:numel(new_ph_choices)
                if ~isempty(new_ph_choices{ph_c_idx}) && new_ph_choices{ph_c_idx} == 1
                    count_choices = count_choices + 1;
                end
            end
            if count_choices == 0 && data{1, 1} == 1
                % All Cells chosen
                for ph_c_idx=1:numel(new_ph_choices)
                    if ~isempty(new_ph_choices{ph_c_idx})
                        new_ph_choices{ph_c_idx} = true;
                    end
                end
            elseif count_choices == 0
                new_ph_choices{1} = true;
            end

            size_smpl = sum(~cellfun('isempty', data(:, 4)));
            new_data = cell(max(size_smpl, numel(new_ph_choices)), 5);

            new_data(1:numel(new_ph_choices), 1) = {1};
            new_data(1:numel(new_ph_choices), 2) = new_ph_choices;
            new_data(1:numel(new_ph),         3) = new_ph;
            new_data(1:size_smpl,             4) = data(1:size_smpl, 3);
            new_data(1:size_smpl,             5) = data(1:size_smpl, 4);
            data = new_data;
        end

        function data = phn_from_populate(data, app, num)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % phn_from_populate - Converts data of table returned from
            %       populate table to version which can be used for
            %       phenotype/sample table.
            %       Inverse to phn_to_populate.
            %
            % Input:
            %   - data - cell - data returned from populate table
            %   - app - Instance of CytoMAP
            %   - num - Figure identifier (For example: 'fig1')
            %
            % Output:
            %   - data - cell - 2D cell corresponding to given data, but in
            %       format that is acceptable by the Phenotype/Sample
            %       table.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n_phenos = sum(~cellfun('isempty', data(:, 2)));
            n_samples = sum(~cellfun('isempty', data(:, 4)));
            size_tab = 1;

            %% Figure out which MFI/Neighborhoods work
            rsn_all = true;
            rsn_neigh = {};

            ccn_all = true;
            ccn_neigh = {};

            for s_idx=1:numel(app.figOpts.smpls.(num))
                if ~isfield(app.data.(app.figOpts.smpls.(num){s_idx}), 'MFIRSN')
                    rsn_all = false;
                end
                if ~isfield(app.data.(app.figOpts.smpls.(num){s_idx}), 'MFICCN')
                    ccn_all = false;
                end
            end
            if rsn_all
                size_tab = size_tab + 1;
                rsn_neigh = Helper.get_neighs(app, app.figOpts.smpls.(num), 'MFIRSN');
                size_tab = size_tab + numel(rsn_neigh);
            end
            if ccn_all
                size_tab = size_tab + 1;
                ccn_neigh = Helper.get_neighs(app, app.figOpts.smpls.(num), 'MFICCN');
                size_tab = size_tab + numel(ccn_neigh);
            end

            %% Make Table to be shown
            new_data = cell(max(n_phenos + size_tab, n_samples), 4);

            % Phenos
            % Add All Cells and typical phenotypes
            new_data(1, 1) = {false};
            new_data(1, 2) = {'All Cells'};

            new_data(2:(n_phenos+1), 1) = data(1:n_phenos, 2);
            new_data(2:(n_phenos+1), 2) = data(1:n_phenos, 3);
            current_idx = n_phenos + 2;

            % Add RSN/CCN
            if rsn_all
                new_data(current_idx, 1) = {false};
                new_data(current_idx, 2) = {'Density/MFI RSN'};
                current_idx = current_idx + 1;
                if ~isempty(rsn_neigh)
                    new_data(current_idx:current_idx + numel(rsn_neigh) - 1, 1) = {false};
                    new_data(current_idx:current_idx + numel(rsn_neigh) - 1, 2) = rsn_neigh;
                    current_idx = current_idx + numel(rsn_neigh);
                end
            end
            if ccn_all
                new_data(current_idx, 1) = {false};
                new_data(current_idx, 2) = {'Density/MFI CCN'};
                current_idx = current_idx + 1;
                if ~isempty(ccn_neigh)
                    new_data(current_idx:current_idx + numel(ccn_neigh) - 1, 1) = {false};
                    new_data(current_idx:current_idx + numel(ccn_neigh) - 1, 2) = ccn_neigh;
                end
            end

            % Samples
            new_data(1:size(data, 1), 3) = data(:, 4);
            new_data(1:size(data, 1), 4) = data(:, 5);

            data = new_data;
        end
    end % End of Methods
end % End of Class Definition
