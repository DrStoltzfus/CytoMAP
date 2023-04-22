function Overlay_Confocal_Images(app)

% FUNC_IMSLOAD Work in progress
%

    %% add the path to the third party imaris load function        
    CurrentPath = fileparts(mfilename('fullpath'));
    CurrentPath = split(CurrentPath,filesep);
    CurrentPath = [CurrentPath(1:(end-1)); {'3rdPartyFunctions'; 'ImarisReader'}];
    CurrentPath = join(CurrentPath,filesep);            
    addpath(CurrentPath{1});  

    %% Load the image file

    % Select the file
    % returns the filename (fnm) and the path name (pnm)
    [fnmParent, pnm]=uigetfile({'*.ims'}, 'Select the .ims file','Multiselect', 'off');

    % Load the image
    IMSDat = struct;
    IMSDat.imsObj = ImarisReader([pnm fnmParent]);
    IMSDat.ChInfo = IMSDat.imsObj.DataSet.ChannelInfo;
    IMSDat.ChName = {IMSDat.ChInfo.Name};
    
    %% Pull the initial frame
    cIdx = 1;
    tIdx = 0;
%     vol = IMSDat.imsObj.DataSet.GetDataVolume(cIdx-1, tIdx);   
    
    % pull the number of channels
    Nch = numel(IMSDat.ChName);
    % Pull the channel colors
    IMSDat.clr = reshape([IMSDat.ChInfo.Color], [Nch,3]);
    
    %% Load in the image data in a downsampled way
    
% % %     Dreduc_opts = cell(1,2);
% % %     % Options for t-SNE
% % %     Dreduc_opts(:,1) = {'Downsample ratio'};
% % %     % Defaults
% % %     Dreduc_opts(:,2) = {'1'};
% % %     
% % %     dlg_title = 'Input Downsample ratio';
% % %     num_lines = 1;
% % %     vAnswer = inputdlg(Dreduc_opts(:,1),dlg_title,num_lines,Dreduc_opts(:,2));
% % %     if isempty(vAnswer)
        ratio = 1;
% % %     else
% % %         ratio = str2double(vAnswer);
% % %     end
    
    %define the end size of the image
    numrows = ceil(1024 / ratio);
    numcols = ceil(1024 / ratio);
    
    % Pull the first volume
    vol = IMSDat.imsObj.DataSet.GetDataVolume(0, tIdx);
    IMSDat.ImageData = zeros(Nch, numrows, numcols, size(vol, 3));
    
    ['Image resized from: ' num2str(size(vol)) ' to: ' num2str([numrows, numcols, size(vol,3)])]

    for ch_i = 1:Nch 
        % pull the ith volume
        if ch_i ~=1
            vol = IMSDat.imsObj.DataSet.GetDataVolume(ch_i-1, tIdx);
        end
        
        for z_i = 1:size(vol, 3)
            IMSDat.ImageData(ch_i, :, :, z_i) = imresize(vol(:,:,z_i),[numrows numcols]); 
        end
    end
    
    % resize the image based on the pixel size
    dX = (IMSDat.imsObj.DataSet.ExtendMaxX-IMSDat.imsObj.DataSet.ExtendMinX)/(numcols);
    dY = (IMSDat.imsObj.DataSet.ExtendMaxY-IMSDat.imsObj.DataSet.ExtendMinY)/(numrows);
%     X = flip((1:numcols)*dX) + IMSDat.imsObj.DataSet.ExtendMinX;
%     Y = flip((1:numrows)*dY) + IMSDat.imsObj.DataSet.ExtendMinY;
    X = flip((1:numcols)*dX);
    Y = flip((1:numrows)*dY);

    IMSDat.XData = X;
    IMSDat.YData = Y;
    
    %% Build the figure
    
    Plotting.func_newfig(app)
    fig = gcf;
    num = [ 'fig' num2str(fig.Number)]; %(Use this to assign handles)
    cla

    % Pull the first channel
    vol = reshape(IMSDat.ImageData(1, :, :, :), ...
        [size(IMSDat.ImageData(1, :, :, :), 2), ...
         size(IMSDat.ImageData(1, :, :, :), 3), ...
         size(IMSDat.ImageData(1, :, :, :), 4)]);
    % Build an image summed along Z    
    frame = cat(3, IMSDat.clr(1, 1).*sum(vol(:, :, :), 3), ...
                   IMSDat.clr(1, 2).*sum(vol(:, :, :), 3), ...
                   IMSDat.clr(1, 3).*sum(vol(:, :, :), 3));
    % normalize the frame
    frame = 1.*(frame./max(max(max(frame))));
    % make the orientation of the image normal
    frame = permute(frame, [2,1,3]);
    frame = flip(frame, 2);
    frame = flip(frame, 1);
    
    % generate the image
    image(IMSDat.XData, IMSDat.YData, frame);
% % %     image(X, Y, frame);

    %% Tweak the figure to be better suited to confocal images
    
    % invert the figure color
    invert = struct;
    invert.State = 'on';
    Plt_Helper.pltIC(invert, app, num);
    
    %Discard unnescesary options
    app.figOpts.xaxTXT.(num).Visible = 'off';
    app.figOpts.xaxIF.(num).Visible = 'off';    
    app.figOpts.yaxTXT.(num).Visible = 'off';
    app.figOpts.yaxIF.(num).Visible = 'off';
    app.figOpts.zaxTXT.(num).Visible = 'off';
    app.figOpts.zaxIF.(num).Visible = 'off';
    app.figOpts.sldrIF.(num).Visible = 'off';
    app.figOpts.sldrIFTXT.(num).Visible = 'off';
    
    %% these might be added back in later if their functionality is fixed
    app.figOpts.opts.(num).Visible = 'off';
    
    app.figOpts.smplPhenoBtn.(num).String = 'Phenotypes Table';
    %app.figOpts.smplPhenoBtn.(num).Position(3:4) = app.figOpts.smplPhenoBtn.(num).Extent(3:4);
    app.figOpts.smplPhenoBtn.(num).Position(3:4) = [100,40];
    
    %% Add overlay option to plot functions
    app.figOpts.smplPhenoBtn.(num).Callback =  @(~, ~) Plotting.show_phn_smpl_table(app, num, fig, 1);
    app.figOpts.caxIF.(num).Callback = @(~,~) Plotting.func_plot(app, num, 1);
% %     app.figOpts.SldrZ.(num).PostSet = @(~, ~) Plotting.func_plot(app, num, 1);
    
    % Create a channels button
	%% make a table that lists the channels
    app.figOpts.channelBtn.(num) = uicontrol('Style', 'pushbutton', ...
                'String', 'Channel Table', ...
                'Position', [Constants.plt_menu_x0, Constants.plt_menu_x0+40, 100, 40], ...
                'Callback', @(~, ~) show_phn_channel_table(app, num, fig,IMSDat));
    %app.figOpts.channelBtn.(num).Position(3:4) = app.figOpts.channelBtn.(num).Extent(3:4);
    Helper.func_SetCLR(app,  app.figOpts.channelBtn.(num), 'UICpushbutton')
    
    %% Create options button
    app.figOpts.Confopts.(num) = uicontrol('Style', 'pushbutton', ...
                'String', 'Options');
    app.figOpts.Confopts.(num).Position =  [Constants.plt_menu_x0 Constants.plt_menu_x0+0 100 40];
    app.figOpts.Confopts.(num).Callback = @(~, ~) overlay_options(app, num, IMSDat);

    %     Helper.func_SetCLR(app,  app.figOpts.Confopts.(num), 'UICpushbutton')
    
    %% Add a slidebar
% % %     SliderH = uicontrol('Parent',fig,'style','slider','position',[10 10 20 500],...
% % %         'min', 1, 'max', zDim, 'value', zDim/2);
% % % 
% % %     addlistener(SliderH, 'Value', 'PostSet', @(source, eventdata) callbackfn(SliderH, vol, color));
% % % 
% % %     function callbackfn(SliderH, vol, color)
% % %         num = round(SliderH.Value);
% % %         frame = cat(3, color(1).*sum(vol(:, :, num), 3), color(2).*sum(vol(:, :, num), 3), color(3).*sum(vol(:, :, num), 3));
% % %         im = image(1.5.*frame./max(max(max(frame))));
% % %         axis equal
% % %         axis tight
% % %     end
% % %     function pick_color
% % %         color = uisetcolor(color);
% % %     end

    %% Add a opacity slider
    
    
    %%
end % end of main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define internal functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function overlay_options(app, num, IMSDat)
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

    if ~isfield(app.figOpts, 'ims_offset')
        app.figOpts.ims_offset = struct;
    end
    
    alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

    ax2 = gca;
% % %     app.figOpts.limits.(num) = cell(4,5);
    tdat = cell(6,5);
    if ~isempty(ax2.Colorbar)
        tdat(1:4, 1) = {ax2.XLim(1), ax2.YLim(1), ax2.ZLim(1), ax2.Colorbar.Limits(1)};
        tdat(1:4, 2) = {ax2.XLim(2), ax2.YLim(2), ax2.ZLim(2), ax2.Colorbar.Limits(2)};
    else
        tdat(1:3, 1) = {ax2.XLim(1), ax2.YLim(1),ax2.ZLim(1)};
        tdat(1:3, 2)= {ax2.XLim(2), ax2.YLim(2),ax2.ZLim(2)};
    end

    tdat(1:6, 3) = {'X-Axis', 'Y-Axis', 'Z-Axis', 'Color-Axis', 'X-shift', 'Y-shift'};
    tdat(5:6, 2) = {0};
    tdat(:, 4) = {false};
    tdat(:, 5) = {false};

    UIfig = uifigure('Name', 'Plot Options');


    UIfig.Position(3:4) = [350 300];
    % Create the table of options
    t = uitable(UIfig);
    t.Data = tdat;
    t.Position = alpha*[0 40 350 260];
    t.ColumnName = {'Min','Max','Axis', 'Auto', 'Log'};
    t.ColumnEditable = [true true false, true, true];
    t.ColumnWidth = {alpha*75 alpha*75 alpha*75, alpha*50, alpha*50};

    % Create a push button
    btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) overlay_options_backend( ...
        app, UIfig, t, IMSDat ...
    ));
%         cell2mat(t.Data(1:4, 4)), ... Whether ranges are set to automode, or manual
%         cell2mat(t.Data(1:4, 1:2)), ... Manual ranges
%         cell2mat(t.Data(1:4, 5)) ... Whether scale if log or linear
        
    btn.Position = alpha*[300/2-25, 5, 50, 30];
    btn.Text = 'Apply';
    Helper.func_SetCLR(app, btn, 'button')
end % end options button

function show_phn_channel_table(app, num, main_fig, IMSDat)
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
    
    if ~isfield(app.figOpts, 'smpls') || ~isfield(app.figOpts.smpls, num) || isempty(app.figOpts.smpls.(num))
        app.figOpts.smpls.(num) = {app.DataN.Value};
    end

    % Get current phenotypes in selected samples if they do not exist
    if ~isfield(app.figOpts, 'phenos') || ~isfield(app.figOpts.phenos, num) || isempty(app.figOpts.phenos.(num))
        smpl = app.figOpts.smpls.(num);
        smpl = smpl{1};
        app.figOpts.phenos.(num) = app.data.(smpl).GateTags(2, 1);
    end
    
    if ~isfield(app.figOpts, 'ims_channel_tab')
        app.figOpts.ims_channel_tab = struct;
    end
    
    % Load or make new ui table
    if ~isfield(app.figOpts.ims_channel_tab, num)
        %% Generate initial table and it's properties
        
        Nch = numel(IMSDat.ChName);
        clr = reshape([IMSDat.ChInfo.Color], [Nch,3]);
        
        tab = cell(max([numel(IMSDat.ChName), 10]), 8);
        tab(1, 1) ={true};
        tab(2:numel(IMSDat.ChName), 1) = {false};
        tab(1:numel(IMSDat.ChName), 2) = IMSDat.ChName;
        tab(1:numel(IMSDat.ChName), 3) = {0};
        tab(1:numel(IMSDat.ChName), 4) = {255};
        % Edit the color axis
        tab(1:numel(IMSDat.ChName), 5) ={false};
        %Red
        tab(1:numel(IMSDat.ChName), 6) = num2cell(clr(:,1));
        %Blue
        tab(1:numel(IMSDat.ChName), 7) = num2cell(clr(:,2));
        %Green
        tab(1:numel(IMSDat.ChName), 8) =num2cell(clr(:,3));

        % UI
        app.figOpts.ims_channel_tab.(num) = uitable(main_fig, ...
            'Data', tab, ...
            'ColumnName', {'','Channel','min', 'max', 'Color', 'R', 'B', 'G'}, ...
            'ColumnEditable', [true, false, true, true, true, true, true, true], ...
            'ColumnWidth', {40 55 40, 40, 45, 15, 15, 15} ...
        );
    else
        app.figOpts.ims_channel_tab.(num).Visible = 'on';
    end

    %% Set Positions
    % Set position based on current main figure.
    width_tab_axis = (main_fig.Position(3) - (2.0 * Constants.plt_menu_x0 + app.figOpts.SldrZ.(num).Position(3))) / 2;
    height_tab_axis = main_fig.Position(4) - (Constants.plt_menu_y0 + app.figOpts.sldrIF.(num).Position(2) + app.figOpts.sldrIF.(num).Position(4));

    app.figOpts.ims_channel_tab.(num).Position = [...
        Constants.plt_menu_x0, ...
        app.figOpts.sldrIF.(num).Position(2) + app.figOpts.sldrIF.(num).Position(4) + Constants.plt_menu_dy, ...
        width_tab_axis, ...
        height_tab_axis
    ];

    ax = gca;
    % Set axis to current position/extent
    ax.Position(1) = (app.figOpts.ims_channel_tab.(num).Position(1) + app.figOpts.ims_channel_tab.(num).Position(3) + 50.0) / main_fig.Position(3);
    ax.Position(2) = (app.figOpts.ims_channel_tab.(num).Position(2) + 40.0) / main_fig.Position(4);
    ax.Position(3) = 1 - ax.Position(1) - ((Constants.plt_menu_x0 + app.figOpts.SldrZ.(num).Position(3) + 0.0) / main_fig.Position(3));
    ax.Position(4) = 1 - ax.Position(2) - ((Constants.plt_menu_y0 + 0.0) / main_fig.Position(4));
    app.figOpts.ims_channel_tab.(num).ColumnWidth = {40 55 40, 40, 45, 15, 15, 15};
    app.figOpts.ims_channel_tab.(num).CellEditCallback = @(dd, p) cell_edit(dd, p, app, num, IMSDat,main_fig);

    %% Make Exit button
    app.figOpts.close_phn_smpl_tab_btn.(num) = uicontrol(main_fig, 'Style', 'pushbutton', ...
        'String', 'Exit Table', ...
        'Callback', @(~, ~) btn_exit(app, num, main_fig));            
    app.figOpts.close_phn_smpl_tab_btn.(num).Position(3:4) = [91 22];
    app.figOpts.close_phn_smpl_tab_btn.(num).Position(1) = ...
        app.figOpts.ims_channel_tab.(num).Position(1) + ...
        app.figOpts.ims_channel_tab.(num).Position(3);
    app.figOpts.close_phn_smpl_tab_btn.(num).Position(2) = ...
        app.figOpts.ims_channel_tab.(num).Position(2) + ...
        app.figOpts.ims_channel_tab.(num).Position(4) - ...
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
    app.figOpts.ims_channel_tab.(num).ColumnWidth{2} = (app.figOpts.ims_channel_tab.(num).Position(3) - 75) * 0.6;
% %     app.figOpts.ims_channel_tab.(num).ColumnWidth{3} = 50;
% %     app.figOpts.ims_channel_tab.(num).ColumnWidth{4} = app.figOpts.ims_channel_tab.(num).Position(3) - sum(cell2mat(app.figOpts.ims_channel_tab.(num).ColumnWidth(1:3)));

end % end call channel table function

%% Helper Functions
function overlay_options_backend(app, UIfig, t, IMSDat)
    %%    
    fig_i = gcf;
    ax_i = gca;
    
    xoffset = t.Data{5, 2};
    yoffset = t.Data{6, 2};
    
%     xoffset = 0;
%     yoffset = 0;
%    min(fig.CurrentAxes.Children(end).XData)
%    max(fig.CurrentAxes.Children(end).XData)
%    
%    min(fig.CurrentAxes.Children(end).YData)
%    max(fig.CurrentAxes.Children(end).YData)

    IMSDat.XData = fig_i.CurrentAxes.Children(end).XData + xoffset;
    IMSDat.YData = fig_i.CurrentAxes.Children(end).YData + yoffset;

    fig_i.CurrentAxes.Children(end).XData = IMSDat.XData;
    fig_i.CurrentAxes.Children(end).YData = IMSDat.YData;
    % Shift the axes limits
    if ~isempty(ax_i.Colorbar)
        ax_i.XLim(1) = t.Data{1, 1};
        ax_i.YLim(1) = t.Data{2, 1};
        ax_i.ZLim(1) = t.Data{3, 1};
        ax_i.Colorbar.Limits(1) = t.Data{4, 1};
        
        ax_i.XLim(2) = t.Data{1, 2};
        ax_i.YLim(2) = t.Data{2, 2};
        ax_i.ZLim(2) = t.Data{3, 2};
        ax_i.Colorbar.Limits(2) = t.Data{4, 2};
    else
        ax_i.XLim(1) = t.Data{1, 1};
        ax_i.YLim(1) = t.Data{2, 1};
        ax_i.ZLim(1) = t.Data{3, 1};
        ax_i.XLim(2) = t.Data{1, 2};
        ax_i.YLim(2) = t.Data{2, 2};
        ax_i.ZLim(2) = t.Data{3, 2};
    end
    


end

function cell_edit(dd, p, app, num, IMSDat,main_fig)
    % try and make this faster
    t0 = tic();

    % Currently does not work for time data
    tIdx = 0;
    % I need to pull this from the slider
    zposi = 1;

    if p.Indices(2) == 1 
        % If a new channel is added

        % Pull which cannels were selected
        cIdx = find([dd.Data{:, 1}]);
        chclr = zeros(numel([dd.Data{:, 1}]), 3);
        chclr(:,1) = [dd.Data{:, 6}]';
        chclr(:,2) = [dd.Data{:, 8}]';
        chclr(:,3) = [dd.Data{:, 7}]';

        N_i = 0;
        for ch_i = cIdx
            % pull the ith volume
% % %             vol = IMSDat.imsObj.DataSet.GetDataVolume(ch_i-1, tIdx);
            vol = reshape(IMSDat.ImageData(ch_i, :, :, :), ...
                 [size(IMSDat.ImageData(ch_i, :, :, :), 2), ...
                  size(IMSDat.ImageData(ch_i, :, :, :), 3), ...
                  size(IMSDat.ImageData(ch_i, :, :, :), 4)]);
              
            frame_i = cat(3, chclr(ch_i, 1).*sum(vol(:, :, :), 3), ... %red
               chclr(ch_i, 2).*sum(vol(:, :, :), 3), ... %blue
               chclr(ch_i, 3).*sum(vol(:, :, :), 3));    %green

            % normalize the frame
            frame_i = frame_i./max(max(max(frame_i)));
            % sum the channels
            if N_i==0
                frame = frame_i;
            else
                frame = (frame + frame_i);
            end
           N_i = N_i+1;
        end

        frame = permute(frame, [2,1,3]);
        frame = flip(frame, 2);
        frame = flip(frame, 1);
        
% % %         %Resize the image
% % %         frame = imresize(frame, 0.5);

% % %         % resize the image based on the pixel size
% % %         dX = (IMSDat.imsObj.DataSet.ExtendMaxX-IMSDat.imsObj.DataSet.ExtendMinX)/(size(frame,2));
% % %         dY = (IMSDat.imsObj.DataSet.ExtendMaxY-IMSDat.imsObj.DataSet.ExtendMinY)/(size(frame,1));
% % %         X = flip((1:size(frame, 2))*dX);
% % %         Y = flip((1:size(frame, 1))*dY);

        main_fig.CurrentAxes.Children(end).CData = frame;
% % %         main_fig.CurrentAxes.Children(end).XData = X;
% % %         main_fig.CurrentAxes.Children(end).YData = Y;
        
% % %         % replace the image frame            
% % %         main_fig.CurrentAxes.Children(end).CData = frame;

    elseif p.Indices(2) == 3 ...
            || p.Indices(2) == 4

    elseif p.Indices(2) == 5
        % user selected the change color button


    elseif p.Indices(2) == 6 ...
            || p.Indices(2) == 7 ...
            || p.Indices(2) == 8 
        % user changed te color of a channel


    end
    
    % try and make this faster
    runtime = toc(t0);
end

function btn_exit(app, num, main_fig)
    app.figOpts.close_phn_smpl_tab_btn.(num).Visible = 'off';
    app.figOpts.order_list_btn.(num).Visible = 'off';
    app.figOpts.ims_channel_tab.(num).Visible = 'off';
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
        uibutton(fig, 'Position', alpha*[0 25 fig.Position(3) / 2.0 25], 'Text', 'Move Up', 'ButtonPushedFcn', @(btn, p) move_up(app, num));
        uibutton(fig, 'Position', alpha*[0 0 fig.Position(3) / 2.0 25], 'Text', 'Move Top', 'ButtonPushedFcn', @(btn, p) move_top(app, num));
        uibutton(fig, 'Position', alpha*[fig.Position(3) / 2.0, 25, fig.Position(3) / 2.0, 25], 'Text', 'Move Down', 'ButtonPushedFcn', @(btn, p) move_down(app, num));
        uibutton(fig, 'Position', alpha*[fig.Position(3) / 2.0, 0, fig.Position(3) / 2.0, 25], 'Text', 'Move Bottom', 'ButtonPushedFcn', @(btn, p) move_bottom(app, num));

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

    function move_up(app, num)
        idx = find(strcmp(app.figOpts.order_items.(num), ...
            app.figOpts.order_list.(num).Value));
        if idx == 1
            return;
        end
        tmp = app.figOpts.order_items.(num)(idx - 1);
        app.figOpts.order_items.(num)(idx - 1) = ...
            app.figOpts.order_items.(num)(idx);
        app.figOpts.order_items.(num)(idx) = tmp;

        % Backend stack, actually used for ordering the plots/scatter
        tmp = app.figOpts.order_stack.(num)(idx - 1);
        app.figOpts.order_stack.(num)(idx - 1) = ...
            app.figOpts.order_stack.(num)(idx);
        app.figOpts.order_stack.(num)(idx) = tmp;
        app.figOpts.order_list.(num).Items = app.figOpts.order_items.(num);
        Plotting.func_plot(app, num);
    end

    function move_down(app, num)
        idx = find(strcmp(app.figOpts.order_items.(num), ...
            app.figOpts.order_list.(num).Value));
        if idx == numel(app.figOpts.order_items.(num))
            return;
        end
        tmp = app.figOpts.order_items.(num)(idx + 1);
        app.figOpts.order_items.(num)(idx + 1) = ...
            app.figOpts.order_items.(num)(idx);
        app.figOpts.order_items.(num)(idx) = tmp;

        % Backend stack, actually used for ordering the plots/scatter
        tmp = app.figOpts.order_stack.(num)(idx + 1);
        app.figOpts.order_stack.(num)(idx + 1) = ...
            app.figOpts.order_stack.(num)(idx);
        app.figOpts.order_stack.(num)(idx) = tmp;
        app.figOpts.order_list.(num).Items = app.figOpts.order_items.(num);
        Plotting.func_plot(app, num);
    end

    function move_top(app, num)
        idx = find(strcmp(app.figOpts.order_items.(num), ...
            app.figOpts.order_list.(num).Value));
        if idx == 1
            return;
        end
        tmp = app.figOpts.order_items.(num)(idx);
        app.figOpts.order_items.(num)(2:idx) = ...
            app.figOpts.order_items.(num)(1:idx-1);
        app.figOpts.order_items.(num)(1) = tmp;

        % Backend stack, actually used for ordering the plots/scatter
        tmp = app.figOpts.order_stack.(num)(idx);
        app.figOpts.order_stack.(num)(2:idx) = ...
            app.figOpts.order_stack.(num)(1:idx-1);
        app.figOpts.order_stack.(num)(1) = tmp;
        app.figOpts.order_list.(num).Items = app.figOpts.order_items.(num);
        Plotting.func_plot(app, num);
    end

    function move_bottom(app, num)
        idx = find(strcmp(app.figOpts.order_list.(num).Items, ...
            app.figOpts.order_list.(num).Value));
        if idx == numel(app.figOpts.order_list.(num).Items)
            return;
        end
        tmp = app.figOpts.order_list.(num).Items(idx);
        app.figOpts.order_list.(num).Items(idx:end - 1) = ...
            app.figOpts.order_list.(num).Items(idx + 1:end);
        app.figOpts.order_list.(num).Items(end) = tmp;

        % Backend stack, actually used for ordering the plots/scatter
        tmp = app.figOpts.order_stack.(num)(idx);
        app.figOpts.order_stack.(num)(idx:end - 1) = ...
            app.figOpts.order_stack.(num)(idx + 1:end);
        app.figOpts.order_stack.(num)(end) = tmp;
        app.figOpts.order_list.(num).Items = app.figOpts.order_items.(num);
        Plotting.func_plot(app, num);
    end
end

