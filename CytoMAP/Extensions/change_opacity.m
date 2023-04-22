function change_opacity(app)

    % Pull the names of the plotted elements
    ax = gca;
    names = cell(numel(ax.Children), 1);
    for ch_i = 1: numel(ax.Children)
        if strcmp(ax.Children(ch_i).Type, 'image')
            names{ch_i} = 'Image';
        else
            names{ch_i} = ax.Children(ch_i).DisplayName;
        end
    end
    
    [indx,tf] = listdlg('PromptString',{'Select elements of plot to adjust.',''},...
    'SelectionMode','multiple','ListString',names);
    if ~tf
        return
    end
    
    % Create the slider control
    slider = uicontrol('style','slider',...
                    'Min',0,...
                    'Value',1,...
                    'Max',1,...
                    'Position',[10 65 15 330],...
                    'Visible','on');
    slider.Callback = @(~,~) hscroll_Callback(slider, names(indx), ax);
% % %     addlistener(app.figOpts.SldrZ.(num), 'Value', 'PostSet', @(~, ~) Plotting.func_plot(app, num));

    
    function hscroll_Callback(slider, names, ax)
        
        for ch_i = 1: numel(ax.Children)
            if strcmp(ax.Children(ch_i).Type, 'image')
                if contains('Image', names)
                    fieldnames(ax.Children(ch_i));
                    alpha(ax.Children(ch_i), slider.Value^3);
                end
            else
                if contains(ax.Children(ch_i).DisplayName, names)
% % %                     fieldnames(ax.Children(ch_i))
%                     ax.Children(ch_i).MarkerEdgeAlpha = slider.Value;
%                     ax.Children(ch_i).MarkerFaceAlpha = slider.Value;
                    alpha(ax.Children(ch_i), slider.Value^3);

                end
            end
        end
% %         alpha(src.Value^3) 
    end
end