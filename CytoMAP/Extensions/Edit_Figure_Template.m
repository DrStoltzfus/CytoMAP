function Edit_Figure_Template(app)
    ax = gca;
    
    %% Apply the following figure settings to the currently selected figure
    
   % Delete all of the buttons and controls in the figure
    h = findobj(gcf,'type','UIControl');
    delete(h)

    % Define your axes limits
    rangex = 10000;
%     ax.XLim(1) = ax.XLim(1)-0.5*rangex;
    ax.XLim(2) = ax.XLim(1)+rangex;
    
    rangey = 7000;
%     ax.YLim(1) = ax.YLim(1)-0.5*rangey;
    ax.YLim(2) = ax.YLim(1)+rangey;
    
    % Define the colors of your axes labels
    ax.XColor = 'k';
    ax.YColor = 'k';

    % Turn the color bar off/on
    colorbar off

    % Get rid of the title
    ax.Title = [];



end