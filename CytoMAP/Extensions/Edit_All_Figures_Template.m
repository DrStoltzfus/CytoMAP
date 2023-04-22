function Edit_All_Figures_Template(app)
    % Pulls the names of all of the currently open figures
    FigN = get(groot, 'Children');
    if ~isempty(FigN)
        FigN = [FigN.Number];
        for i=1:numel(FigN)
            fig = figure(FigN(i));
            ax=gca;
%%%%%%%%%%%
%%%%%%%%%%%
            %% Apply these figure settings to all open figures
            
%             % Delete all of the buttons and controls in the figure
%             h = findobj(gcf,'type','UIControl');
%             delete(h)

            % Define your axes limits
            % Define your axes limits
            rangex = 15000;
            ax.XLim(2) = ax.XLim(1)+rangex;
            ax.XLim(1) = ax.XLim(1)-0.1*rangex;


            rangey = 15000;
            ax.YLim(2) = ax.YLim(1)+rangey;
            ax.YLim(1) = ax.YLim(1)-0.1*rangey;


% % %             box off
% % %             % Define the colors of your axes labels
% % %             ax.XColor = 'w';
% % %             ax.YColor = 'w';
% % %             
% % %             % Turn the color bar off/on
% % %             colorbar off
% % %             
% % %             % Get rid of the title
% % %             ax.Title = [];

%% Just change limits
%             ax.YLim = [0 500]
%             ax.XLim = [0 8]

        end
    end    
end
