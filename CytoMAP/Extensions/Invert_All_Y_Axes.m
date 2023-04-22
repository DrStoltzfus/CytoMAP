function Invert_All_Y_Axes(app)
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
            ax.YDir = 'reverse';

        end
    end    
end
