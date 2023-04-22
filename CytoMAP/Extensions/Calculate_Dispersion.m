function Calculate_Dispersion(app)
    %% Calculate the NeighborhoodDistribution

% % %     smpls = fieldnames(app.data);
% % %     % r = 30;
% % %     % Define a surface containing the area of the tissue
% % %     smplarea = 'Tissue';
% % % 
% % %     % results = cell(numel(smpls)+1, 1)'
% % %     results = {'Sample', 'CellType', 'Density', 'NN_Distance', 'Dispersion'};
% % %     for i = 1:numel(smpls)
% % % 
% % %         chAll = Helper.get_others(app, smpls{i});
% % %         chDen = chAll(contains(chAll, 'LocalDensityOf_'));
% % %         chDist = chAll(contains(chAll, 'DistTo'));
% % %         chDisp = [chDist, chDen];
% % % 
% % %         chCell = strrep(chDisp, Constants.other_tag, '');
% % %         chCell = strrep(chCell, 'LocalDensityOf_', '');
% % %         chCell = strrep(chCell, 'DistTo_', '');
% % %         chCell = unique(chCell, 'stable');
% % % 
% % %         dat = app.data.(smpls{i}).MFIRSN;
% % % 
% % %     %     gatenames = Helper.get_gates(app, smpls{i})';
% % %         tags = app.data.(smpls{i}).GateTags.Properties.VariableNames';
% % %         gatenames = app.data.(smpls{i}).GateTags{2, :}';
% % %         gatenames = Helper.valid_var(gatenames);
% % %         for j = 1:numel(chCell)
% % %             % Pull the gate tag of the selected cell
% % %             nm = strsplit( chCell{j}, '_');
% % %             INDg = contains(gatenames, nm{1});
% % %             if numel(nm)>1
% % %                 for k = 2:numel(nm)
% % %                     INDg = INDg.*contains(gatenames, nm{k});
% % %                 end
% % %             end
% % %             nm = gatenames{INDg==1};
% % %             tag = tags{INDg==1};
% % % 
% % %             tagLogic = dat.(tag);
% % % 
% % %             NCells = sum(app.data.(smpls{i}).AllCells.(tag));
% % %             NNghALL = sum(dat.NCells > 0);
% % %             NNgh_j = sum(dat.(tag) > 3);
% % %             % find the expected number of neighborhoods with cells
% % %             if NCells < NNghALL
% % %                 NNexp = NCells;
% % %             else
% % %                 NNexp = NNghALL;
% % %             end
% % %             %%%% want some sort of expected number of neighborhoods with cells
% % %             %%%% compared to actual number of neighborhoods with cells
% % %             dispers  = NNexp/NNgh_j;
% % % 
% % %             results = [results; {smpls{i}, chCell{j}, NNgh_j, NNexp, dispers}];
% % %         end
% % %     end
    %% Calculate the dispersion

    smpls = fieldnames(app.data);
    % r = 30;
    % Define a surface containing the area of the tissue
    smplarea = 'Tissue';   
    
    %generate random points
    for i = 1:numel(smpls)
         % results = cell(numel(smpls)+1, 1)'
        results = {'Sample', 'CellType', 'Dispersion1', 'Dispersion2'};

        chAll = Helper.get_others(app, smpls{i});
        chDen = chAll(contains(chAll, 'LocalDensityOf_'));
        chDist = chAll(contains(chAll, 'DistTo'));
        chDisp = [chDist, chDen];
        dat = app.data.(smpls{i}).AllCells(:, chDisp);

        chCell = strrep(chDisp, Constants.other_tag, '');
        chCell = strrep(chCell, 'LocalDensityOf_', '');
        chCell = strrep(chCell, 'DistTo_', '');
        chCell = unique(chCell, 'stable');
    %     gatenames = Helper.get_gates(app, smpls{i})';
        tags = app.data.(smpls{i}).GateTags.Properties.VariableNames';
        gatenames = app.data.(smpls{i}).GateTags{2, :}';
        gatenames = Helper.valid_var(gatenames);
        
        for j = 1:numel(chCell)
            % Pull the gate tag of the selected cell
            nm = strsplit( chCell{j}, '_');
            INDg = contains(gatenames, nm{1});
            % Restore the name
            if numel(nm)>1
                for k = 2:numel(nm)
                    INDg = INDg.*contains(gatenames, nm{k});
                end
            end
            nm = gatenames{INDg==1};
            tag = tags{INDg==1};
            tagLogic = app.data.(smpls{i}).AllCells.(tag);
            chDistCellj = chDist{contains(chDist, chCell{j})};
            % average distance to neares neighbor for cell type j
            NNDist = mean(table2array(dat(tagLogic==1, chDistCellj)));
            % Average density of cells
            %%%%Neeeds to be total cells by the volume of the smaple
            %app.data.(smpls{i}).GateTags.(chCell{j});
            chCellj = strrep(chCell{j}, 'POS', '+');
            chCellj = strrep(chCellj, 'NEG', '-');
            chCellj = strrep(chCellj, '_', '/');
            chCellj = strrep(chCellj, '/dDC', ' dDC');
            chCellj = strrep(chCellj, '/Macs', ' Macs');
            chCellj = strrep(chCellj, '/Spots', ' Spots');
            chCellj = strrep(chCellj, 'RDP/Ungated', 'RDP_Ungated');
            chCellj = strrep(chCellj, 'All Spots', 'All/Spots');
            
            idx = find(strcmp({chCellj}, table2cell(app.data.(smpls{i}).GateTags(2, :))));
            R = sum(app.data.(smpls{i}).AllCells.(strcat('Gate_',num2str(idx))));
            area = app.data.(smpls{i}).Surfaces.UDS.(smplarea).Surf.area;
            Density =sum(tagLogic)*log(R)/area;
    %         Density = mean(table2array(dat(:, chDenCellj)))/(hght*pi*r^2);
            dispersion = 1/(Density^(1/2)*log(NNDist));
            % Eberhardts metric of clustering
            chDistRDP = chDist{contains(chDist, 'RDP')};
            nPoints = sum(tagLogic);
            NNDistsq = sum(table2array(dat(tagLogic==1, chDistRDP)).^2);
            NNDist2 = (sum(table2array(dat(tagLogic==1, chDistRDP)))/nPoints)^2;
            dispE = NNDistsq/(nPoints*NNDist2);
            results = [results; {smpls{i}, chCell{j}, dispersion, dispE}];
        end
        
        if size(results, 1)>1                  
           %% Plot stuff           
            % First figure             
            fig = figure;
            fig.Color = 'w';
            clf
            title(['Dispersion by method 1 ',smpls{i}]);
            hold on    
            for j = 1:numel(chCell)    
                logic = strcmp(results(:,2),chCell{j});
                plot(zeros(sum(logic),1)+j, [results{logic,3}], 'o', 'DisplayName', chCell{j});
            end
            yl = ylim; % Get current limits.
            ylim([0, yl(2)]); % Replace lower limit only with a y of 0.
            ax = gca;
            ax.XTick = 1:1:numel(chCell);
            ax.XTickLabel = chCell;
            ax.XTickLabelRotation = 30;
            ax.TickLabelInterpreter  = 'none';          
            
            
            % Second figure
            fig = figure;
            fig.Color = 'w';
            clf
            hold on
            title(['Dispersion by method 2 ',smpls{i}]);
            for j = 1:numel(chCell)    
                logic = strcmp(results(:,2),chCell{j});
                plot(zeros(sum(logic),1)+j, [results{logic,4}], 'o', 'DisplayName', chCell{j});
            end
            yl = ylim; % Get current limits.
            ylim([1, yl(2)]); % Replace lower limit only with a y of 0.
            ax = gca;
            ax.XTick = 1:1:numel(chCell);
            ax.XTickLabel = chCell;
            ax.XTickLabelRotation = 30;
            ax.TickLabelInterpreter  = 'none';
        end        
    end
end