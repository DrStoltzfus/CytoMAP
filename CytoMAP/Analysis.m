 classdef Analysis

% Analysis defines a suite of functions used perform
% calculations on your data
%
% Last edited by CS 2019/01/24

    methods(Static)

        function Import_Definitions_Func
% % %             import Helper.*;
% % %             import Plotting.*;
% % %             import Constants.*;
        end

        function dist_func(app, web)
            % DIST_FUNC Front-end to distance function
            % Makes table for user to choose on which phenotype/sample
            % combinations to calculate distance to.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %
            % Modifies:
            %   - app - Adds a distance to things chosen by user in samples
            %       chosen by user.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            if ~Helper.any_sample(app)
                return;
            end

            %%  Pull the phenotype names from the data table
            smpls = fieldnames(app.data);
            special_names = {'Spots', 'Regions', 'Polygons', 'User Defined Surfaces', 'Select All'};
            tmp = Helper.populate_table(app, 'smpls', smpls{1});
            n_pheno = sum(~cellfun('isempty', tmp(:, 3)));
            if n_pheno < size(tmp, 1)
                % More samples than phenotypes
                tbldat = cell(max(size(tmp, 1) + 1, n_pheno + 5), 4);

                tbldat(1:n_pheno, 1)= {false};
                tbldat(end-4:end, 2) = {false};

                tbldat(1:n_pheno, 2) = tmp(1:n_pheno, 3);
                tbldat(end - 4:end, 1) = {false};
                tbldat(end - 4:end, 2) = special_names';

                tbldat(1:numel(smpls), 3) = {false};
                tbldat(1, 3) = {true};
                tbldat(1:size(tmp, 1), 4) = tmp(:, 5);
                tbldat(end, 3) = {false};
                tbldat(end, 4) = {'Select All'};
            else
                % Number of phenotypes equal to size of table
                tbldat = cell(size(tmp, 1) + 5, 4);

                tbldat(:, 1)= {false};

                tbldat(1:size(tmp, 1), 2) = tmp(:, 3);
                tbldat(end - 4:end, 2) = special_names';

                tbldat(1:numel(smpls), 3) = {false};
                tbldat(1, 3) = {true};
                tbldat(end, 3) = {false};
                tbldat(1:size(tmp, 1), 4) = tmp(:, 5);
                tbldat(end, 4) = {'Select All'};
            end


            % Find the names of the polygons
            if max(contains(fieldnames(app), 'polygons'))
                PolyNames = fieldnames(app.polygons);
                IND = strcmp(PolyNames, 'TMP')+strcmp(PolyNames,'PNTSx')+strcmp(PolyNames, 'PNTSy')+strcmp(PolyNames, 'n');
                PolyNames = PolyNames(~IND);
                PolyNames = [{'Polygons', 'none'}, PolyNames(:)'];
                if isempty(PolyNames)
                    PolyNames = {'Polygons', 'none'};
                end
            else
                PolyNames = {'Polygons', 'none'};
            end

            % Find the names of the user defined surfaces
            if max(contains(fieldnames(app.data.(smpls{1}).Surfaces), 'UDS'))
                UDSNames = fieldnames(app.data.(smpls{1}).Surfaces.UDS);
                UDSNames = [{'Surfaces', 'none'}, UDSNames(:)'];
                if isempty(UDSNames)
                    UDSNames = {'Surfaces', 'none'};
                end
            else
                UDSNames = {'Surfaces', 'none'};
            end

            UIfig = uifigure('Name', 'Calculate Distance', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end

            % Select User defined surfaces drop down
            SurfN = uidropdown(UIfig);
            SurfN.Items = UDSNames(2:end);
            SurfN.Value = UDSNames{2};
            SurfN.Position = alpha*[10+100 1 100 30];
            Helper.func_SetCLR(app, SurfN, 'UIField')
            % Select Surface label
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Surface:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+100 32 200 15];
            
            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 800 800];
            t = uitable(UIfig);
            t.ColumnFormat = ({[], PolyNames, [], []});
            t.Data = tbldat;
            t.Position = alpha*[0 50 800 750];
            t.ColumnName = {'Calculate distance to:', 'Phenotype', 'For', 'Sample'};
            t.ColumnEditable = [true, false, true, false];
            t.ColumnWidth = {alpha*150 alpha*350 alpha*50 alpha*250};
            t.CellEditCallback = @(dd, p) switched_smpl(app, special_names, SurfN, dd, p);
            Helper.func_SetCLR(app, t, 'table')

            % Select Polygon drop down
            ponlyN = uidropdown(UIfig);
            ponlyN.Items = PolyNames(2:end);
            ponlyN.Value = PolyNames{2};
            ponlyN.Position = alpha*[10 1 100 30];
            Helper.func_SetCLR(app, ponlyN, 'UIField')
            % Select Polygon label
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Polygon:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10 32 100 15];

            % Select Model drop down
            ModelN = uidropdown(UIfig);
            ModelN.Items = [{'none'}; Helper.full_var(fieldnames(app.net))];
            ModelN.Value = ModelN.Items{1};
            ModelN.Position = alpha*[10+100+100 1 100 30];
            Helper.func_SetCLR(app, ModelN, 'UIField')
            % Select Model label
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Model:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+100+100 32 200 15];

            % Make a radius of neighborhood numeric option
            lbl = uilabel(UIfig); lbl.Text = 'Radius of local Density Calculation:';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[235+10+100 32 300 15];
            Den = uieditfield(UIfig,'numeric');
            Den.Value = app.rwindowRSN;
            Den.Position = alpha*[235+10+200 1 80 30];
            % Create a Density button
            DensButton = uibutton(UIfig,'state');
            DensButton.Position = alpha*[235+10+100, 1, 100, 30];
            DensButton.Text = 'Include Density';
            DensButton.Value = false;
            Helper.func_SetCLR(app, DensButton, 'button');

            % Create a Start button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap( ...
                app, t, UIfig, ponlyN, DensButton, Den, SurfN, ModelN, special_names ...
            ));
            btn.Position = alpha*[800-60, 10, 50, 30];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')

            function switched_smpl(app, special_names, SurfN, dd, p)
                if p.Indices(2) == 3
                    % Make sure that at least one thing is selected
                    if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 3)), 3)))
                        dd.Data{p.Indices(1), 3} = true;
                        return;
                    end
                    
                    % Process 'Select All' button
                    fill_select_all = false;
                    if p.Indices(1) == find(strcmp(dd.Data(:, 4), 'Select All'))
                        ind = ~cellfun('isempty', dd.Data(:, 3));
                        dd.Data(ind, 3) = {logical(p.NewData)};
                        % Make sure that at least one thing is selected
                        if ~logical(p.NewData)
                            dd.Data{strcmp(dd.Data(:, 4), app.DataN.Value), 3} = true;
                        end
                        fill_select_all = logical(p.NewData);
                    end
                    ind = ~cellfun('isempty', dd.Data(:, 3));
                    ind(end) = false;
                    ind(ind == 1) = logical(cell2mat(dd.Data(ind, 3)));

                    % Extract special ones (broadcast logic)
                    logic = ~cellfun('isempty', dd.Data(:, 2));
                    logic(logic) = ismember(dd.Data(logic, 2), special_names);
                    keep_special = dd.Data(logic, 1:2);

                    tmp_s = cell(size(dd.Data, 1) - 5, 5);

                    tmp_s(:, 1) = {1};
                    tmp_s(:, 2) = dd.Data(1:end - 5, 1);
                    tmp_s(:, 3) = dd.Data(1:end - 5, 2);

                    tmp_s(:, 4) = dd.Data(1:end - 5, 3);
                    tmp_s(:, 5) = dd.Data(1:end - 5, 4);
                    tmp_s = Helper.populate_table(app, ...
                        'smpls', dd.Data(ind, 4), ...
                        'prev_table', tmp_s ...
                    );

                    non_empty_pheno = sum(~cellfun('isempty', tmp_s(:, 2)));
                    if  non_empty_pheno < numel(app.DataN.Items)
                        % Less phenotypes than samples.
                        dd.Data = cell(max(non_empty_pheno + size(keep_special, 1), numel(app.DataN.Items) + 1), 4);

                        dd.Data(1:non_empty_pheno, 1) = tmp_s(1:non_empty_pheno, 2);
                        dd.Data(1:non_empty_pheno, 2) = tmp_s(1:non_empty_pheno, 3);
                        dd.Data(end - size(keep_special, 1) + 1:end, 1:2) = keep_special;

                        dd.Data(1:size(tmp_s, 1), 3) = tmp_s(:, 4);
                        dd.Data(1:size(tmp_s, 1), 4) = tmp_s(:, 5);
                        dd.Data(end, 3:4) = {fill_select_all, 'Select All'};
                    else
                        % Regular case.
                        %Only re-draw the table if the size changes
                        %TODo put something like this this in all of the table generation
                        %sections to speed up table remaking
                        if size(dd.Data, 1) ~= (size(tmp_s, 1) + 5)
                            dd.Data = cell(size(tmp_s, 1) + 5, 4);
                        end

                        dd.Data(1:end - 5, 1) = tmp_s(:, 2);
                        dd.Data(1:end - 5, 2) = tmp_s(:, 3);
                        dd.Data(end - size(keep_special, 1) + 1:end, 1:2) = keep_special;

                        dd.Data(1:end - 5, 3) = tmp_s(:, 4);
                        dd.Data(1:end - 5, 4) = tmp_s(:, 5);
                        dd.Data(end, 3:4) = {fill_select_all, 'Select All'};
                    end
                    
                    % Find the names of the user defined surfaces
                    tmpsmpls = dd.Data(ind, 4);
                    if max(contains(fieldnames(app.data.(tmpsmpls{1}).Surfaces), 'UDS'))
                        UDStmp = fieldnames(app.data.(tmpsmpls{1}).Surfaces.UDS);
                        UDStmp = [{'Surfaces', 'none'}, UDStmp(:)'];
                        if isempty(UDStmp)
                            UDStmp = {'Surfaces', 'none'};
                        end
                    else
                        UDStmp = {'Surfaces', 'none'};
                    end
                    SurfN.Items = UDStmp(2:end);
                    SurfN.Value = UDStmp{2};
                    
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
                    ind = ~cellfun('isempty', dd.Data(:, 2));
                    dd.Data(ind, 1) = {logical(p.NewData)};
                end
            end

            function backend_wrap(app, t, UIfig, ponlyN, DensButton, Den, SurfN, ModelN, special_names)
                mod_name = ModelN.Value;
                if ~strcmp(mod_name, 'none')
                    mod_name = [Constants.other_tag, Helper.valid_var(mod_name)];
                end
                pheno = logical([t.Data{:, 1}]);
                pnames = t.Data(:, 2);
                pnames = pnames(~cellfun(@isempty,pnames));
                cspot = pheno(strcmp(pnames, 'Spots'));
                cpoly = pheno(strcmp(pnames, 'Polygons'));
                cred = pheno(strcmp(pnames, 'Regions'));
                cuds = pheno(strcmp(pnames, 'User Defined Surfaces'));
                pheno(pheno) = ~ismember(pnames(pheno), special_names);
                
                Analysis.dist_func_backend( ...
                    app, ...
                    t.Data([t.Data{:, 3}]'==1, 4), ... Samples
                    pnames(pheno), ... Phenotypes
                    'calc_spots', logical(cspot), ... Spots
                    'calc_poly', logical(cpoly), ... Polygons
                    'calc_reg', logical(cred), ... Regions
                    'calc_surf', logical(cuds), ... User Defined Surface
                    'poly_name', ponlyN.Value, ... Name of Polygon
                    'calc_dens', DensButton.Value, ... Density boolean
                    'dens_radius', Den.Value, ... Radius of density
                    'surf_name', SurfN.Value, ... Name of Surface
                    'model_name', mod_name ... Name of Model
                );
                if isa(UIfig, 'matlab.ui.Figure') && isvalid(UIfig)
%                     close(UIfig)
                end
            end

        end

        function dist_func_backend(app, smplnms, pnms, varargin)
            % DIST_FUNC_BACKEND Back-end to distance function
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smplnms - cell array of names of samples
            %   - pnms - cell array of phenotypes to which distances should
            %       be calculated. In full form.
            %
            % Key-word (Optional) Inputs:
            %   - calc_spots  - default: false - Whether to calculate
            %       distances to spots.
            %   - calc_poly   - default: false - Whether to calculate
            %       distances to polygons.
            %   - poly_name   - default: 'none' - Name of polygon to which
            %       distance should be calculated, if calc_poly is true.
            %   - calc_reg    - default: false - Whether to calculate
            %       distances to regions.
            %   - model_name  - default: 'none' - Name of model which
            %       defines regions to calculate distances to, if calc_reg 
            %       is true.
            %   - calc_surf   - default: false - Whether to calculate
            %       distances to surfaces.
            %   - surf_name   - default: 'none' - Name of surface to which
            %       distance should be calculated, if calc_surf is true.
            %   - calc_dens   - default: false - Whether to calculate
            %       density of distances of phenotypes.
            %   - dens_radius - default: 50 - Radius of density inside of
            %       which distances will be included, if calc_dens is true.
            %
            % Modifies:
            %   - app - Adds a distance to things chosen by user in samples
            %       chosen by user.
            
            defaults = struct;
            defaults.calc_spots  = false;
            defaults.calc_poly   = false;
            defaults.calc_reg    = false;
            defaults.calc_surf   = false;
            defaults.poly_name   = 'none';
            defaults.calc_dens   = false;
            defaults.dens_radius = 50;
            defaults.surf_name   = 'none';
            defaults.model_name  = 'none';

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
            
            if ~iscell(smplnms)
                smplnms = {smplnms};
            end
            if ~iscell(pnms)
                pnms = {pnms};
            end
            
            % Make sure all of the elements of samples are actually sample names
            INDEMPTY = cellfun(@isempty, smplnms);
            smplnms(INDEMPTY) = [];

            % Calculate the distance for the selected samples
            for smp_i=1:numel(smplnms)
                
                % This should't be necesary but without it I get errors..
                app.DataN.Items = fieldnames(app.data);
                app.DataN.Value = smplnms{smp_i};

                %%% Calculate the distances for the different cases
                % if the user selected anything
                vPD = waitbar(0, 'Preparing to calculate things', ...
                    'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));
                vPD.Position(4) = 85;
                
                fordat = app.data.(smplnms{smp_i}).AllCells(:, {'X', 'Y', 'Z'});

                % Pull the cells to calculate distance for
                %% Calculate distance to populations
                %%%% If user selected distance to any populations
                if ~isempty(pnms)
                    % Find the names of the to data
                    ToNames = Helper.gate_full2tag(app, pnms, smplnms{smp_i});

                    % Calculate the distance to each population selected
                    Nto = numel(ToNames);
                    for i=1:Nto
                        % Calculate distance for all data
                        todatLOGIC = app.data.(smplnms{smp_i}).AllCells.(ToNames{i})==1;
%                         smplnms{smp_i}
%                         sum(todatLOGIC)
                        todat = fordat(todatLOGIC, {'X', 'Y', 'Z'});
                        % Break the data into smaller chunks so I
                        % don't break my computer
                        SzLim = 500;
                        Nloop = ceil(size(fordat, 1)/SzLim);
                        for loopi = 1:Nloop
                            if ~isvalid(vPD)
                                return;
                            end
                            st =((loopi-1)*SzLim+1);
                            en = min((st+SzLim-1), size(fordat, 1));

                            waitbar((loopi*i)/(Nloop*Nto), vPD, sprintf(['Calculating distances to populations\n' ...
                                    'Matrix# ' num2str(i) ' of ' num2str(Nto) '\n' ...
                                    num2str(3*numel(fordat.X)*numel(todat.X)) ' Calculations']));
                            %%% Do the actual distance calculation
                            if size(todat.X, 1)==1
                                Dist = sqrt( (fordat.X(st:en)-todat.X).^2 + (fordat.Y(st:en)-todat.Y).^2 + (fordat.Z(st:en)-todat.Z).^2 );
                            elseif size(todat.X, 1)>1
                                Dist = sqrt((fordat.X(st:en)' - todat.X).^2 + (fordat.Y(st:en)' - todat.Y).^2 + (fordat.Z(st:en)' - todat.Z).^2);
                                % Find the indeces of distances to
                                % "self" and calculate distance to
                                % neighbor instead of self
                                if sum(todatLOGIC(st:en))>=1
                                    DistSelf = Dist(:, todatLOGIC(st:en));
                                    DistSelf(DistSelf == 0) = NaN;
                                    Dist(:, todatLOGIC(st:en)) = DistSelf;
                                end
                                if defaults.calc_dens
                                    Dens = sum(Dist <= defaults.dens_radius);
                                end
                                Dist = min(Dist, [],'omitnan');

                                Dist = Dist';
                            else % if there are no points to calculate distance to
                                Dist = fordat.X(st:en)*0;
                            end
                            %%% add the distance as a column in the data
                            ToName = Helper.gate_tag2full(app, ToNames{i});
                            ToName = Helper.valid_var(ToName);

                            if defaults.calc_dens
                                app.data.(smplnms{smp_i}).AllCells.(Helper.valid_var([Constants.other_tag 'LocalDensityOf_' ToName]))(st:en) = Dens;
                            end
                            app.data.(smplnms{smp_i}).AllCells.(Helper.valid_var([Constants.other_tag 'DistTo_' ToName]))(st:en) = Dist;
                        end % end of loop matrix resizing
                    end
                end
                %% Find distance to spots
                if defaults.calc_spots
                    todat = table;
                    if ~isfield(app, 'points') && any(~isfield(app.points, {'PNTSx', 'PNTSy'}))
                        warndlg('No Spots found. Skipping this option.', smplnms{smp_i});
                    else
                        todat.X = app.points.PNTSx;
                        todat.Y = app.points.PNTSy;
                        if isfield(app.points, 'PNTSz')  % Sanity check/Compatibility with prev versions.
                            todat.Z = app.points.PNTSz;
                        else
                            todat.Z = 0 .* app.points.PNTSx;
                        end
                        % Calculate distance for each data
                        if ~isvalid(vPD)
                            return;
                        end
                        waitbar(0.5, vPD, 'Calculating distances to spots');

                        %%% Do the actual distance calculation (This is in 2D if spot definition is 2D)
                        if size(todat.X, 1)==1
                            if sum(todat.Z)==0
                                % Do this in 2D
                                Dist = sqrt( (fordat.X-todat.X).^2 + (fordat.Y-todat.Y).^2);
                            else
                                % Do this in 3D
                                Dist = sqrt( (fordat.X-todat.X).^2 + (fordat.Y-todat.Y).^2 + (fordat.Z-todat.Z).^2 );
                            end
                        elseif size(todat.X, 1)>1
                            if sum(todat.Z)==0
                                % Do this in 2D
                                Dist = min(sqrt( (fordat.X'-todat.X).^2 + (fordat.Y'-todat.Y).^2 ));
                            else
                                % Do this in 3D
                                Dist = min(sqrt( (fordat.X'-todat.X).^2 + (fordat.Y'-todat.Y).^2 + (fordat.Z'-todat.Z).^2 ));
                            end
                            Dist = Dist';
                        end
                        %%% add the distance as a column in the data
                        app.data.(smplnms{smp_i}).AllCells.([Constants.other_tag 'DistToSpots']) = Dist;
                    end
                end
                %% Find distance to Polygons
                if defaults.calc_poly
                    % find the names of the polygons
                    if ~isfield(app.polygons, defaults.poly_name)
                        warndlg('No Polygons with that name found. Skipping this option.', smplnms{smp_i});
                    else
                        polynms = fieldnames(app.polygons.(defaults.poly_name));
                        n = numel(polynms); % number of polygons
                        DIMPoly = 0;
                        for j=1:n
                            DIMPoly = DIMPoly + size(app.polygons.(defaults.poly_name).(polynms{j}).pos, 1);
                        end

                        %%% Do the actual distance calculation
                        %%%%% Calculate the distance in 3D because cross product %%%%%
                          for k=1:n
                              % pull the position of the nth plygon
                              POLY = app.polygons.(defaults.poly_name).(polynms{k});
                              todat = table;
                              todat.X = [POLY.pos(:, 1); POLY.pos(1, 1)];
                              todat.Y = [POLY.pos(:, 2); POLY.pos(1, 2)]; % Polygons does not work in 3D yet
                              DistLocal = zeros(size(fordat.X, 1), numel(todat.X));
                              Nlines = 1:(numel(todat.X));
                              for l=2:(numel(todat.X))
                                    if ~isvalid(vPD)
                                        return;
                                    end
                                    waitbar((k*l)/(n*numel(Nlines)), vPD, 'Calculating distances to polygons');
                                    % pt should be nx3
                                    % v1 and v2 are vertices on the line (each 1x3)
                                    % dil is a nx1 vector with the orthogonal distances
                                    pt = [fordat.X, fordat.Y, 0.*fordat.Z];  % Ignore the Z axis because polygons are 2D
                                    v1 = [todat.X(l-1), todat.Y(l-1), 0]; % Start of the line
                                    v2 = [todat.X(l), todat.Y(l), 0];     % End of the line
                                    v1 = repmat(v1,size(pt,1),1);
                                    v2 = repmat(v2,size(pt,1),1);
                                    a = v1 - v2;
                                    b = pt - v2;
                                    % dsp = distance between the start point of the line and the cell
                                    dsp = sqrt(sum((pt-v1).^2,2));
                                    % dep = distance between the end point of the line and the current point
                                    dep = sqrt(sum((pt-v2).^2,2));
                                    % dse = distance between the end point and the start point
                                    dse = sqrt(sum((v2-v1).^2,2));
                                    % Distance between the cell and the infinite line
                                    dil = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));
                                    %Type1: sqrt(dse^2+dep^2)>=dsp, if Type 1 Dist=dep
                                    ind1 = sqrt(dse.^2+dep.^2)<=dsp;
                                    DistLocal(:,l-1) = DistLocal(:,l-1) + dep.*ind1;
                                    %Type2: sqrt(dse^2+dsp^2)>=dep, if Type 2 Dist = dsp
                                    ind2 = sqrt(dse.^2+dsp.^2)<=dep;
                                    DistLocal(:,l-1) = DistLocal(:,l-1) + dsp.*ind2;
                                    %Type3: all other points
                                    ind3 = ind1==0 & ind2==0;
                                    DistLocal(:,l-1) = DistLocal(:,l-1) + dil.*ind3;
                              end % end local polygon distance loop
                              DistLocal(:,sum(DistLocal)==0) = [];
                              DistLocal = min(DistLocal,[],2);
                              % Find the points inside/outside the polygons
                              in = inpolygon(fordat.X,fordat.Y,todat.X,todat.Y);
                              DistLocal(~in) = -1.*DistLocal(~in);
                              if k==1
                                  Dist = DistLocal;
                              else
                                  Dist = -1.*min([-1.*Dist, -1.*DistLocal], [], 2);
                              end
                          end
                        % add the distance as a column in the data
                        app.data.(smplnms{smp_i}).AllCells.([Constants.other_tag 'DistTo' defaults.poly_name '_Polygon']) = Dist;
                    end
                end
                %%%
                %% Find distance to Regions
                if logical(defaults.calc_reg)
                    
                    switch app.net.(defaults.model_name((numel(Constants.other_tag)+1):end)).userdata.DataType
                        case 'Raster Scanned Neighborhoods'
                            type = 'MFIRSN';
                        case 'Cell Centered Neighborhoods'
                            type = 'MFICCN';
                        case 'Individual Cells'
                            type = 'AllCells';
                    end
                    
                    if ~isfield(app.data.(smplnms{smp_i}), type) || ...
                            ~contains(defaults.model_name, fieldnames(app.data.(smplnms{smp_i}).(type)))
                        warndlg('No Neighborhoods/Models with this name found. Skipping Regions option.', smplnms{smp_i});
                    else
                        ROWDat = app.data.(smplnms{smp_i}).(type).(defaults.model_name);
                        for clist_i=0:max(ROWDat) % for each region
                            if ~isvalid(vPD)
                                return;
                            end
                            waitbar((clist_i)/(max(ROWDat)), vPD, ['Calculating distances to Region: ' num2str(clist_i) ' sample: ' smplnms{smp_i}]);
                            % Pull out the x and y data from the
                            % cellularity table corresponding to the
                            % current region
                            todat = app.data.(smplnms{smp_i}).(type)(ROWDat==clist_i, {'X', 'Y', 'Z'});
                            % Break the matrices up a little bit for memory
                            SzLim = 500;
                            Nloop = ceil(size(fordat, 1)/SzLim);
                            for loopi = 1:Nloop
                                st =((loopi-1)*SzLim+1);
                                en = min((st+SzLim-1), size(fordat, 1));
                                %%% Do the actual distance calculation
                                if size(todat.X, 1)==1
                                    Dist = sqrt( (fordat.X(st:en)-todat.X).^2 + (fordat.Y(st:en)-todat.Y).^2 + (fordat.Z(st:en)-todat.Z).^2 );
                                elseif size(todat.X, 1)>1
                                    Dist = min(sqrt( (fordat.X(st:en)'-todat.X).^2 + (fordat.Y(st:en)'-todat.Y).^2 + (fordat.Z(st:en)'-todat.Z).^2 ));
                                    Dist = Dist';
                                elseif isempty(todat)
                                    % There are no neighborhoods of this cluster
                                    if clist_i==0
                                        break
                                    else
                                        Dist = nan(size(fordat.X(st:en)));
                                    end
                                end
                                % add the distance as a column in the data
                                model_name = defaults.model_name;
                                if startsWith(model_name, Constants.other_tag)
                                    model_name = model_name(numel(Constants.other_tag) + 1:end);
                                end
                                app.data.(smplnms{smp_i}).AllCells.(Helper.valid_var([Constants.other_tag 'DistToReg_' num2str(clist_i) 'Model_' model_name]))(st:en) = Dist;
                            end
                        end
                    end
                end
                %% Find distance to user defined surfaces
                if  defaults.calc_surf
                    if ~isfield(app.data.(smplnms{smp_i}), 'Surfaces') || ...
                            ~isfield(app.data.(smplnms{smp_i}).Surfaces, 'UDS') || ...
                            ~isfield(app.data.(smplnms{smp_i}).Surfaces.UDS, defaults.surf_name)
                        warndlg('User Defined Surfaces with that name not found. Skipping this option.', smplnms{smp_i});
                    else
                        surf = app.data.(smplnms{smp_i}).Surfaces.UDS.(defaults.surf_name).Surf;
                        if ~isvalid(vPD)
                            return;
                        end
                        waitbar(0.5, vPD, 'Calculating distances to Surfaces');
                        % Break the matrices up a little bit for memory
                        SzLim = 500;
                        Nloop = ceil(size(fordat, 1)/SzLim);
                        for loopi = 1:Nloop
                            st =((loopi-1)*SzLim+1);
                            en = min((st+SzLim-1), size(fordat, 1));

                            % Do the actual distance calculation
                            % account for 2D or 3D surfaces
                            if size(surf.Points,2)==3
                                % 3D
                                [~,Dist] = nearestNeighbor(surf,fordat.X(st:en),fordat.Y(st:en),fordat.Z(st:en));
                                % Find the points that are in the surface
                                in = inShape(surf,fordat.X(st:en),fordat.Y(st:en),fordat.Z(st:en));
                            elseif size(surf.Points,2)==2
                                % 2D
                                [~,Dist] = nearestNeighbor(surf,fordat.X(st:en),fordat.Y(st:en));
                                % Find the points that are in the surface
                                in = inShape(surf,fordat.X(st:en),fordat.Y(st:en));
                            end

                            % Find the points that are in the surface
                            Dist(~in) = -1.*Dist(~in);

                            % Add the distance as a column in the data
                            app.data.(smplnms{smp_i}).AllCells.([Constants.other_tag 'DistToUDS_' defaults.surf_name])(st:en) = Dist;
                        end
                    end
                end
                Helper.func_closeVPD(app, vPD)
            end
            
            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end
        end

        function gen_neigh_func(app, web)
            % CELL_CTR_FUNC Front-End to Cell Centered Neighborhoods Calculation
            % Makes table for user to choose on which samples to calculate
            % cell centered scan on.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %
            % Modifies:
            %   - app - Adds a distance to things chosen by user in samples
            %       chosen by user.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

            if ~Helper.any_sample(app)
                return;
            end
            %% Make some options

            tmp = Helper.populate_table(app);

            tData = cell(size(tmp, 1) + 1, 4);
            tData(1:size(tmp, 1), 1) = tmp(:, 2);
            tData{end, 1} = true;
            tData(1:size(tmp, 1), 2) = tmp(:, 3);
            tData{end, 2} = 'Select All';
            tData(1:size(tmp, 1), 3) = tmp(:, 4);
            tData{end, 3} = false;
            tData(1:size(tmp, 1), 4) = tmp(:, 5);
            tData{end, 4} = 'Select All';

            % Build the clustering options menus
            UIfig = uifigure('Name', 'Neighborhood Calculation Options', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end

            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 800 800];
            % Create the table of options
            t = uitable(UIfig);
            t.Data = tData;
            t.Position = alpha*[0 80 940 730];
            t.ColumnName = {'Include', 'Phenotypes to Center Neighborhoods on', 'Include', 'Sample'};
            t.ColumnEditable = [false false true false];
            t.ColumnWidth = {alpha*125, alpha*300, alpha*100, alpha*275};
            t.CellEditCallback = @(dd, p) switched_smpl(app, dd, p);
            
            %% Options common to both types of neighborhoods

            % Select Neighborhood type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Neighborhood Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10 60 400 15];
            NeighType = uidropdown(UIfig);
            NeighType.Items = {'Raster Scanned Neighborhoods', 'Cell Centered Neighborhoods'};
            NeighType.Value = NeighType.Items{1};
            NeighType.Position = alpha*[10 30 175 30];
            Helper.func_SetCLR(app, NeighType, 'button')
            
            % Make a radius of neighborhood numeric option
            lbl = uilabel(UIfig); lbl.Text = 'Neighborhood Radius:';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+175+5 60 130 15];
            edt = uieditfield(UIfig,'numeric');
            edt.Value = app.rwindowCCN;
            edt.Position = alpha*[10+175+5 30 110 25];
            
            
            %% Options if Raster Scanned Neighborhoods are selected
            
            % Create a loop or vectorized raster scan button
            bgRSN = uibuttongroup(UIfig, 'Visible','off');
            bgRSN.Position = alpha*[10+175+130, 25, 80, 50];
            Helper.func_SetCLR(app, bgRSN, 'UICpopup')
            % Create two radio buttons in the button group.
            r1 = uiradiobutton(bgRSN, 'Text','Fast Way',...
                'Position',alpha*[5, 25, 75, 20]);
            Helper.func_SetCLR(app, r1, 'label')
            r2 = uiradiobutton(bgRSN, 'Text','Slow Way',...
                'Position',alpha*[5, 5, 75, 20]);
            Helper.func_SetCLR(app, r2, 'label')
            
            bgRSN.Visible = 'on';
            
            %% Options if Cell-Centered Neighborhoods are selected

            % Make batch size option
            batchlbl = uilabel(UIfig); 
            batchlbl.Text = 'Batch Size (Memory Usage):';
            Helper.func_SetCLR(app, batchlbl, 'label')
            batchlbl.Position = alpha*[10+175+130 60 160 15];
            batchopt = uieditfield(UIfig,'numeric');
            batchopt.Value = get_batch(app);
            batchopt.Position = alpha*[10+175+130 30 140 25];
            
            batchlbl.Visible = 'off';
            batchopt.Visible = 'off';

%             NeighType.ValueChangedFcn = @(btn, event) switched_smpl(app, t, NeighType);
            %% functions for backend
            
            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) backend_wrap(app, t, edt, bgRSN, batchopt, NeighType.Value));
            btn.Position = alpha*[800-115, 25, 75, 50];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')
            
            NeighType.ValueChangedFcn = @(~, ~) ChangeDataType(app, t, NeighType, bgRSN, batchlbl, batchopt);
            
            function ChangeDataType(app, t, NeighType, bgRSN, batchlbl, batchopt)
                % changes visiblity of manual choice numerical input.
                p.Indices = [find(strcmp(t.Data(:, 2), 'Select All')) 1];
                p.EditData = 0;
                switch NeighType.Value
                    case 'Raster Scanned Neighborhoods'
                        batchlbl.Visible = 'off';
                        batchopt.Visible = 'off';
                        bgRSN.Visible = 'on';
                        t.ColumnEditable = [false false true false];
                        t.Data{end, 1} = true;
                    case 'Cell Centered Neighborhoods'
                        batchlbl.Visible = 'on';
                        batchopt.Visible = 'on';
                        bgRSN.Visible = 'off';
                        t.ColumnEditable = [true false true false];
                end
                switched_smpl(app, t, p);
            end
            
            function switched_smpl(app, dd, p)
                if p.Indices(2) == 3
                    % Make sure that at least one thing is selected
                    if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 3)), 3)))
                        dd.Data{p.Indices(1), 3} = true;
                        return;
                    end
                    % Process Select All button
                    fill_select_all = false;
                    if p.Indices(1) == find(strcmp(dd.Data(:, 4), 'Select All'))
                        ind = ~cellfun('isempty', dd.Data(:, 3));
                        dd.Data(ind, 3) = {logical(p.NewData)};
                        % Make sure that at least one thing is selected
                        if ~logical(p.NewData)
                            dd.Data{strcmp(dd.Data(:, 4), app.DataN.Value), 3} = true;
                        end
                        fill_select_all = logical(p.NewData);
                    end
                    ind = ~cellfun('isempty', dd.Data(:, 4));
                    ind(end) = false;
                    smpls = dd.Data(ind, 4);
                    ind = cell2mat(dd.Data(ind, 3));
                    smpls = smpls(ind);
                    tmpDat = cell(size(dd.Data, 1) - 1, 5);
                    tmpDat(:, 1) = {1};
                    tmpDat(:, 2:5) = dd.Data(1:end-1, :);

                    tmpDat = Helper.populate_table(app, 'smpls', smpls, 'prev_table', tmpDat);

                    if size(tmpDat, 1) + 1 ~= size(dd.Data, 1)
                        dd.Data = cell(size(tmpDat, 1) + 1, 4);
                    end
                    dd.Data(1:size(tmpDat, 1), 1) = tmpDat(:, 2);
                    dd.Data{end, 1} = false;
                    dd.Data(1:size(tmpDat, 1), 2) = tmpDat(:, 3);
                    dd.Data{end, 2} = 'Select All';
                    dd.Data(1:size(tmpDat, 1), 3) = tmpDat(:, 4);
                    dd.Data{end, 3} = fill_select_all;
                    dd.Data(1:size(tmpDat, 1), 4) = tmpDat(:, 5);
                    dd.Data{end, 4} = 'Select All';
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
                    ind = ~cellfun('isempty', dd.Data(:, 2));
                    dd.Data(ind, 1) = {logical(dd.Data{end, 1})};
                end
            end

            function backend_wrap(app, t, edt, bgRSN, batchopt, NeighType)
                tmpData = t.Data(1:end-1, :);
                switch NeighType
                    case 'Raster Scanned Neighborhoods'
                        app.rwindowRSN = edt.Value;
                        Analysis.raster_scan_backend( ...
                            app, ...
                            tmpData([tmpData{:,3}]==1, 4), ...   Samples
                            edt.Value, ...                       neighborhood radius
                            strcmp(bgRSN.SelectedObject.Text, 'Slow Way') ... Vectorization option
                        );
                    case 'Cell Centered Neighborhoods'
                        app.rwindowCCN = edt.Value;
                        Analysis.cell_ctr_backend( ...
                            app, ...
                            tmpData([tmpData{:,3}]==1, 4), ... Samples
                            tmpData([tmpData{:,1}]==1, 2), ... Phenotypes
                            edt.Value, ...                     neighborhood radius
                            batchopt.Value ...
                        );
                end
            end

            function batch_size = get_batch(app)
                max_all = 0;
                smpls = fieldnames(app.data);
                for smpl_idx=1:numel(smpls)
                    smpl = smpls{smpl_idx};
                    all_size = size(app.data.(smpl).AllCells, 1);
                    if all_size > max_all
                        max_all = all_size;
                    end
                end

                mem_use = max_all * 3 * 8;
                mem_free = Helper.get_free_memory();

                batch_size = 0.75 * mem_free / mem_use;
                batch_size = max(floor(batch_size), 1);
            end
        end

        function raster_scan_backend(app, smplnms, radius, looped)
            % RASTER_SCAN_BACKEND Back-End to Raster Scanned Neighborhoods Calculation
            % Performs calculation of raster scan and adds it to samples.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smplnms - Names of samples to perform raster scan on
            %   - looped - Whether to perform looped version of raster scan
            %       (slow, low memory usage), or vectorized (fast, high
            %       memory usage).
            %
            % Modifies:
            %   - app - Adds a MFIRSN field app.data.(sample name) for each
            %       of the sample, which name was given.
%%%%%
t0 = tic();
%%%%%
            for i = 1:numel(smplnms)
                %% Select sample
                app.DataN.Value = smplnms{i};

                %% Do some neighborhood analysis
                vPD = waitbar(0, 'Calculating # of cells / neighborhood', ...
                    'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));
                vPD.Position(4) = 85;
                
                clear CellularityM cells cellN x y Voxel voxels Logical datCell
                MinX = round(min(app.data.(smplnms{i}).AllCells.X));
                MaxX = round(max(app.data.(smplnms{i}).AllCells.X));
                DStepX = radius/2;
                MinY = round(min(app.data.(smplnms{i}).AllCells.Y));
                MaxY = round(max(app.data.(smplnms{i}).AllCells.Y));
                DStepY = radius/2;
                MinZ = round(min(app.data.(smplnms{i}).AllCells.Z));
                MaxZ = round(max(app.data.(smplnms{i}).AllCells.Z));
                DStepZ = radius/2;

                pnms = Helper.get_gate_tags(app, smplnms{i});
                % Pull just the channel names
                chnms = Helper.get_channels(app, smplnms{i});
                % Don't include position of cells
                chnms = chnms(~ismember(chnms, {'X', 'Y', 'Z'}));
                % Include distances etc. in the neighborhoods
                otnms = Helper.get_others(app, smplnms{i});
                datMFI = app.data.(smplnms{i}).AllCells(:,[chnms, otnms]);
                
                % Exclude strings Arrays
                INDstr = cellfun(@ischar,table2cell(datMFI),'un',0);
                INDstr = ~any(cell2mat(INDstr));
                datMFI = datMFI(:, INDstr);
                
                chnms = datMFI.Properties.VariableNames;
                % If the there is an AllCells "0" level gate, remove it (this is the same as AllCells for plotting)
                pnms = pnms(~strcmp(pnms, [Constants.gate_tag '0']));
                % Pull just the data's binary cell type annotations:
                datCell = app.data.(smplnms{i}).AllCells(:, pnms);

                %% Deal with NaN values (treat them as 0, i.e. non-contributing
                IND_nan = isnan(table2array(datCell));
                % Pull the columns with NaNs
                IND_nancol = find(any(IND_nan));
                if ~isempty(IND_nancol)
                    for col_i = IND_nancol
                        datCell{IND_nan(:, col_i), col_i} = 0;
                    end
                end
                
                IND_nan = isnan(table2array(datMFI));
                % Pull the columns with NaNs
                IND_nancol = find(any(IND_nan));
                if ~isempty(IND_nancol)
                    for col_i = IND_nancol
                        datMFI{IND_nan(:, col_i), col_i} = 0;
                    end
                end
                
                %% If there isn't any thickness to the tissue, Set z=MaxZ
                dsphere = 0;
                if MinZ==MaxZ
                    zRange = MaxZ;
                    Volume = pi*(radius^2);
                    V = 'Neigh_Area';
                    % Convert from um^2 to mm^2
                    Volume = Volume/(1000^2);
                elseif (MaxZ-MinZ) < radius
                    % If there isn't substantial thickness to the tissue
                    zRange = ((MaxZ-MinZ)/2) + MinZ; % Put Z in the middle of the tissue
                    % Call it a column
                    Volume = pi*(radius^2)*(MaxZ-MinZ);
                    V = 'Neigh_Volume';
                    % Convert from um^3 to mm^3
                    Volume = Volume/(1000^3);
                elseif (MaxZ-MinZ) >= radius
                    % If there is thickness to the tissue
                    dsphere = 1;
                    zRange = MinZ:DStepZ:MaxZ;
                    Volume = (4/3)*pi*(radius^3);
                    Volume = Volume/(1000^3);
                    V = 'Neigh_Volume';
                end
                %% Calculate CellularityM in 3D (In for loops)
                if ~looped
                    try
                        %% Calculate MFIRSN in 3D (vectorized) x,y,z vectors are all different sizes                    
                        % Figure out the ranges so number along each axis is the same
                        % Step through x positions
                        x = MinX:DStepX:MaxX;
                        % Step through Y positions
                        y = MinY:DStepY:MaxY;
                        % Step through z positions
                        z = zRange;
                        %%% TO DO: Memory Check. Otherwise use the for loop approach (paralellized?)
                        % We might also be able to break this into chunks, like
                        % I do in claculate distance... i.e. for it in loops of
                        % ~100 rows of MFIRSN at a time
                        [x, y, z] = meshgrid(x, y, z);
                        x = reshape(x, 1, []);
                        y = reshape(y, 1, []);
                        z = reshape(z, 1, []);

                        %Find the number of voxels you are subdividing your image into
                        voxels = size(x, 2);
                        %Pre-allocate your arrays
                        MFIRSN = zeros(voxels, (numel(pnms)+numel(chnms)+6));

                        % Find cellularity in a running circle
                        if ~isvalid(vPD)
                            return;
                        end
                        waitbar(0.5, vPD, sprintf(['Calculating # of cells / neighborhood for\n' strrep(smplnms{i}, '_', ' ')]));
                        if dsphere==1
                            % Find the distance between neighborhoods and all cells in the sphere
                            d = (app.data.(smplnms{i}).AllCells.X-x).^2 + (app.data.(smplnms{i}).AllCells.Y-y).^2 + (app.data.(smplnms{i}).AllCells.Z-z).^2;
                        else
                            % Find the distance between neighborhoods and all cells in a circle or cyllinder
                            d = (app.data.(smplnms{i}).AllCells.X-x).^2 + (app.data.(smplnms{i}).AllCells.Y-y).^2;
                        end
                        Logical = d < radius^2;
                        % Make a matrix with all of the cells/MFI
                        datCell = Logical' * table2array(datCell);
                        % add the number of cells of each phenotype into the
                        % main table
                        MFIRSN(:,numel(chnms)+6+1:end) = datCell;
                        datMFI = Logical' * table2array(datMFI);
                        MFIRSN(:,(6+1):(numel(chnms)+6)) = datMFI;

                        % Put in the x, y, z, volume, and number of cells (ignore Veff for now) values
                        MFIRSN(:,1:6)  = [x', y', z', Volume.*ones(numel(x),1), Volume.*ones(numel(x),1), sum(Logical)'];
                    catch
                        warning('This data is too big to do this the fast way, attempting to do it the slow way.');
                        looped = true;
                    end
                end
                %% If the vectorized attempt failed run the raster scan in a loop
                if looped
                    %Find the number of voxels you are subdividing your image into
                    voxels = numel(MinY:DStepY:MaxY)*numel(MinX:DStepX:MaxX)*numel(zRange);
                    %Pre-allocate your arrays
                    MFIRSN = zeros(voxels, (numel(pnms)+numel(chnms)+6));
                    % Pull out the position of your cells
                    Position = table2array(app.data.(smplnms{i}).AllCells(:, {'X', 'Y', 'Z'}));
                    n=0;
                    % Step through z positions if there is z data
                    for z = zRange
                        % Step through x positions
                        for x = MinX:DStepX:MaxX
                            %Step through Y positions
                            for y = MinY:DStepY:MaxY
                                n=n+1;
                                % Find cellularity in a running circle
                                if ~isvalid(vPD)
                                    return;
                                end
                                waitbar((0.95*(n/voxels)), vPD, sprintf(['Calculating # of cells / neighborhood for\n' strrep(smplnms{i}, '_', ' ')]));

                                if dsphere==1
                                    % Find the cells in the sphere;
                                    Logical = sum((Position-[x, y, z]).^2, 2) < radius^2;
                                else
                                    % Find the cells in a circle or cyllinder sphere;
                                    Logical = sum((Position(:, 1:2)-[x, y]).^2, 2) < radius^2; 
                                end

                                % Calculate the volume of the sphere that is inside the tissue
                                if MinZ==MaxZ
                                    [~, Veff] = boundary(Position(Logical, 1:2), 0);
                                    Veff = Veff/(1000^2);
                                elseif (MaxZ-MinZ) < radius  % Reduntant with next if
                                    [~, Veff] = boundary(Position(Logical, :), 0);
                                    Veff = Veff/(1000^3);
                                elseif (MaxZ-MinZ) >= radius
                                    [~, Veff] = boundary(Position(Logical, :), 0);
                                    Veff = Veff/(1000^3);
                                end
                                if sum(Logical)<3
                                    Veff = Volume;
                                end
                                if sum(Logical)==0      % There are 0 cells in this voxel do nothing
                                elseif sum(Logical)==1  % There is only one cell in this voxel
                                    MFIRSN(n,(7+numel(chnms)):end) = table2array(datCell(Logical, :));
                                    MFIRSN(n,7:(6+numel(chnms))) = table2array(datMFI(Logical, :));
                                else
                                    MFIRSN(n,(7+numel(chnms)):end) = sum(table2array(datCell(Logical, :)));
                                    % Calculate the total fluorescent intensity across all cells in this neighborhood
                                    MFIRSN(n,7:(6+numel(chnms))) = mean(table2array(datMFI(Logical, :)));
                                end
                                %Put in the x, y values into the nth row
                                MFIRSN(n,1:6)  = [x, y, z, Volume, Veff, sum(Logical)];
                            end
                        end
                    end
                end
                %% Add Cellularity M back to the data structure
                MFIRSN = array2table(MFIRSN);
                MFIRSN.Properties.VariableNames = [{'X', 'Y', 'Z', V, ['Effective_' V], 'NCells'}, chnms, pnms];
                app.data.(smplnms{i}).MFIRSN = MFIRSN;
                app.data.(smplnms{i}).ScanType = 'MFIRSN';
                
                Helper.func_closeVPD(app, vPD);
                
%%%%%
'Elapsed time define neighborhoods';
b = toc(t0);
%%%%%
            end
            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end
        end

        function cell_ctr_backend(app, smplnms, pnms, radius, batch_size)
            % CELL_CTR_BACKEND Back-End to Cell Centered Neighborhoods Calculation
            % Performs calculation of cell centered scan 
            % and adds it to samples.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smplnms - Names of samples to perform cell centered 
            %       scan on
            %   - pnms - Names of phenotypes to perform cell centered 
            %       scan on
            %   - radius - Radius of the neighborhood around the cell.
            %
            % Modifies:
            %   - app - Adds a MFICCN field app.data.(sample name) for each
            %       of the sample, which name was given.
%%%%%
t0 = tic;
%%%%%
            batch_size = floor(batch_size);
            if batch_size < 1
                return;
            end

            for i = 1:numel(smplnms)
                %% Select sample and pull the data
                vPD = waitbar(0, ['Loadinng ' smplnms{i} ' ...'], ...
                    'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));
                vPD.Position(4) = 85;

                app.DataN.Value = smplnms{i};
                Phenotypes = Helper.gate_full2tag(app, pnms, smplnms{i});

                pnms_sub = Helper.get_gate_tags(app, smplnms{i});
                % Pull just the channel names
                chnms_sub = Helper.get_channels(app, smplnms{i});
                % Don't include position of cells
                chnms_sub = chnms_sub(~ismember(chnms_sub, {'X', 'Y', 'Z'}));
                % Include distances etc. in the neighborhoods
                otnms = Helper.get_others(app, smplnms{i});
                datMFI = app.data.(app.DataN.Value).AllCells(:,[chnms_sub, otnms]);
                
                % Exclude strings Arrays
                INDstr = cellfun(@ischar,table2cell(datMFI),'un',0);
                INDstr = ~any(cell2mat(INDstr));
                datMFI = datMFI(:, INDstr);
                
                chnms_sub = datMFI.Properties.VariableNames;
                % If the there is an AllCells "0" level gate, remove it (this is the same as AllCells for plotting)
                pnms_sub = pnms_sub(~strcmp(pnms_sub, [Constants.gate_tag '0']));
                % Pull just the data's binary cell type annotations:
                datCell = app.data.(app.DataN.Value).AllCells(:, pnms_sub);
                
                %% Deal with NaN values (treat them as 0, i.e. non-contributing
                IND_nan = isnan(table2array(datCell));
                % Pull the columns with NaNs
                IND_nancol = find(any(IND_nan));
                if ~isempty(IND_nancol)
                    for col_i = IND_nancol
                        datCell{IND_nan(:, col_i), col_i} = 0;
                    end
                end
                
                IND_nan = isnan(table2array(datMFI));
                % Pull the columns with NaNs
                IND_nancol = find(any(IND_nan));
                if ~isempty(IND_nancol)
                    for col_i = IND_nancol
                        datMFI{IND_nan(:, col_i), col_i} = 0;
                    end
                end
                %% Do some neighborhood analysis
                if ~isvalid(vPD)
                    return;
                end
                waitbar(0, vPD, 'Calculating # of cells / neighborhood', ...
                    'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));

                % Number of header variables
                ind2 = 6;
                ind1 = ind2+1;
                    
                % Find the cell positions of interest
                % Pull the cell columns from sorted dat
                Dat_PH = app.data.(app.DataN.Value).AllCells(:, Phenotypes);
                % Find indeces of all selected cells
                selected_idxs = find(all(table2array(Dat_PH)==0,2)==0);
                %Find the number of voxels you are subdividing your image into
                n_batches = ceil(size(selected_idxs, 1) / batch_size);
                MFICCN = zeros(size(selected_idxs, 1), (numel(pnms_sub)+numel(chnms_sub)+ind2));

                Position = app.data.(app.DataN.Value).AllCells(:, {'X', 'Y', 'Z'});
                warning = false;
                % If there isn't any thickness to the tissue, Set z=MaxZ
                dsphere = 0;
                if numel(unique(Position.Z))==1
                    Volume = pi*(radius^2);
                    V = 'Neigh_Area';
                    % Convert from um^2 to mm^2
                    Volume = Volume/(1000^2);
                % If there isn't substantial thickness to the tissue
                elseif (max(Position.Z)-min(Position.Z)) < radius
                    % Call it a column
                    Volume = pi*(radius^2)*(max(Position.Z)-min(Position.Z));
                    V = 'Neigh_Volume';
                    % Convert from um^3 to mm^3
                    Volume = Volume/(1000^3);
                % If there is thickness to the tissue
                elseif (max(Position.Z)-min(Position.Z)) >= radius
                    Volume = (4/3)*pi*(radius^3);
                    Volume = Volume/(1000^3);
                    V = 'Neigh_Volume';
                    dsphere = 1;
                end
                Position = table2array(Position);

                % Step through all cells
                for batch_idx=1:n_batches
                    if ~isvalid(vPD)
                        return;
                    end
                    start = 1 + batch_size * (batch_idx - 1);
                    stop = start + batch_size - 1;
                    stop = min(stop, size(selected_idxs, 1));  % Avoid out of bounds array error.
                    selected_batch = selected_idxs(start:stop);
                    % Center your neighborhood around each selected cell
                    x=app.data.(app.DataN.Value).AllCells.X(selected_batch)';
                    y=app.data.(app.DataN.Value).AllCells.Y(selected_batch)';
                    z=app.data.(app.DataN.Value).AllCells.Z(selected_batch)';
                    
% % %                     % I can't remember why I explicitly pull this out and
% % %                     % put it in the front
% % %                     % Sometimes the Ch ID parameter has different names
% % %                     ChNms = fieldnames(app.data.(app.DataN.Value).AllCells);
% % %                     % Find if the ID is named ChchID or ChID
% % %                     IDIND = contains(ChNms, 'ChID') | contains(ChNms, 'ChchID') ...
% % %                         | contains(ChNms, 'Ch_ID') | contains(ChNms, 'chID');
% % %                     % If the data contains the cell IDs
% % %                     if sum(IDIND)~=0
% % %                         % Get the ID of the cell this neighborhood is
% % %                         % centered on
% % %                         ID = app.data.(app.DataN.Value).AllCells.(ChNms{IDIND})(selected_batch);
% % %                     else
% % %                         warning = true;
% % %                         ID = (start:stop) - 1;
% % %                     end
                    %% Find cellularity in a running circle
                    if ~isvalid(vPD)
                        return;
                    end
                    waitbar((0.95*(batch_idx/n_batches)), vPD, sprintf(['Calculating # of cells / cell neighborhood for\nSample: ' strrep(smplnms{i}, '_', ' ') '\nNeighorhood for cell #: ' num2str(batch_idx) ' of ' num2str(n_batches) ' Batches']));
                    if dsphere
                        % Find the cell distnace in the sphere
                        Logical = (app.data.(app.DataN.Value).AllCells.X-x).^2 + (app.data.(app.DataN.Value).AllCells.Y-y).^2 + (app.data.(app.DataN.Value).AllCells.Z-z).^2 < radius^2;
                    else
                        % Find the cell distance in a circle or cyllinder
                        Logical = (app.data.(app.DataN.Value).AllCells.X-x).^2 + (app.data.(app.DataN.Value).AllCells.Y-y).^2  < radius^2;
                    end
                    batch_sum = sum(Logical, 1);
                    % Calculate the volume of the sphere that is inside the tissue
                    Veff = zeros(stop - start + 1, 1);
                    small_logical = batch_sum < 3;
                    large_num = find(~small_logical);
                    Veff(small_logical) = Volume;
                    Veff(~small_logical) = arrayfun(@(x) boundary_wrap(x, Position, Logical, V), large_num);

                    sum_zero = batch_sum == 0;
                    sum_one  = batch_sum == 1;
                    sum_more = ~sum_zero & ~sum_one;
                    sum_one = find(sum_one);
                    sum_more = find(sum_more);

                    small_MFI = zeros(stop - start + 1, size(MFICCN, 2));
                    
                    % Logical = 0 - Cells
                    small_MFI(sum_zero, (ind1+numel(chnms_sub)):end) = zeros(sum(sum_zero), size(datCell, 2));

                    % Logical = 0 - MFI
                    small_MFI(sum_zero, ind1:(ind2+numel(chnms_sub))) = zeros(sum(sum_zero), size(datMFI, 2));

                    % Logical = 1 - Cells
                    in = cell2mat(arrayfun(@(x) single_dat_cell(x, datCell, Logical), sum_one, 'UniformOutput', false));
                    n_cols = size(small_MFI, 2) - (ind1+numel(chnms_sub)) + 1;
                    in = reshape(in, numel(in) / n_cols, n_cols);
                    small_MFI(sum_one, (ind1+numel(chnms_sub)):end) = in;

                    % Logical = 1 - MFI
                    in = cell2mat(arrayfun(@(x) single_dat_mfi(x, datMFI, Logical), sum_one, 'UniformOutput', false));
                    n_cols = numel(chnms_sub);
                    in = reshape(in, numel(in) / n_cols, n_cols);
                    small_MFI(sum_one, ind1:(ind2+numel(chnms_sub))) = in;

                    % Logical > 1 - Cells
                    in = cell2mat(arrayfun(@(x) mult_dat_cell(x, datCell, Logical), sum_more, 'UniformOutput', false));
                    n_cols = size(small_MFI, 2) - (ind1+numel(chnms_sub)) + 1;
                    in = reshape(in, numel(in) / n_cols, n_cols);
                    small_MFI(sum_more, (ind1+numel(chnms_sub)):end) = in;

                    % Logical > 1 - MFI
                    in = cell2mat(arrayfun(@(x) mult_dat_mfi(x, datMFI, Logical), sum_more, 'UniformOutput', false));
                    n_cols = numel(chnms_sub);
                    in = reshape(in, numel(in) / n_cols, n_cols);
                    small_MFI(sum_more, ind1:(ind2+numel(chnms_sub))) = in;

% % %                     % Put in the x, y values into the nth row
% % %                     small_MFI(1:end,1:7)  = [x', y', z', ID', ones(stop - start + 1, 1) * Volume, Veff, batch_sum'];
% % %                     MFICCN(start:stop, :) = small_MFI;

                    % Put in the x, y values into the nth row
                    small_MFI(1:end,1:ind2)  = [x', y', z', ones(stop - start + 1, 1) * Volume, Veff, batch_sum'];
                    MFICCN(start:stop, :) = small_MFI;
                end

                MFICCN = array2table(MFICCN);
% % %                 MFICCN.Properties.VariableNames = [{'X', 'Y', 'Z', 'ID', V, ['Effective_' V], 'NCells'}, chnms_sub, pnms_sub];
                MFICCN.Properties.VariableNames = [{'X', 'Y', 'Z', V, ['Effective_' V], 'NCells'}, chnms_sub, pnms_sub];
                app.data.(app.DataN.Value).MFICCN = MFICCN;
                app.data.(app.DataN.Value).ScanType = 'MFICCN';
                Helper.func_closeVPD(app, vPD);
% % %                 if warning
% % %                     warning(sprintf( [smplnms{i} '\nData does not contain cell ID Channel.']), 'Error');
% % %                 end
                
%%%%%
'Elapsed time define Cell Centered neighborhoods';
b = toc(t0);
%%%%%
            end

            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end

            % Functions for fast calculations on batches
            function single_cell = single_dat_cell(idx, datCell, Logical)
                single_cell = table2array(datCell(Logical(:, idx), :));
            end

            function single_mfi = single_dat_mfi(idx, datMFI, Logical)
                single_mfi = table2array(datMFI(Logical(:, idx), :));
            end

            function mult_cell = mult_dat_cell(idx, datCell, Logical)
                mult_cell = sum(table2array(datCell(Logical(:, idx), :)));
            end
            
            function mult_mfi = mult_dat_mfi(idx, datMFI, Logical)
                mult_mfi = mean(table2array(datMFI(Logical(:, idx), :)));
            end

            function Veff = boundary_wrap(idx, Position, Logical, V)
                if strcmp(V, 'Neigh_Area')
                    [~, Veff] = boundary(Position(Logical(:, idx), 1:2), 0);
                    Veff = Veff/(1000^2);
                elseif strcmp(V, 'Neigh_Volume')
                    [~, Veff] = boundary(Position(Logical(:, idx), :), 0);
                    Veff = Veff/(1000^3);
                end
            end
        end

        function TNN_func(app, type, web)
            % TNN_FUNC Front-end of Train Model Function (aka Classify Neighborhoods
            % into Regions function)
            %
            % Input:
            %   - app - Instance of CytoMAP
            %
            % Modifies:
            %   - app - Adds Region with name specified by user to either
            %       MFIRSN or MFICCN (also specified by user) 
            %       in samples chosen by user, and adds a model with the
            %       same name under app.data
            
            if nargin<3
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            if ~Helper.any_sample(app)
                return;
            end
            %% Build options for sorting

            % Pull the channel and sample names
            [mfi, smpl] = Helper.find_MFI(app);
            if isempty(mfi)
                mfi = 'AllCells';
%                 errordlg("In order to run this operation you must define neighborhoods on at least one of your samples.");
%                 return;
            end
            switch type
                case 'Cells'
                    mfi = 'AllCells';
            end

            % Build the clustering options menus
            UIfig = uifigure('Name', 'Sorting weights', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end

            Helper.func_SetCLR(app, UIfig, 'UIfigure');
            UIfig.Position = alpha*[10 10 1020 800];

            % Select Neighborhood type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Input Data Type:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+235+400 75 150 15];
            DataType = uidropdown(UIfig);
            % if there are any manually defined neighborhoods add the option to use them        
            DataType.Items = {'Raster Scanned Neighborhoods', ...
                            'Cell Centered Neighborhoods', ...
                            'Individual Cells'};
            switch mfi
                case 'MFIRSN'
                    DataType.Value = {'Raster Scanned Neighborhoods'};
                    UIfig.Name = 'Cluster RSN Neighborhoods';
                case 'MFICCN'
                    DataType.Value = {'Cell Centered Neighborhoods'};
                    UIfig.Name = 'Cluster CCN Neighborhoods';
                case 'AllCells'
                    DataType.Value = {'Individual Cells'};
                    UIfig.Name = 'Cluster Cells';
            end
            DataType.Position = alpha*[10+235+400 45 150 30];
            Helper.func_SetCLR(app, DataType, 'button');

            % Create the table of options
            t = uitable(UIfig);
            tmp = Helper.populate_table(app, 'smpls', {smpl}, 'MFI', mfi, 'fill_checkbox', false);
            t.Data = cell(size(tmp, 1) + 1, size(tmp, 2));
            t.Data(1:end-1, :) = tmp;
            t.Data{end, 2} = false;
            t.Data{end, 3} = 'Select All';
            t.Data{end, 5} = false;
            t.Data{end, 6} = 'Select All';
            t.Data{end, 7} = false;
            t.Data{end, 8} = 'Select All';
            t.Position = alpha*[0 90 1000 710];
            t.ColumnName = {'Weight', 'Use for sorting', 'Phenotype (names must be consistent across samples)', ...
                            'Weight', 'Use for sorting', 'Channel MFI', 'Use for sorting', 'Sample'};
            t.ColumnEditable = [true true false true true false true false];
            t.ColumnWidth = {alpha*50, alpha*100, alpha*250, alpha*50, alpha*100, alpha*150, alpha*100, alpha*150};

            % Create initial number of regions
            lbl = uilabel(UIfig); lbl.Text = 'Number of Regions:';
            lbl.HorizontalAlignment = 'right';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+235 40 150 15];
            edt = uieditfield(UIfig,'numeric','ValueChangedFcn',@(edt,event) funcNReg(app,edt));
            edt.Value = app.NReg;
            edt.Position = alpha*[400+165 40 80 20];
            edt.Visible = 'off';

            % Automatic choices for number of classes
            NofClassesAlg = uidropdown(UIfig, ...
                'ValueChangedFcn', @(NofClassesAlg, manualBtn) ...
                                    ChangeNofClasses(NofClassesAlg, edt));
            NofClassesAlg.Position = alpha*[400 35 240 25];
            NofClassesAlg.Items = {'Davies Bouldin (Default)', ...
                                   'Calinski-Harabasz', ...
                                   'Gap', ...
                                   'Manual choice of number of regions'};
            NofClassesAlg.Value = NofClassesAlg.Items{1};
            Helper.func_SetCLR(app, NofClassesAlg, 'button');

            % If there are manually define regions add the option to
            % cluster those
% % %             if any(startsWith(fieldnames(app.data.(t.Data{[t.Data{:,7}]==1, 8}).(mfi)), Constants.neigh_tag))
            DataType.Items = [DataType.Items, {'Manually Defined Regions (RSN)', 'Manually Defined Regions (CCN)'}];
% % %             end
            
            % Make a color scheme option
            lbl = uilabel(UIfig); lbl.Text = sprintf('Color scheme:');
            lbl.HorizontalAlignment = 'right';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+235 60 150 20];
            ClrMap = uieditfield(UIfig,'text');
            ClrMap.Value = 'sum(y,2)';
            ClrMap.Position = alpha*[165+235 60 200 20];

            % Select Clusterin algorithm
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Algorithm:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+235+400+150 75 150 15];
            ClustAlgorithm = uidropdown(UIfig);
            ClustAlgorithm.Items = { ...
                'NN Self Organizing Map',...
                'k-means',...
                'Gaussian Distribution Model', ...
                'DBSCAN'};
            ClustAlgorithm.Value = {'NN Self Organizing Map'};
            ClustAlgorithm.Position = alpha*[10+235+400+150 45 150 30];
            Helper.func_SetCLR(app, ClustAlgorithm, 'button');

            lbl = uilabel(UIfig);
            ShowNetBtn = uibutton(UIfig);
            lbl.Text = sprintf('Model name:');
            lbl.HorizontalAlignment = 'right';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+235 15 150 20];
            NetName = uidropdown(UIfig);
            NetName.Items = [{Constants.new_model}, Helper.full_var(fieldnames(app.net))'];
            NetName.Value = NetName.Items(1);
            NetName.Position = alpha*[165+235 10 240 25];
            Helper.func_SetCLR(app, NetName, 'button')            

            ShowNetBtn.Text = 'Show Model Configuration';
            Helper.func_SetCLR(app, ShowNetBtn, 'button')
            ShowNetBtn.Position = alpha*NetName.Position;
            ShowNetBtn.Position(1) = alpha*(ShowNetBtn.Position(1) + ShowNetBtn.Position(3) + 5);
            ShowNetBtn.Position(3) = alpha*150;
            ShowNetBtn.Visible = ~strcmp(NetName.Value, {Constants.new_model});
            ShowNetBtn.ButtonPushedFcn = @(btn, event) show_train(app, NetName.Value);
            Helper.func_SetCLR(app, ShowNetBtn, 'button')

            % Select Data Preperation
            DataPrep = uidropdown(UIfig);
            DataPrep.Items = Constants.neigh_cell_norm;
            DataPrep.Value = {'Composition: Number of Cells / Total cells in Neighborhood'};
            DataPrep.Position = alpha*[10 50 175 30];
            Helper.func_SetCLR(app, DataPrep, 'button')
            
            % Select Data Preperation for MFI
            DataPrepMFI = uidropdown(UIfig);
            DataPrepMFI.Items = Constants.neigh_mfi_norm;
            DataPrepMFI.Value = {'MFI normalized to max MFI per neighborhood'};
            DataPrepMFI.Position = alpha*[10 20 175 30];
            Helper.func_SetCLR(app, DataPrepMFI, 'button')
            
            if strcmp(DataType.Value, 'Individual Cells')
                DataPrep.Visible = 'off';
                DataPrepMFI.Items = Constants.cell_mfi_norm;
                DataPrepMFI.Value = DataPrepMFI.Items(1);
            else
                DataPrep.Visible = 'on';
                DataPrepMFI.Items = Constants.neigh_mfi_norm;
                DataPrepMFI.Value = DataPrepMFI.Items(3);
            end
            
            % Create a Norm per sample or per dataset label
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Normalize per:');
            lbl.HorizontalAlignment = 'left';
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[175+10 20+45 150 20];
            % Create a Norm per sample or per dataset button
            NormPer = uibuttongroup(UIfig, 'Visible','off');
            NormPer.Position = alpha*[175+10, 20, 100, 45];
            Helper.func_SetCLR(app, NormPer, 'UICpopup')
            % Create two radio buttons in the button group.
            r1 = uiradiobutton(NormPer, 'Text','Sample',...
                'Position',alpha*[5, 25, 75, 20]);
            Helper.func_SetCLR(app, r1, 'label')
            r2 = uiradiobutton(NormPer, 'Text','Dataset',...
                'Position',alpha*[5, 5, 75, 20]);
            Helper.func_SetCLR(app, r2, 'label')
            NormPer.Visible = 'on';

            % Create Menu bar
            FileMenu = uimenu(UIfig);
            FileMenu.Text = 'File';
            
            % Create a load model option
            loadmdl = uimenu(FileMenu);
            loadmdl.MenuSelectedFcn = @(btn,event) Analysis.Load_NN(app, NetName);
            loadmdl.Text = 'Load Model';
            
            % Create a save model option
            savemdl = uimenu(FileMenu);
            savemdl.MenuSelectedFcn = @(btn,event)  Analysis.Save_NN(app);
            savemdl.Text = 'Save Model';

            % Deal with changing data types
            NetName.ValueChangedFcn = @(dd, p) net_changed(dd, p, ShowNetBtn, app, DataType, t, UIfig, NetName.Value);
            t.CellEditCallback = @(dd, p) switched_sample(app, DataType, dd, p, UIfig, NetName.Value);
            DataType.ValueChangedFcn = @(btn, event) ChangeDataType(app, t, DataType, DataPrep, DataPrepMFI, UIfig, NetName.Value);

            % Create a plot option
            plt = uiswitch(UIfig, 'Orientation', 'vertical');
            plt.Items = {'on','Plot Output'};
            plt.Value = plt.Items{2};
            plt.Position = alpha*[10+235+400+150+150+32 25 150 30];
            plt.ValueChangedFcn = @(btn, event) SelectPlotFunction(btn, event);
 
            
            % Create a push button
            btn = uibutton(UIfig, 'push', 'ButtonPushedFcn', @(btn,event) backend_wrap( ...
                t, app, ClrMap, DataPrep, DataType, ClustAlgorithm, NofClassesAlg, NetName, NormPer, DataPrepMFI, plt.Value));
            btn.Position = alpha*[10+235+400+150 10 150 30];
            btn.Text = 'Ok';
            Helper.func_SetCLR(app, btn, 'button')
            
            function SelectPlotFunction(btn, event)
                switch event.Value
                    case 'on'
                        btn.Items = {'off','Plot Output'};
                        btn.Value = btn.Items{1};
                    case 'Plot Output'
                        btn.Items = {'on','Plot Output'};
                        btn.Value = btn.Items{2};
                end
            end
            
            function switched_sample(app, DataType, dd, p, UIfig, NetName)
                switch DataType.Value
                    case 'Manually Defined Regions (RSN)'
                        UIfig.Name = ['Cluster Custom RSN Neighborhoods: ' NetName];
                    case 'Manually Defined Regions (CCN)'
                        UIfig.Name = ['Cluster Custom CCN Neighborhoods: ' NetName];
                    case 'Raster Scanned Neighborhoods'
                        UIfig.Name = ['Cluster RSN Neighborhoods: ' NetName];
                    case 'Cell Centered Neighborhoods'
                        UIfig.Name = ['Cluster CCN Neighborhoods: ' NetName];
                    case 'Individual Cells'
                        UIfig.Name = ['Cluster Cells: ' NetName];
                end
                
                if p.Indices(2) == 7
                    switch DataType.Value
                        case 'Manually Defined Regions (RSN)'
                            scan_type = 'MFIRSN';
                        case 'Manually Defined Regions (CCN)'
                            scan_type = 'MFICCN';
                        case 'Raster Scanned Neighborhoods'
                            scan_type = 'MFIRSN';
                        case 'Cell Centered Neighborhoods'
                            scan_type = 'MFICCN';
                        case 'Individual Cells'
                            scan_type = 'AllCells';
                    end

                    % Make sure that at least one thing is selected
                    if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 7)), 7)))
                        dd.Data{p.Indices(1), 7} = true;
                        return;
                    end
                    
                    % Process Select All button
                    fill_select_all = false;
                    if p.Indices(1) == find(strcmp(dd.Data(:, 8), 'Select All'))
                        ind = dd.Data(:, 7);
                        ind = ~cellfun('isempty', ind);
                        dd.Data(ind, 7) = {logical(p.NewData)};
                        if ~logical(p.NewData)
                            dd.Data{strcmp(dd.Data(:, 8), app.DataN.Value), 7} = true;
                        end
                        fill_select_all = logical(p.NewData);
                    elseif ~isfield(app.data.(dd.Data{p.Indices(1), 8}), scan_type) && p.NewData
                        dd.Data(p.Indices(1), 7) = {false};
                        warndlg(strcat('Sample ', dd.Data(p.Indices(1), 8), ' does not have ', DataType.Value, '.', newline, ...
                            'Please run corresponding scan before choosing that sample.'));
                        return;
                    end

                    ind = ~cellfun('isempty', dd.Data(:, 7));
                    ind(end) = false;
                    ind = dd.Data(ind, 7);
                    ind = cell2mat(ind);
                    tmp_data = Helper.populate_table(app, ...
                        'smpls', dd.Data(ind, 8), ...
                        'mfi', scan_type, ...
                        'prev_table', dd.Data);
                    
                    % if nothing changed in the table don't rewrite it
                    if ~isempty(tmp_data) % Sanity check.
                        if ~Helper.setequal(dd.Data(1:end-1, :), tmp_data)
                            dd.Data = cell(size(tmp_data, 1) + 1, size(tmp_data, 2));
                            dd.Data(1:end-1, :) = tmp_data;
                            dd.Data{end, 2} = false;
                            dd.Data{end, 3} = 'Select All';
                            dd.Data{end, 5} = false;
                            dd.Data{end, 6} = 'Select All';
                            dd.Data{end, 7} = fill_select_all;
                            dd.Data{end, 8} = 'Select All';
                        end
                    end
                elseif p.Indices(2) == 2 && p.Indices(1) == size(dd.Data, 1)
                    ind = dd.Data(:, 3);
                    ind = ~cellfun('isempty', ind);
                    dd.Data(ind, 2) = {logical(p.NewData)};
                elseif p.Indices(2) == 5 && p.Indices(1) == size(dd.Data, 1)
                    ind = dd.Data(:, 5);
                    ind = ~cellfun('isempty', ind);
                    dd.Data(ind, 5) = {logical(p.NewData)};
                end
            end

            function backend_wrap(t, app, ClrMap, DataPrep, DataType, ClustAlgorithm, NofClassesAlg, NetName, NormPer, DataPrepMFI, if_plt)
                switch if_plt
                    case 'off'
                        if_plt = false;
                    case 'Plot Output'
                        if_plt = true;
                end
                tData = t.Data(1:end - 1, :);
                phenotypes    = tData([tData{:,2}]==1, 3);
                pheno_weights = tData([tData{:,2}]==1, 1);
                pheno_weights = cell2mat(pheno_weights)';
                % Find the index of the samples, which are going to be used to train the NN
                INDSmpls = [tData{1:numel(app.DataN.Items), 7}]==1;
                % Find the names of the samples
                smplnms = tData(INDSmpls, 8);
                % Find MFI Names
                MFIList     = tData([tData{:,5}]==1, 6);
                MFI_weights = tData([tData{:,5}]==1, 4);
                MFI_weights = cell2mat(MFI_weights)';
                
                % Convert to channel tags/valid ignore names
                MFIList = Helper.valid_channel(MFIList);
                MFIList(ismember(MFIList, Constants.ignore_names)) = Helper.valid_var(MFIList(ismember(MFIList, Constants.ignore_names)));
                Analysis.TNN_backend(app, ...
                    smplnms, ...
                    phenotypes, ...
                    pheno_weights, ...
                    MFIList, ...
                    MFI_weights, ...
                    ClrMap.Value, ...
                    DataPrep.Value, ...
                    DataType.Value, ...
                    ClustAlgorithm.Value, ...
                    NofClassesAlg.Value, ...
                    NetName, ...
                    NormPer.SelectedObject.Text, ...
                    DataPrepMFI.Value, ...
                    if_plt)
            end

            % Change number of regions
            function funcNReg(app,edt)
                app.NReg = edt.Value;
            end

            function ChangeDataType(app, t, DataType, DataPrep, DataPrepMFI, UIfig, NetName)
                % changes visiblity of manual choice numerical input.
                p.Indices = [1 7];
                p.EditData = 0;
                switched_sample(app, DataType, t, p, UIfig, NetName);
                if strcmp(DataType.Value, 'Individual Cells')
                    DataPrep.Visible = 'off';
                    DataPrepMFI.Items = Constants.cell_mfi_norm;
                    DataPrepMFI.Value = DataPrepMFI.Items(1);
                else
                    DataPrep.Visible = 'on';
                    DataPrepMFI.Items = Constants.neigh_mfi_norm;
                    DataPrepMFI.Value = DataPrepMFI.Items(3);
                end
            end

            function show_train(app, NetName)

                apha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

                show_fig = uifigure('Name', strcat('Training configuration for: ', NetName), 'Scrollable', 'on');

                show_fig.Position = apha*[10 10 1000 600];
                udata = app.net.(Helper.valid_var(NetName)).userdata;
                
                past_table = cell(max([numel(udata.SampleNames), numel(udata.MFI), numel(udata.pheno), 4]), 6);
                past_table(1:numel(udata.pheno), 1)       = num2cell(udata.pheno_weights)';
                past_table(1:numel(udata.pheno), 2)       = udata.pheno;
                past_table(1:numel(udata.MFI)  , 3)       = num2cell(udata.MFI_weights)';
                past_table(1:numel(udata.MFI)  , 4)       = udata.MFI;
                past_table(1:numel(udata.SampleNames), 5) = udata.SampleNames;
                past_table(1, 6) = {udata.DataPrep};
                past_table(2, 6) = {udata.DataType};
                if any(contains(fieldnames(udata), 'NofClassesAlg'))
                    past_table(3, 6) = {udata.NofClassesAlg};
                end
                if any(contains(fieldnames(udata), 'NormOpt'))
                    past_table(4, 6) = {udata.NormOpt};
                end
                
                uitable(show_fig, 'Data', past_table, 'Position', apha*[10, 10, 1000, 600]);
                
            end

            function net_changed(~, p, Btn, app, DataType, t, UIfig, NetName)
                Btn.Visible = ~strcmp(p.Value, {Constants.new_model});
                
                p2.Indices = [1 7];
                p2.EditData = 0;
                switched_sample(app, DataType, t, p2, UIfig, NetName)
            end

            function ChangeNofClasses(dd,p)
                if strcmp(dd.Value, 'Manual choice of number of regions')
                    p.Visible = 'on';
                    dd.Position = alpha*[400 35 160 25];
                else
                    p.Visible = 'off';
                    dd.Position = alpha*[400 35 240 25];
                end
            end
        end

        function TNN_backend(app, smplnms, phenotypes, pheno_weights, MFIList, MFI_weights, ClrMap, DataPrep, ...
                DataType, ClustAlgorithm, NofClassesAlg, ModelName, NormOpt, DataPrepMFI, if_plt)
            % TNN_BACKEND Back-end of Train Model Function (aka Classify Neighborhoods
            % into Regions function)
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smplnms - cell of strings/chars containing samples to 
            %       train model on.
            %   - phenotypes - cell of strings/chars containing phenotypes
            %       (short name version) to train model on.
            %   - pheno_weights - numeric array containing weights of those
            %       phenotypes. Weights denote how important given
            %       phenotype is for this model.
            %   - MFIList - cell of strings/chars containing MFIs (aka
            %       channels, ignore names etc.) (valid version) to train
            %       model on.
            %   - MFI_weights - numeric array containing weights of those
            %       MFIs. Weights denote how important given MFI is for
            %       this model.
            %   - ClrMap - matlab expresion in terms of 'y' (double array
            %       containing data from all phenotype/sample combinations
            %       chosen)
            %   - DataPrep - Type of Data Preparation Algorithm. For more
            %       details look at documentation for Helper.Func_DataPrep
            %   - DataType - Whether to pull data from:
            %       'Raster Scanned Neighborhoods',
            %       'Cell Centered Neighborhoods', "Individual Cells', or
            %       'Manually defined Neighborhoods'
            %       Those are also the only valid options.
            %   - ClustAlgorithm - Type of clustering algorithm. Possible
            %       choices are:
            %           - 'NN Self Organizing Map' - For Clustering Neural
            %               Network
            %           - 'Gaussian Distribution Model'
            %           - 'k-means'
            %   - NofClassesAlg - Type of algorithm to choose number of
            %       regions to cluster data into. It can be either:
            %           - 'Manual choice of number of regions' - For user
            %               defined number of regions.
            %           - Any from 'Davies Bouldin', 'Calinski-Harabasz',
            %               'Gap', where the following algorithm will
            %               determine number of region.
            %   - ModelName - Name of the model to be trained, and name
            %       under which regions will be saved.
            %   - NormOpt - This selects weather to normalize each sample
            %       independently or the whole dataset together
            %   - DataPrepMFI - This selects the type of normalization to
            %      apply yo the MFI data, if there is any
            %   - ifplt - This toggles weather or not to plot the resulting
            %      clustered output
            %
            % Modifies:
            %   - app - Adds Region with name specified by user to either
            %       MFIRSN or MFICCN (also specified by user) 
            %       in given samples, and adds a model with the same name 
            %       under app.data
            
%             if all(cellfun('isempty', phenotypes))
%                 errordlg('No Phenotypes Available. Aborting.');
%                 return;
%             end
            NetName = ModelName;
            ModelName = ModelName.Value;

            % If user wants to create new model ask them to name it
            if strcmp(ModelName, Constants.new_model)
                prompt = {'Enter name of new model:'};
                title = 'Name new model';
                dims = [1 35];
                definput = {strcat('Model', num2str(numel(fieldnames(app.net)) + 1))};
                ModelName = inputdlg(prompt,title,dims,definput);
                if isempty(ModelName)
                    return;
                end
            end
            
%%%%%
t0 = tic();
%%%%%

            % Make sure model is a valid variable (dropdowns use full_var version of names)
            ModelName = Helper.valid_var(ModelName);
            if iscell(ModelName)
                ModelName = ModelName{1};
            end
            if ~ismember(ModelName, fieldnames(app.net))
                app.net.(ModelName) = struct;
                new_model = true;
            else
                retrain = questdlg( ...
                    strcat('Network already found in the CytoMAP.', newline, ...
                        'Do you want to reuse or retrain it?'), ...
                    'Network Found', ...
                    'Retrain','Reuse','Cancel');
                % Handle response
                switch retrain
                    case 'Retrain'
                        % Override and reset the model
                        app.net.(ModelName) = struct;
                        new_model = true;
                    case 'Reuse'
                        new_model = false;
                    case 'Cancel'
                        return;
                    otherwise
                        return;
                end
            end
            
            NetName.Items = [{Constants.new_model}, Helper.full_var(fieldnames(app.net))'];

            vPD = waitbar(0, 'Preparing data for sorting', ...
                'CreateCancelBtn', @(h, ~) cancel_waitbar_callback(h));
            vPD.Position(4) = 85;
            %% Load and prep the data for sorting
            rmvzeros = 1;
            Man_Reg = [];
            
            if strcmp(DataType, 'Individual Cells')
                % Change the data preperation to not normalize the cell
                % number, since we don't technically have a cell number in
                % te case of individual cels
                DataPrep = 'Cellularity: Number of Cells / Neighborhood';
            elseif strcmp(DataType, 'Manually Defined Regions (RSN)')
                DataType = 'Raster Scanned Neighborhoods';
                Man_Reg = 1;
            elseif strcmp(DataType, 'Manually Defined Regions (CCN)')
                DataType = 'Cell Centered Neighborhoods';
                Man_Reg = 1;
            end
            
            [Dat_Pre, Dat_PreALL, INDDatPre, INDzrs, INDons] = Helper.func_loaddata( ...
                app, ...
                smplnms, ...
                phenotypes, ...
                MFIList, ...
                pheno_weights, ...
                MFI_weights, ...
                DataType, ...
                DataPrep, ...
                NormOpt, ...
                rmvzeros, ...
                DataPrepMFI);
            
            switch DataType
                case 'Raster Scanned Neighborhoods'
                    type = 'MFIRSN';
                case 'Cell Centered Neighborhoods'
                    type = 'MFICCN';
                case 'Individual Cells'
                    type = 'AllCells';
                    % remove the binary cell type classifiers from the data table
                    Dat_Pre = Dat_Pre(:, (numel(phenotypes)+1):end);
            end
            % if manual regions were selected pull just those
            if Man_Reg==1          
                
                % Pull all available neighborhood names
                Man_names = fieldnames(app.data.(smplnms{1}).(type));
                indnms = startsWith(Man_names, Constants.neigh_tag);
                Man_names = Man_names(indnms);
                % Make sure those names are in all samples
                for smpl_i=2:numel(smplnms)
                    % Pull all neighborhood names
                    Man_names_tmp = fieldnames(app.data.(smplnms{smpl_i}).(type));
                    indnms = startsWith(Man_names_tmp, Constants.neigh_tag);
                    Man_names_tmp = Man_names_tmp(indnms);
                    % if the sets are not equal take the more restrictive one
                    if any(~contains(Man_names, Man_names_tmp))
                        Man_names = Man_names_tmp;
                    end
                end
                % have user select manual regions
                [indx,tf] = listdlg('ListString',Man_names);
                if ~tf
                    return
                end
                
                ManRegs = Man_names(indx);
                
% %                 % Pull the manual region names
% %                 fulllist = [phenotypes; MFIList];
% %                 IND_reg = startsWith(fulllist, Constants.neigh_tag);
% %                 ManRegs = fulllist(IND_reg);
% %                 % remove the binary cell type classifiers from the data table
% %                 % This assumes they are alwayse the last rows in the table
% %                 Dat_Pre = Dat_Pre(:, ~IND_reg);
                
                nmissing = 1;
                for smpl_i=1:numel(smplnms)
                    for reg_i = 1:numel(ManRegs)
                        if reg_i==1
                            Reg_Logic = app.data.(smplnms{smpl_i}).MFIRSN.(ManRegs{reg_i});
                        else
                            reglogic = app.data.(smplnms{smpl_i}).MFIRSN.(ManRegs{reg_i});
                            Reg_Logic = Reg_Logic | reglogic;
                        end
                    end
                    % Pull the overlap of non-zero elements and region logic
                    Reg_Logic = Reg_Logic & INDons{smpl_i};
                    % set all elements not in sample to 0
                    dattmp = Dat_Pre(INDDatPre{smpl_i}{1}(1):INDDatPre{smpl_i}{1}(2),:);
%                     dattmp(~Reg_Logic(INDons{smpl_i}), :) = 0.*dattmp(~Reg_Logic(INDons{smpl_i}), :);
                    dattmp(~Reg_Logic(INDons{smpl_i}), :) = NaN;

                    Dat_Pre(INDDatPre{smpl_i}{1}(1):INDDatPre{smpl_i}{1}(2),:) = dattmp;
                    
                    % fix INDons{smpl_i}
                    INDons{smpl_i} = Reg_Logic;
                    
                    % fix INDzrs{smpl_i}
                    INDzrs{smpl_i} = ~INDons{smpl_i};
                    

                    % fix INDDatPre
                    if INDDatPre{smpl_i}{1}(1)~=1
                        nmissing = nmissing + sum(isnan(dattmp(:,1)));
                        INDDatPre{smpl_i}{1}(1) = INDDatPre{smpl_i-1}{1}(2)+1;
                        INDDatPre{smpl_i}{1}(2) = INDDatPre{smpl_i}{1}(2) - nmissing;
                    else
                        nmissing = nmissing + sum(isnan(dattmp(:,1)))-1;
                        INDDatPre{smpl_i}{1}(2) = INDDatPre{smpl_i}{1}(2) - nmissing;
                    end
                end
                % remove zeros
                Dat_Pre = Dat_Pre(~isnan(Dat_Pre(:,1)), :);
            end
            
            SortNames = Helper.gate_full2tag(app, phenotypes, smplnms{1});
            app.data.(smplnms{1}).SortNames = SortNames;

            if numel(smplnms)>1
                for i=2:numel(smplnms)
                    SortNames_i = Helper.gate_full2tag(app, phenotypes, smplnms{i});
                    app.data.(smplnms{i}).SortNames = SortNames_i;
                end
            end

            % Use this notation to pull indeces etc. begin
            % begin:  INDDatPre{1}{1}(1)
            % end:    INDDatPre{1}{1}(2)
            % Sample: INDDatPre{1}{2}
            %% Automatically decide how many regions to break tissue into
            if new_model
                if ~isvalid(vPD)
                    return;
                end
                if ~strcmp(NofClassesAlg, 'Manual choice of number of regions')
                    waitbar(0.1, vPD, 'Determining Number of Clusters');
                    
                    % Initialize the figure
                    fig = figure;
                    % add an export data option
                    % Create Menu bar
                    ExportMenu = uimenu(fig);
                    ExportMenu.Text = 'Export';
                    % Create an export data option
                    ExportPltDat = uimenu(ExportMenu);
                    ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                    ExportPltDat.Text = 'Export Plot Data to .csv';
                    
                    fig.Visible = false;
                    fig.Color = 'w';
                    
                    % Add an export data option
                    % Create Menu bar
                    ExportMenu = uimenu(fig);
                    ExportMenu.Text = 'Export';
                    % Create an export data option
                    ExportPltDat = uimenu(ExportMenu);
                    ExportPltDat.MenuSelectedFcn = @(btn,event) Plt_Helper.export_plot_data(app, 'Line_Point_Plot');
                    ExportPltDat.Text = 'Export Plot Data to .csv';
                    
                    % Determine the number of clusters
                    if strcmp(ClustAlgorithm, 'DBSCAN')
                        NR = 2;
                        NRegCValues = 2;
                    else
                        NCAlg = strrep(strrep(strrep(NofClassesAlg, ' (Default)', ''), ' ', ''), '-', '');
                        NR = evalclusters(Dat_Pre, @(data, idx) Helper.func_cluster(idx, ClustAlgorithm, data, []), NCAlg, 'KList', 1:24);
                        figure(fig)
                        plot(NR.CriterionValues, 'b-o');
                        fig.Visible = true;
                        xlabel('Number of Clusters');
                        ylabel([NofClassesAlg ' values']);
                        axis square
                        box off
                        NRegCValues = NR.CriterionValues;
                        NR = NR.OptimalK;
                    end
                    
                    if isnan(NR)
                        NR = 2;
                    end
                else
                    NRegCValues = app.NReg;
                    NR = app.NReg;
                end
             
                app.net.(ModelName).name = ModelName; % Seems redundant, but needed for loading and saving network
                app.net.(ModelName).cmap = jet(NR + 1);
                app.net.(ModelName).cmap(1, :) = [1,1,1];
                app.net.(ModelName).NR = NR;
            else
                NR = app.net.(ModelName).NR;
            end
            if ~isvalid(vPD)
                return;
            end
            
%%%%%
'Elapsed time to determine number of clusters';
b = toc(t0);

t0 = tic();
%%%%%
            waitbar(0.15, vPD, 'Clustering Data');
            %% Sort the data with the specified algorithm

            % Pull out any NaN Rows
            INDNaN = sum(isnan(Dat_Pre), 2)~=0;

            % Initialize cluster vector
            ROWAllDat = zeros(size(Dat_Pre, 1), 1);
                
            if ~isvalid(vPD)
                return;
            end
            
            if new_model
                waitbar(0.2, vPD, 'Training New Model');
                app.net.(ModelName).type = ClustAlgorithm;
                [ROWAllDat(~INDNaN), model] = ...
                    Helper.func_cluster(NR, ClustAlgorithm, Dat_Pre(~INDNaN, :), []);
                app.net.(ModelName).Network = model;
            else
                waitbar(0.2, vPD, 'Clustering Data');
                [ROWAllDat(~INDNaN), ~] = ...
                    Helper.func_cluster(NR, ClustAlgorithm, Dat_Pre(~INDNaN, :), app.net.(ModelName).Network);
            end
            
            if strcmp(ClustAlgorithm, 'DBSCAN')
                ROWAllDat(ROWAllDat==-1)=0;
                ROWAllDat(~INDNaN) = ROWAllDat(~INDNaN)+1;
                NR = numel(unique(ROWAllDat));
                
                app.net.(ModelName).cmap = jet(NR + 1);
                app.net.(ModelName).cmap(1, :) = [1,1,1];
                app.net.(ModelName).NR = NR;
                
            end
            
            if ~isvalid(vPD)
                return;
            end
            
            ROW = cell(numel(smplnms), 1);
            % Make a classifier matrix for each data
            for i=1:numel(smplnms)
                % Pull out the classification data for each sample
                ROW{i} = ROWAllDat(INDDatPre{i}{1}(1):INDDatPre{i}{1}(2));

                % If there were no 0 elements in CellularityM then ROW
                % and CellularityM will have the same length and ROW
                % can just be tacked on as a new column

                % Pull the number of elements
                NNeighbor = size(app.data.(smplnms{i}).(type), 1);

                if NNeighbor==size(ROW{i}, 1)
                    Class=ROW{i};

                % If ROW{i} and NNeighbor do not have the same
                % length padd with zeros where there were zero elements
                % in CellularityM
                else
                    Class = zeros(NNeighbor, 1);
                    Class(INDzrs{i})=0.*Class(INDzrs{i});
                    Class(INDons{i})=ROW{i};

                end

                % Redefine ROW so it is the proper length size, etc.
                ROW{i} = Class;
            end
            % Define the plot title
            ttl = cell(numel(smplnms), 1);

            %% re-order color scheme based on cellularity (Do this only for the first sample)
            ROW2 = ROW{1}(ROW{1}~=0);
            y = Dat_PreALL(ROW{1}~=0, :);
            % Use this to determine index of colorscheme
% %                 y = table2array(y);
            Selected = find(ismember(Dat_PreALL.Properties.VariableNames, [SortNames; MFIList]));
            y = table2array(y(:,Selected));
            y = eval([ClrMap ';']);
            
            [~, IND] = sort(y, 'ascend');
            ORDER = unique(ROW2(IND), 'stable');
            % Add back in any regions not in this sample
            ORDER = unique([ORDER' 1:NR], 'stable');
            RowOrder = zeros(1,2);
            % do this for 1 to the number of regions
            for i=1:NR
                if ORDER(i) ~= i
                    RowOrder = [RowOrder; i, ORDER(i)];
                    ROW2(ROW2==i) = nan;
                    ROW2(ROW2==ORDER(i)) = i;
                    ROW2(isnan(ROW2)) = ORDER(i);
                end
                ORDER = unique(ROW2(IND), 'stable');
                ORDER = unique([ORDER' 1:NR], 'stable');
            end
            RowOrder(1,:) = [];
            ROW{1}(ROW{1}~=0) = ROW2;

            clear ROW2 y
            %% Put the other samples in the same order
            if numel(smplnms)>1
                for i=2:numel(smplnms)
                    ROW2 = ROW{i}(ROW{i}~=0);
                    for j=1:size(RowOrder, 1)
                            ROW2(ROW2==RowOrder(j,1)) = nan;
                            ROW2(ROW2==RowOrder(j,2)) = RowOrder(j,1);
                            ROW2(isnan(ROW2)) = RowOrder(j,2);
                    end
                    ROW{i}(ROW{i}~=0) = ROW2;
                end
            end

            %% Add the ROW data to the ungated cell objects
            if ~isvalid(vPD)
                return;
            end
            waitbar(0.8, vPD, 'Adding region information to cell objects');
            for i=1:numel(smplnms)
                ttl{i} = smplnms{i};
                % Put the region data into the samples data table
                % If neighborhoods are cell centered
                if strcmp(DataType, 'Cell Centered Neighborhoods')
                    app.data.(smplnms{i}).MFICCN.(strcat(Constants.other_tag, ModelName)) = ROW{i};
                    % Pull the unique neighborhoods
                    [PositionsN, INDN] = unique(table2array(app.data.(smplnms{i}).MFICCN(:, {'X', 'Y', 'Z'})), 'rows', 'stable');
                    % Pull out the classification of the unique neighborhoods
                    Class =  ROW{i};
                    Class = Class(INDN);

                    % Pull the positions of all cells
                    PositionsCells = table2array(app.data.(smplnms{i}).AllCells(:, {'X', 'Y', 'Z'}));

                    % Break the data into smaller chunks so I
                    % don't break my computer
                    SzLim = 5000;
                    Nloop = ceil(size(PositionsCells, 1)/SzLim);
                    for loopi = 1:Nloop
                        st =((loopi-1)*SzLim+1);
                        en = min((st+SzLim-1), size(PositionsCells, 1));
                        %%% Do the actual distance calculation
                        Dist = sqrt( (PositionsCells(st:en,1)'-PositionsN(:,1)).^2 ...
                                +(PositionsCells(st:en,2)'-PositionsN(:,2)).^2 ...
                                +(PositionsCells(st:en,3)'-PositionsN(:,3)).^2 );
                        % Finds the index for closest X-Y-Z coordinate in PositionsN for each cell
                        % INDPosN is the index, in PositionN for the Nth Cell in PositionsCells
                        [~, INDPosN] = min(Dist);
                        %%% add the classification as a column in the data
                        app.data.(smplnms{i}).AllCells.(strcat(Constants.other_tag, ModelName))(st:en) = Class(INDPosN);
                    end % end of loopi matrix resizing
                    clear Dist
                % If neighborhoods are raster scanned
                elseif strcmp(DataType, 'Raster Scanned Neighborhoods')
                    app.data.(smplnms{i}).MFIRSN.(strcat(Constants.other_tag, ModelName)) = ROW{i};
                    % Put the region info in the AllCells matrix
                    % The ROW{i} vector should have the same x-y coordinates as
                    % the CellularityM matrix
                    % Find the unique coordinates for the Neighborhoods
                    x = unique(app.data.(smplnms{i}).MFIRSN.X);
                    y = unique(app.data.(smplnms{i}).MFIRSN.Y);
                    z = unique(app.data.(smplnms{i}).MFIRSN.Z);
                    [~, XD] = min(abs(app.data.(smplnms{i}).AllCells.X-x'), [],2);
                    [~, YD] = min(abs(app.data.(smplnms{i}).AllCells.Y-y'), [],2);
                    [~, ZD] = min(abs(app.data.(smplnms{i}).AllCells.Z-z'), [],2);
                    % 2D
                    if sum(strcmp(app.data.(smplnms{i}).MFIRSN.Properties.VariableNames, 'Z'))==0
                        NXSteps = numel(app.data.(smplnms{i}).MFIRSN.X)/numel(unique(app.data.(smplnms{i}).MFIRSN.X));
                        Region = ROW(NXSteps.*(XD-1)+YD);
                        app.data.(smplnms{i}).AllCells.(strcat(Constants.other_tag, ModelName)) = Region;
                    end
                    % Add the ROW data to the ungated cell objects 3D
                    if sum(strcmp(app.data.(smplnms{i}).MFIRSN.Properties.VariableNames, 'Z'))~=0
                        % Position of the nearest neighborhood for each cell
                        Position = [x(XD), y(YD), z(ZD)];
                        % Index in CellularityM of the closest neighborhood for each cell
                        [~, INDReg] = ismember(Position, table2array(app.data.(smplnms{i}).MFIRSN(:, {'X', 'Y', 'Z'})), 'rows');
                        Region = ROW{i};
                        Region = Region(INDReg);
                        app.data.(smplnms{i}).AllCells.(strcat(Constants.other_tag, ModelName)) = Region;
                    end
                elseif strcmp(DataType, 'Individual Cells')
                    app.data.(smplnms{i}).AllCells.(strcat(Constants.other_tag, ModelName)) = ROW{i};
                end % For RSN or CCN
            end % For All Samples
            %%%%%

            % Put some information about the neural network into the
            % userdata field, only if it was trained
            if new_model
                userdata = struct;
                userdata.FileNames = app.data.(smplnms{1}).FNames;
                userdata.SampleNames = smplnms;
                userdata.RowOrder = RowOrder;
                userdata.MFI = MFIList;
                userdata.MFI_weights = MFI_weights;
                userdata.pheno = phenotypes;
                userdata.pheno_weights = pheno_weights;
                userdata.SortNames = SortNames;
                userdata.DataPrep = DataPrep;
                userdata.DataType = DataType;
                userdata.NRegCValues = NRegCValues;
                userdata.NofClassesAlg = NofClassesAlg;
                userdata.NormOpt = NormOpt; 

                app.net.(ModelName).userdata = userdata;
            end


            %% Plot network sorted image
            if ~isvalid(vPD)
                return;
            end
            waitbar(0.95, vPD, 'Plotting');
            for i=1:numel(smplnms)
                Class = ROW{i};
                % 3D Data
                if any(strcmp(Dat_PreALL.Properties.VariableNames, 'Z'))
                    % If the neighborhoods are a raster scan
                    if strcmp(DataType, 'Raster Scanned Neighborhoods')
                        Class = reshape(Class, [numel(unique(app.data.(smplnms{i}).MFIRSN.Y)), ...
                             numel(unique(app.data.(smplnms{i}).MFIRSN.X)), numel(unique(app.data.(smplnms{i}).MFIRSN.Z))]);
                        Class = median(Class, 3);

                        %% Build Surfaces around each region
                        % options = {Radius, Region Threshold, Hole Threshold, Color, Alpha, SurfaceTyp, SurfaceName, Phenotype, SampleName}
                        % Empty NN-Surfaces Data
                        app.data.(smplnms{i}).Surfaces.RSN = struct;
                        for j=unique(ROW{i})'
                            options = {2.*(app.rwindowRSN), 50, 50, 'r', 0.5, 'RSN', ['RSN_' num2str(j)],'RSNeighborhoods', smplnms{i}};
                            alpha = options{1};
                            region_threshold = options{2};
                            hole_threshold = options{3};
                            type = options{6};
                            name = options{7};
                            phenotypes = options{8};
                            sample = options{9};
                            dat = app.data.(smplnms{i}).MFIRSN(ROW{i}==j, {'X', 'Y', 'Z'});
                            Plt_Helper.surf_create(app, sample, phenotypes, dat, 'X', 'Y', 'Z', alpha, type, name, region_threshold, hole_threshold,true);
                        end
                    elseif strcmp(DataType, 'Cell Centered Neighborhoods')
                        %% Build Surfaces around each region
                        Position = app.data.(smplnms{i}).MFICCN(:, {'X', 'Y', 'Z'});

                        % Empty NN-Surfaces Data
                        app.data.(smplnms{i}).Surfaces.CCN = struct;
                        for j=unique(ROW{i})'
                            options = {2.*(app.rwindowCCN), 50, 50, 'r', 0.5, 'CCN', ['CCN_' num2str(j)],'CCNeighborhoods', smplnms{i}};
                            alpha = options{1};
                            region_threshold = options{2};
                            hole_threshold = options{3};
                            type = options{6};
                            name = options{7};
                            phenotypes = options{8};
                            sample = options{9};
                            dat = Position(ROW{i}==j, :);
                            Plt_Helper.surf_create(app, sample, phenotypes, dat, 'X', 'Y', 'Z', alpha, type, name, region_threshold, hole_threshold, true);
                        end
                        Position = table2array(Position);
                    end
                end

                if strcmp(DataType, 'Raster Scanned Neighborhoods') && if_plt
                    Plotting.func_newfig(app,'P', {'Density/MFI RSN'}, 'S', smplnms{i}, 'C', Helper.full_var(ModelName));
                    fig = gcf;
                    cla
                    fig.Color = app.GUIOPTS.bgclr;
                    fig.InvertHardcopy = 'off';

                    imagesc('XData', unique(app.data.(smplnms{i}).MFIRSN.X), 'YData', unique(app.data.(smplnms{i}).MFIRSN.Y), 'CData', Class, [0 round(max(max(Class)))])

                    axis tight
                    axis equal
                    ax = gca;
                    ax.Color = app.GUIOPTS.bgclr;
                    ax.YColor = app.GUIOPTS.txtclr;
                    ax.XColor = app.GUIOPTS.txtclr;
                    ax.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
                    xlabel('x, \mum')
                    ylabel('y, \mum')
                    ax.Title.String = smplnms{i};
                    % Define the Colorbar
                    colormap(app.net.(ModelName).cmap);
                    caxis([0 app.net.(ModelName).NR])
                    c = colorbar;
                    % c.Color = app.GUIOPTS.txtclr;
                    c.Ticks = 0:1:app.net.(ModelName).NR;
                    c.Limits = [0, app.net.(ModelName).NR];
                    ylabel(c, 'Region number');
                elseif strcmp(DataType, 'Cell Centered Neighborhoods') && if_plt

                    Plotting.func_newfig(app, 'P', {'Density/MFI CCN'}, 'S', smplnms{i},'C', Helper.full_var(ModelName));
                    fig = gcf;
                    cla
                    fig.Color = app.GUIOPTS.bgclr;
                    fig.InvertHardcopy = 'off';

                    plot3(app.data.(smplnms{i}).AllCells.X,app.data.(smplnms{i}).AllCells.Y,app.data.(smplnms{i}).AllCells.Z, ...
                        '.', 'Color', [0.9 0.9 0.9],'MarkerSize',3)
                    grid on
                    hold on
                    plot3(Position(:,1),Position(:,2),Position(:,3),'or','MarkerSize',3)
% % %                             title(smplnms{i}, 'Color', app.GUIOPTS.txtclr)
                    axis tight
                    axis equal
                    ax = gca;
                    ax.Color = app.GUIOPTS.bgclr;
                    ax.YColor = app.GUIOPTS.txtclr;
                    ax.XColor = app.GUIOPTS.txtclr;
                    ax.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
                    xlabel('x, \mum')
                    ylabel('y, \mum')
                    zlabel('z, \mum')
                    ax.Title.String = smplnms{i};
                    % Change the view to be face on
                    view(0,90)
%                                 trisurf(Bndry,Position(Logical,1),Position(Logical,2),Position(Logical,3),'Facecolor','b','FaceAlpha',0.1)
                    [Xs,Ys,Zs] = sphere;
                    for j = 1:numel(Position(:, 1))
%                                     surf(Xs*app.rwindowCCN+Position(j,1), Ys*app.rwindowCCN+Position(j,2), Zs*app.rwindowCCN+Position(j,3), 'Facecolor',app.map(Class(j)+1, :), 'EdgeColor', 'none','FaceAlpha',0.3);
                        surf(Xs*app.rwindowCCN+Position(j,1), Ys*app.rwindowCCN+Position(j,2), Zs*app.rwindowCCN+Position(j,3), 'Facecolor',app.net.(ModelName).cmap(Class(j)+1, :), 'EdgeColor', 'none','FaceAlpha',0.3);
                    end
                    colormap(app.net.(ModelName).cmap);
                    c = colorbar;
                    c.Color = app.GUIOPTS.txtclr;
                    caxis([0 app.net.(ModelName).NR])
                    c.Ticks = 0:1:app.net.(ModelName).NR;
                    c.Limits = [0, app.net.(ModelName).NR];
                    ylabel(c, 'Region number')
                elseif strcmp(DataType, 'Individual Cells') && if_plt
                    if isempty(SortNames)
                        Plotting.func_newfig(app,'P', 'All Cells', 'S', smplnms{i}, 'C', Helper.full_var(ModelName));
                        ax = gca;
                        ax.Title.String = smplnms{i};
                    else
                        Plotting.func_newfig(app,'P', SortNames, 'S', smplnms{i}, 'C', Helper.full_var(ModelName));
                        ax = gca;
                        ax.Title.String = smplnms{i};
                    end
                end
            end % end of plots

            Helper.func_closeVPD(app, vPD);
%%%%%
userdata = app.net.(ModelName).userdata;
'Elapsed time to cluster'
b = toc(t0)
userdata.RunTime = b;
app.net.(ModelName).userdata = userdata;
%%%%%
            function cancel_waitbar_callback(hObject)
                delete(ancestor(hObject, 'figure'));
            end

        end

        function Save_NN(app)
            % SAVE_NN Front-End of the Save Model function
            %
            % Input:
            %   - app - Instance of CytoMap
            %
            % Creates:
            %   - File containing models that were choosen to be saved
            
            [file,path] = uiputfile({'*.mat'},'Save file name');
            list = Helper.full_var(fieldnames(app.net));
            [indx, ~] = listdlg('ListString',list, ...
                'Name', 'Select Models to save');

            if isempty(indx)
                return;
            end
            models = fieldnames(app.net);
            models = models(indx);
            Analysis.Save_NN_backend(app, file, path, models);
        end
        
        function Save_NN_backend(app, file, path, models)
            % SAVE_NN_BACKEND Back-End of the Save Model function
            %
            % Input:
            %   - app - Instance of CytoMap
            %   - file - name of file that models should be saved to
            %   - path - path to that directory which should contain that
            %       file
            %   - models - names of models which were chosen to be saved
            %       (in form which allows to read from app.net)
            %
            % Creates:
            %   - File containing models that were choosen to be saved
            
            curr_path = pwd;
            cd(path);

            dat = struct;
            
            for m = 1:numel(models)
                dat.(models{m}) = app.net.(models{m});
            end

            save(file,'dat');

            cd(curr_path);
        end

        function Load_NN(app, NetName)
            % LOAD_NN Front-End to Load Neural Network function
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - fnm - name of file that models should be loaded from
            %   - pnm - path to that directory which contains that file
            %
            % Modifies:
            %   - app - Adds a loaded in net under app.net
            
            [fnm , pnm]=uigetfile({'*.mat', '*.mat'}, ...
                      'Load Trained Neural Network SOM','Multiselect', 'off');
            Analysis.Load_NN_backend(app, pnm, fnm, NetName);
        end

        function Load_NN_backend(app, pnm, fnm, NetName)
            % LOAD_NN_BACKEND Back-End to Load Neural Network function
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - pnm - 
            %
            % Modifies:
            %   - app - Adds a loaded in net under app.net
            
            curr_path = pwd;
            cd(pnm);
            load(fnm, 'dat')

            if sum(ismember({'net', 'NReg'}, fieldnames(dat))) == 2  % Old version of net (if someone uses it)
                name = strcat('LoadedModel', num2str(numel(app.net) + 1));
                app.net.(name) = struct;
                app.net.(name).name = name;
                app.net.(name).type = 'SOM';
                app.net.(name).Network = dat.net;
                app.net.(name).NR = dat.NReg;
                app.net.(name).cmap = colormap(jet(dat.NReg));
                app.net.(name).userdata = dat.userdata;
            else
                models = fieldnames(dat);
                for m=1:numel(models)
                    if ~Helper.add_model(app, dat, models{m})
                        break;
                    end
                end
            end
            cd(curr_path);
            
            if ~isempty(NetName)
                NetName.Items = [{Constants.new_model}, Helper.full_var(fieldnames(app.net))'];
            end

        end
    end
 end

 