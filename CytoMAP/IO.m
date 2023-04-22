classdef IO
    %IO Summary of this class goes here
    %   Detailed explanation goes here

    methods (Static)
        function Import_Definitions_Func
% % %             import tree.*;
% % %             import Constants.*;
% % %             import Helper.*;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main load functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function func_load(app, NSamp, type)
            % FUNC_LOAD it is front end of loading in .wsp and .mat files 
            %
            % INPUT:
            %   - app - Instance of CytoMAP
            %   - NSamp - int - either 0 or 1. If 0 mutliple files from one
            %       folder will be loaded. If files from multiple folders
            %       will be loaded.
            %   -type the expected file types, .csv, xlsx, or .mat
            %
            % MODIFIES:
            %   - app - If csv files are loaded then only app.data is
            %       modified. If .mat is loaded then almost any field in
            %       app can be modified.
            
            current_folder = pwd;

            % Select multiple files
            if NSamp==0
                %Input the Folder paths and File paths
                switch type
                    case 'type1'
                    [fnm , pnm]=uigetfile({...
                        '*.csv;*.xlsx;*xls', 'Compatible Files';...
                        '*.csv', 'Comma Seperated Table of Cells (*.csv)';...
                        '*.xls', 'Table of Cells (*.xls)';...
                        '*.xlsx', 'Table of Cells (*.xlsx)'...
                        }, ...
                        'Select Files to load, individual files are treated as individual cell types','Multiselect', 'on');
                    case 'type2'
                    [fnm , pnm]=uigetfile({...
                        '*.mat;', 'Compatible Files';...
                        '*.mat', 'CytoMAP Workspace Files (*.mat)'...
                        }, ...
                        'Select CytoMAP Worskapce file (.mat) you want to load','Multiselect', 'on');
                        
                end
                %Convert character array to cell array
                if ischar(fnm)
                    fnm = {fnm};
                end
                % If no file was selected cancel function
                if ~ischar(pnm)
                    return
                end

                pnm = {pnm};
                fnm = {fnm};
                cd(pnm{1});
            end

            % Select multiple folders
            if NSamp==1
                %Input the Folder paths and File paths
                pnm = IO.uigetdir2('', ...
                    'Shift + Select folders with .csv tables of cells in them');
                % If no file was selected cancel function
                if isempty(pnm)
                    return
                end
                
                fnm = cell(numel(pnm), 1);
                for pnm_i = 1:numel(pnm)
                    % If the user selected files instead of folders
                    if endsWith(pnm{pnm_i}, '.csv')
                        pnm_split = strsplit(pnm{pnm_i}, filesep);
                        pnm{pnm_i} = strjoin(pnm_split(1:(end-1)), filesep);
                        fnm_i = pnm_split{end};
                        fnm{pnm_i} = strcat(filesep,{fnm_i});
                    else
                        % If the user selected folders with .csv fils in them
                        fnm_i = dir(fullfile(pnm{pnm_i}, '*.csv'));
                        % Pull out just the names
                        fnm{pnm_i} = strcat(filesep,{fnm_i.name});
                    end
                    % Access filenames using fnm{pnm_i}{fnm_i} syntax
                end
            end
            IO.func_load_backend(app, fnm, pnm, current_folder);
        end

        function func_load_backend(app, fnm, pnm, current_folder)
            % FUNC_LOAD_BACKEND it is a back end of loading in .wsp and .mat
            % files.
            % Basically wrap around func_load_mat and func_load_csv.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - fnm - cell - Filesnames which to load into CytoMAP.
            %   - pnm - cell - Folders in which those files are
            %   - current_folder - folder to which function should return.
            %
            % Modifies:
            %   - app - Loads up everything from given filenames.
            
            if strcmp(fnm{1}{1}((end-3):end), '.mat')
                IO.func_load_mat(app, fnm, current_folder);
% % %             elseif strcmp(fnm{1}{1}((end-3):end), '.csv')
            else
                IO.func_load_csv(app, fnm, pnm, current_folder);
            end
        end

        function func_load_mat(app, fnm, current_folder)
            % FUNC_LOAD_MAT Loads data from .mat into current instance of
            % CytoMAP
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - fnm - cell - Filesnames which to load into CytoMAP.
            %   - current_folder - folder to which function should return.
            %
            % Modifies:
            %   - app - Loads up everything from given filenames.
            %
            % Note: THIS CURRENTLY OVERRIDES ANY EXISTING DATA, except for
            % app.data

            vPD = waitbar(0, 'Initializing');
            waitbar(0.1, vPD, 'Loading data');

            if ~all(endsWith(fnm{1}, '.mat'))
                errordlg('All files have to be the same type');
                return;
            end

            for f_i = 1:numel(fnm{1})
                dat = table;
                load(fnm{1}{f_i});

                % Merge structures
                new_fields = fieldnames(dat.data);
                for field_n = 1:numel(new_fields)
                    app.data.(new_fields{field_n}) = dat.data.(new_fields{field_n});
                end

                if ~isempty(fieldnames(dat.net))
                    mat_net2net(app, dat.net);
                end

                if sum(contains(fieldnames(dat), 'CellInfo'))~=0
                    app.CellInfo = dat.CellInfo;
                end

                if ~isstruct(dat.map)
                    dat.map = struct;
                    dat.map.jet = jet;
                    % Put in cases for RSN and CCN
                    app.map = dat.map;
                else
                    app.map = dat.map;
                end

                app.points = dat.points;
                app.polygons = dat.polygons;

                % Support for new and old versions of DataN Items
                if isfield(dat, 'DataN')
                    new_items = dat.DataN.Items;
                else
                    new_items = dat.DataNItems;
                end

                % Concatanate already existing Items, and new ones. Without duplicates.
                for field_n = 1:numel(new_items)
                    if ~ismember(new_items{field_n}, app.DataN.Items)
                        app.DataN.Items{end + 1} = new_items{field_n};
                    end
                end

                if length(app.DataN.Items) > 1 && strcmp(app.DataN.Items{1}, 'No Data Loaded')
                    app.DataN.Items = app.DataN.Items(2:end);
                end

                app.DataN.Value = app.DataN.Items{1};

                % Convert to Tree Structure
                % Set scan type, with bias to Raster Scan Neighborhoods
                if isfield(app.data.(app.DataN.Value), 'MFIRSN')
                    app.data.(app.DataN.Value).ScanType = 'MFIRSN';
                elseif isfield(app.data.(app.DataN.Value), 'MFICCN')
                    app.data.(app.DataN.Value).ScanType = 'MFICCN';
                end

                if isfield(dat, 'NewGate')
                    % "Folowing Gates were Found"
                    dat.NewGate
                    if ~isfield(app, 'NewGate')
                        app.NewGate = dat.NewGate;
                        % "NewGate is found, no existing gates"
                    else
                        overlap = ~ismember(fieldnames(dat.NewGate), fieldnames(app.NewGate));
                        overlap = overlap & ~strcmp(fieldnames(dat.NewGate), 'internal__curr');
                        new_fields = fieldnames(dat.NewGate);
                        % "Gates were found New Gates to Load = "
                        new_fields = new_fields(overlap);
                        for field_idx = 1:numel(new_fields)
                            app.NewGate.(new_fields{field_idx}) = dat.NewGate.(new_fields{field_idx});
                        end
                    end
                else
                    % "No Gates were found"
                end
                
                % Convert to the tree structure
                % THIS DOES NOT CONVERT MFIRSN OR MFICCN VARIABLE NAMES
                for smpl_idx = 1:numel(app.DataN.Items)
                    smpl_i = app.DataN.Items{smpl_idx};
                    app.DataN.Value = smpl_i;
                    app.data.(smpl_i) = mat2tree(app.data.(smpl_i));
                end
            end

            for smpl_idx = 1:numel(fieldnames(app.data))
                smpl_i = fieldnames(app.data);
                smpl_i = smpl_i{smpl_idx};
                if isfield(app.data.(smpl_i), 'Gates')
                    if ~ismember('NewGate', fields(app)) || isempty(fieldnames(app.NewGate))
                        app.data.(smpl_i) = rmfield(app.data.(smpl_i), 'Gates');
                    else
                        % Be warned loading gates does not have working
                        % error cheking
%                         gates_i = app.data.(smpl_i).Gates;
%                         app.data.(smpl_i).Gates = gates_i(ismember(append('Gate_', gates_i), fieldnames(app.NewGate)));
                    end
                end
            end
            % Select the sample you just uploaded
            waitbar(1, vPD, 'Done!');
            close(vPD);
            cd(current_folder);

            function treed_struct = mat2tree(mat_struct)
                % TREED_STRUCT Takes a app.data.X struct from just loaded
                % .mat file, and converts it to style that .wsp files are
                % loaded in.
                % Additionally checks for version, so that if .mat file was
                % already new and in style of .wsp, then no redundant jobs
                % are done.
                %
                % Input:
                %   - mat_struct - app.data.X struct from just loaded .mat
                %       file.
                %
                % Output:
                %   - treed_struct - app.data.X struct which can be put
                %       into a current verison of CytoMAP.
                %
                % If it's not legacy version, then use it's data format.
                % This assumes all versions after legacy have structure similar to the .wsp/.mat files.
                
                if isfield(mat_struct, 'MetaData') && isfield(mat_struct.MetaData, 'Ver') && ~strcmp(mat_struct.MetaData.Ver, Constants.LEGACY_VER)
                    treed_struct = mat_struct;
                else
                    % If it's old, then create new struct to keep newly-formated data.
                    treed_struct = struct;
                    treed_struct.FNames = mat_struct.FNames;
                    treed_struct.Path = mat_struct.Path;
                    treed_struct.AllCells = mat_struct.AllCells;

                    if isfield(mat_struct, 'MetaData')
                        treed_struct.MetaData = mat_struct.MetaData;
                    else
                        treed_struct.MetaData = struct;
                        treed_struct.MetaData.Group = 1;
                        treed_struct.MetaData.Annotation = 0;
                    end
                    treed_struct.MetaData.Ver = Constants.CURR_VER;
                end

                % Check for non-basic fields. If they do not exist, either create stubs, or ignore.
                if isfield(mat_struct, 'Surfaces')
                    treed_struct.Surfaces = mat_struct.Surfaces;
                else
                    treed_struct.Surfaces = struct;
                    treed_struct.Surfaces.UDS = struct;
                    treed_struct.Surfaces.CCN = struct;
                    treed_struct.Surfaces.RSN = struct;
                end

                if isfield(mat_struct, 'StatsTBL')
                    treed_struct.StatsTBL = mat_struct.StatsTBL;
                end

                if isfield(mat_struct, 'ScanType')
                    treed_struct.ScanType = mat_struct.ScanType;
                end

                % If tree and GateTags already exist, there is no need to repeat the process.
                if isfield(treed_struct, 'tree') && isfield(treed_struct, 'GateTags')
                    return
                end

                treed_struct.GateTags = table;

                % For older data the SortedDat column names need to be
                % changed to Gate_1 Gate_2 etc.

                gate_names = mat_struct.SortedDat.Properties.VariableNames(startsWith(mat_struct.SortedDat.Properties.VariableNames, Constants.gate_tag));
                if isempty(gate_names)  % Backwards compatibility to .mat which did not have the gate tag.
                    gate_names = mat_struct.SortedDat.Properties.VariableNames(~ismember(mat_struct.SortedDat.Properties.VariableNames, {'X', 'Y', 'Z'}));
                end
                gate_tags = cell(size(gate_names));
                main_tree = tree('All', 'tag', [Constants.gate_tag '0']);
                tag = 1;

                for gate_idx = 1:numel(gate_names)
                    gate_kid = gate_names(gate_idx);
                    kid = tree(gate_kid{1, 1}, 'gate_type', 'logic', 'gate_points', mat_struct.SortedDat(:, gate_kid), 'tag', sprintf([Constants.gate_tag '%d'], tag));
                    main_tree = main_tree.add_kid(kid);
                    % Add the full name
                    FullName = strcat(main_tree.name, '/', kid.name);
                    % Add only the immediate parent and kid name
                    ShortName = strsplit(main_tree.name, '/');
                    ShortName = [ShortName{end} '/' kid.name];
                    ShortName = Helper.full_gate(ShortName);
                    treed_struct.GateTags.(sprintf([Constants.gate_tag '%d'], tag)) = {FullName; ShortName};
                    gate_tags(gate_idx) = {sprintf([Constants.gate_tag '%d'], tag)};
                    tag = tag + 1;
                end

                % Load MFI scans. If ScanType exists load it. If not priority is given to MFIRSN.
                if isfield(mat_struct, 'MFICCN')
                    treed_struct.MFICCN = mat_struct.MFICCN;
                    treed_struct.ScanType = 'MFICCN';
                end

                if isfield(mat_struct, 'MFIRSN')
                    treed_struct.MFIRSN = mat_struct.MFIRSN;
                    treed_struct.ScanType = 'MFIRSN';
                end
                
                if isfield(mat_struct, 'Gates')
                    treed_struct.Gates = mat_struct.Gates;
                end

                treed_struct.AllCells = Helper.reorder_cols(treed_struct.AllCells);
                treed_struct.AllCells.Properties.VariableNames = Helper.valid_channel(treed_struct.AllCells.Properties.VariableNames);
                treed_struct.AllCells.Properties.VariableNames = Helper.valid_other(treed_struct.AllCells.Properties.VariableNames);
                [to_drop, new_var_names] = Helper.remove_other_duplicates(treed_struct.AllCells.Properties.VariableNames);

                treed_struct.AllCells = treed_struct.AllCells(:, ~to_drop);

                treed_struct.AllCells.Properties.VariableNames = new_var_names;
                treed_struct.AllCells = main_tree.gate_cells(treed_struct.AllCells);
                treed_struct.tree = main_tree;


                % TODO: When importing older versions of .mat, make sure that MFI... column names are also transformed
                if isfield(treed_struct, 'MFIRSN')
                    for gate_idx = 1:numel(gate_names)
                        treed_struct.MFIRSN.Properties.VariableNames(ismember(treed_struct.MFIRSN.Properties.VariableNames, gate_names(gate_idx))) = gate_tags(gate_idx);
                        treed_struct.MFIRSN.Properties.VariableNames = Helper.valid_channel(treed_struct.MFIRSN.Properties.VariableNames);
                    end
                end
                if isfield(treed_struct, 'MFICCN')
                    for gate_idx = 1:numel(gate_names)
                        treed_struct.MFICCN.Properties.VariableNames(ismember(treed_struct.MFICCN.Properties.VariableNames, gate_names(gate_idx))) = gate_tags(gate_idx);
                        treed_struct.MFICCN.Properties.VariableNames = Helper.valid_channel(treed_struct.MFICCN.Properties.VariableNames);
                    end
                end
            end

            function mat_net2net(app, mat_net_struct)
                % MAT_NET2NET Converts an old nets into the new version of
                % net in current CytoMAP
                %
                % Input:
                %   - app - Instance of CytoMAP
                %   - mat_net_struct - struct corresponding to app.net from
                %       loaded .mat file.
                %
                % Modifies:
                %   - app - Specifically app.net which will be populated
                %       with nets from the given struct.
                
                if sum(ismember({'net', 'NReg'}, fieldnames(mat_net_struct))) == 2  % Old version of net (if someone uses it)
                    name = strcat('LoadedModelFromMat', num2str(numel(app.net) + 1));
                    app.net.(name) = struct;
                    app.net.(name).name = name;
                    app.net.(name).type = 'SOM';
                    app.net.(name).Network = mat_net_struct.net;
                    app.net.(name).NR = mat_net_struct.NReg;
                    app.net.(name).cmap = colormap(jet(mat_net_struct.NReg));
                else
                    fnames = fieldnames(mat_net_struct);
                    for net = 1:numel(fnames)
                        if ~Helper.add_model(app, mat_net_struct, fnames{net})
                            return;
                        end
                    end
                end
            end
        end

        function func_load_csv(app, fnm, pnm, current_folder, web)
            % FUNC_LOAD_CSV Loads data from .csv files into current
            % instance of CytoMAP
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - fnm - cell - Filesnames which to load into CytoMAP.
            %   - pnm - cell - Folders in which those files are
            %   - current_folder - folder to which function should return.
            %
            % Modifies:
            %   - app - Specifically app.data, in which samples are put in.
            if nargin<6
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            vPD = waitbar(0, 'Initializing');
            waitbar(0.1, vPD, 'Loading data');
                       
% % %             if ~all(endsWith(fnm{1}, '.csv')) || ~all(endsWith(fnm{1}, '.xls')) || ~all(endsWith(fnm{1}, '.xlsx'))
% % %                 errordlg('All files have to be the same type');
% % %                 return;
% % %             end

            SampleNames = IO.get_sample_names(app, pnm, fnm);
            try
                dat = IO.get_dat(pnm, fnm, SampleNames);
            catch
                disp(['ERROR in path: ' pnm ...
                    '  File: ' fnm ...
                    ' Sample: ' SampleNames])
                waitbar(0, vPD, 'ERROR loading data, check column headers!');
                return
            end

            close(vPD);

            % Make a table to check the cell names etc.
            if isempty(SampleNames)
                % Break or do some stuff
                error('SampleNames empty');
            end

            % Build the options menus
            UIfig = uifigure('Name', 'File Import Options');
            if web==1
                UIfig.Visible='OFF';
            end
            UIfig.Position = alpha*[10 10 1000 800];
            % Create the table of options
            t = uitable(UIfig);
            t.Data = IO.get_table_data(fnm, dat, SampleNames);
            t.Position = alpha*[0 70 1000 730];
            t.ColumnName = { ...
                'Sample Groups', 'Sample Annotation','Sample Names', ...
                'File Names', 'Import Cells', 'Cell Names','String to Replace', 'Replace with', 'Final Cell Names', ...
                'Channel Names', 'ChNReplace', 'ChNReplacements', 'Final Channel Names' ...
            };
            t.ColumnEditable = [true true true false true false true true, true, false, true, true, true];
            t.ColumnWidth = {alpha*50, alpha*100, alpha*350, alpha*50, alpha*100, alpha*125, alpha*100, alpha*125};

            lbl = uilabel(UIfig); lbl.Text = sprintf('Replace Missing Values With:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = alpha*[10+100+100+30 32 200 15];
            swch = uiswitch(UIfig);
            swch.Items = {'NaN','Value:'};
            swch.Position = [10+100+100+60 10 50 30];
            
            numrep = uieditfield(UIfig,'numeric');
            numrep.Value = 0;
            numrep.Position = alpha*[10+100+100+60+90 10 50 22];
            
            
            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) func_load_rename(t));
            btn.Position = alpha*[10, 5, 100, 50];
            btn.Text = 'ReName';

            % Create a push button
            btn = uibutton(UIfig,'push', ...
                'ButtonPushedFcn', @(btn,event) IO.func_load_execute(UIfig, app, dat, fnm{end}, pnm, SampleNames{end}, t.Data, swch, numrep.Value));
            btn.Position = alpha*[10+100, 5, 100, 50];
            btn.Text = 'Load';

            cd(current_folder);

            % rename channels and cells
            function func_load_rename(t)
                % FUNC_LOAD_RENAME Renames both channels and cells from a
                % given table, based on how user wants to change those
                % names.
                %
                % Input:
                %   - t - table - table, in which columns 9 and 13 are
                %       names to be renamed, and where columns 7 and 11 are
                %       characters that should be replaced, and columns 8
                %       and 12 are characters which those characters will
                %       be replaced with
                %
                % Modifies:
                %   - Modifies that given table on columns 9 and 13 (which
                %   are channels and cells)
                
                fnms = t.Data(:, 9);
                IND = find(~cellfun(@isempty,fnms));
                for i=1:numel(t.Data(:, 7))
                    if ~isempty(t.Data{i, 7})
                        if ~strcmp(t.Data{i, 8}, 'Characters to replace it with')
                            fnms(IND) = strrep(fnms(IND), t.Data{i, 7}, t.Data{i, 8});
                        end
                    end
                end
                t.Data(:, 9) = fnms;

                fnms = t.Data(:, 13);
                IND = find(~cellfun(@isempty,fnms));
                for i=1:numel(t.Data(:, 11))
                    if ~isempty(t.Data{i, 11})
                        if ~strcmp(t.Data{i, 12}, 'Characters to replace it with')
                            fnms(IND) = strrep(fnms(IND), t.Data{i, 11}, t.Data{i, 12});
                        end
                    end
                end
                t.Data(:, 13) = fnms;
            end
        end

        function func_IMSload(app)
            % FUNC_IMSLOAD Work in progress
            %
            % Add the path to the IMaris reader stuff
            % addpath('ImarisReader-master');
            
            addpath('ImarisReader-master');
            % Select the file
            % returns the filename (fnm) and the path name (pnm)
            [fnmParent, pnm]=uigetfile({'*.ims'}, 'Select the .ims file','Multiselect', 'off');
            % Open the file
            smpl = app.DataN.Value;
            app.data.(smpl).IMSDat.imsObj = ImarisReader([pnm fnmParent]);
            app.data.(smpl).IMSDat.ChInfo = app.data.(app.DataN.Value).IMSDat.imsObj.DataSet.ChannelInfo;
            app.data.(smpl).IMSDat.ChName = {app.data.(smpl).IMSDat.ChInfo.Name};
% % %             %% Pull a single frame
% % %             cIdx = 0;
% % %             tIdx = 0;
% % %             vol = imsObj.DataSet.GetDataVolume(cIdx, tIdx);
% % %             % vol = sum(vol(:, :, 1:20), 3);
% % %             %% Crop
% % %             % vol = vol(1:500, 1:500, :);
% % %             %% Build the figure
% % %             fig = figure(1);
% % %             clf
% % %             fig.Color = 'w';
% % %             fig.InvertHardcopy = 'off';
% % %             fig.Position(3:4) = [725 450];
% % %
% % %             % ZAdj = uicontrol('Parent',fig,'Style','slider');
% % %             % ZAdj.Position = [81,54,419,23];
% % %             % ZAdj.value = zDim/2;
% % %             % ZAdj.min = 0;
% % %             % ZAdj.max = zDim;
% % %             % Get dimensions of image
% % %             xDim = size(vol, 2);
% % %             yDim = size(vol, 1);
% % %             zDim = size(vol, 3);
% % %
% % %             % put in RGB color
% % %             % % % rgbImage = ind2rgb(grayImage, colormap);
% % %             color = [0 1 1];
% % %             % frame = cat(3, color(1).*sum(vol(:, :, 1:20), 3), color(2).*sum(vol(:, :, 1:20), 3), color(3).*sum(vol(:, :, 1:20), 3));
% % %             % imshow(frame./max(max(max(frame))), [0, 1])
% % %             % imshow(frame, [0, 255])
% % %
% % %             num = 3;
% % %             frame = cat(3, color(1).*sum(vol(:, :, num), 3), color(2).*sum(vol(:, :, num), 3), color(3).*sum(vol(:, :, num), 3));
% % %             im = image(frame./max(max(max(frame))));
% % %             axis equal
% % %             axis tight
% % %             im.CLim = [0 1];
% % %             %% Add a table
% % %             tbldat = cell(max([numel(ChName), 10]), 6);
% % %             tbldat(1, 1) ={true};
% % %             tbldat(2:numel(ChName), 1) ={false};
% % %             tbldat(1:numel(ChName), 2) =ChName;
% % %             % tbldat(1:numel(ChName), 3) ={[0 1 1]};
% % %             % tbldat(1, 3) ={[0 1 1]};
% % %
% % %             tbldat(1:4, 4) = {xDim, yDim, zDim, 0};
% % %             tbldat(1:4, 5) = {xDim, yDim, zDim, 255};
% % %             tbldat(1:4, 6) = {'X', 'Y', 'Z', 'C'};
% % %
% % %             %Experimental
% % %             % % %     % Create a push button
% % %             % % %     tbldat(1:numel(ChName), 2)
% % %             % % %     btn = uibutton(f,'push', 'ButtonPushedFcn', @(btn,event) BtnPush(app,t,f));
% % %             % % %     btn.Position = [550/2-25, 0, 100, 50];
% % %             % % %     btn.Text = 'Ok';
% % %             %
% % %
% % %             t = uitable(fig);
% % %             t.Data = tbldat;
% % %             t.Position = [500 50 550 400];
% % %             t.ColumnName = {'','Channel','min', 'max', ''};
% % %             t.ColumnEditable = [true, true, true, true, false];
% % %             t.ColumnWidth = {50 100 50, 50, 50, 20};
% % %
% % %             %% Add a slidebar
% % %             SliderH = uicontrol('Parent',fig,'style','slider','position',[10 10 20 500],...
% % %                 'min', 1, 'max', zDim, 'value', zDim/2);
% % %
% % %             addlistener(SliderH, 'Value', 'PostSet', @(source, eventdata) callbackfn(SliderH, vol, color));
% % %
% % %             function callbackfn(SliderH, vol, color)
% % %                 num = round(SliderH.Value);
% % %                 frame = cat(3, color(1).*sum(vol(:, :, num), 3), color(2).*sum(vol(:, :, num), 3), color(3).*sum(vol(:, :, num), 3));
% % %                 im = image(1.5.*frame./max(max(max(frame))));
% % %                 axis equal
% % %                 axis tight
% % %             end
% % %             function pick_color
% % %                 color = uisetcolor(color);
% % %             end

        end

        function func_FCSLoad(app)
            % FUNC_FCSLOAD Work in progress
            %
            % this will eventually be able to pull data directly from an .fcs instead of a .csv file without going through a .wsp
            return
            path = 'C:\Users\calebst\Downloads\ResubmissionFullSession.fcs';
            smpl = 'Sample1';

            [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(path);
            chnames = {fcshdr.par.name}';
            dataset = array2table(fcsdat, ...
                        'VariableNames', chnames);

            dat_t = struct;

            dat_t.Surfaces = struct;
            dat_t.Surfaces.UDS = struct;
            dat_t.Surfaces.CCN = struct;
            dat_t.Surfaces.RSN = struct;

            dat_t.MetaData = struct;
            dat_t.MetaData.Ver = Constants.CURR_VER;
            dat_t.MetaData.Annotation = 0;


            dat_t.AllCells = dataset;
            dat_t.AllCells.gate1 = ones(size(dataset,1),1);
            dat_t.GateTags = table;
            dat_t.GateTags.gate1 = {'AllCells'; 'AllCells'};
            dat_t.Path = path;
            dat_t.NewGate = struct;

            app.data.(smpl) = dat_t;
            app.DataN.Items = fieldnames(app.data);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % .csv Loading Helpers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function tab_data = get_table_data(fnm, dat, SampleNames)
            % GET_TABLE_DATA Returns a table which can be used for csvs
            % reading. Basically a map from samples, data to table
            %
            % Input:
            %   - fnm - cell - Names of files data was loaded from.
            %   - dat - struct - Struct to temporary version of csv struct
            %       which corresponds to app.data in CytoMAP.
            %   - SampleNames - cell - Samples which should be loaded in as
            %       csvs.
            %
            % Output:
            %   - tab_data - cell - Data of Table for loading csvs,
            %       corresponding to given parameters.
            %
            % Get a cell names and channel names of the first sample
            
            CellNames = fieldnames(dat.(SampleNames{1}).GatesTMP);
            ChannelNames = dat.(SampleNames{1}).GatesTMP.(CellNames{1}).Properties.VariableNames;
            
            % Get all unique cell names
            for smpl_idx=1:numel(SampleNames)
                smpl_c_names = fieldnames(dat.(SampleNames{smpl_idx}).GatesTMP);
                % Redefine the cell names to be all of the unique cell
                % names within all of the samples
                CellNames =  unique([CellNames; smpl_c_names], 'stable');
            end
            
            % Get sample with max number of channel names
            for smpl_idx=1:numel(SampleNames)
                smpl_c_names = fieldnames(dat.(SampleNames{smpl_idx}).GatesTMP);
                for s_c_n=1:numel(smpl_c_names)
                    if numel(dat.(SampleNames{smpl_idx}).GatesTMP.(smpl_c_names{s_c_n}).Properties.VariableNames) ...
                            > numel(ChannelNames)
                        ChannelNames = ...
                            dat.(SampleNames{smpl_idx}).GatesTMP.(smpl_c_names{s_c_n}).Properties.VariableNames;
                    end
                end
            end
            
            % Populate TableData
            tab_data = cell(max([size(SampleNames, 1), numel(CellNames), numel(ChannelNames)]), 13);
            tab_data(1:size(SampleNames), 1)  = {1}; % SamepleGroups
            tab_data(1:size(SampleNames), 2)  = {0}; % SamepleAnnotation
            tab_data(1:size(SampleNames), 3)  = SampleNames;   %SamepleNames

            tab_data(1:numel(fnm{end}), 4)  = fnm{end}; % FileNames

            tab_data(1:numel(CellNames), 5)   = {true}; %ImportCells
            tab_data(1:numel(CellNames), 6)   = CellNames; %CellNames
            tab_data(1, 7)                    = {'Characters to be replaced'}; %CellNReplace
            tab_data(1, 8)                    = {'Characters to replace it with'}; %CellNReplacements
            tab_data(1:numel(CellNames), 9)   = CellNames; %NewCeNames

            tab_data(1:numel(ChannelNames), 10) = ChannelNames; %ChannelNames
            tab_data(1, 11) = {'Characters to be replaced'}; %ChNReplace
            tab_data(1, 12) = {'Characters to replace it with'}; %ChNReplacements
            tab_data(1:numel(ChannelNames), 13) = ChannelNames; %NewChNames
        end

        function SampleNames = get_sample_names(app, pnm, fnm)
            % GET_SAMPLE_NAMES Returns names of all samples which will be
            % created from given folders with csv files.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - pnm - cell - Paths to folders containing .csv files which
            %       will be loaded into CytoMAP.
            %
            % Output:
            %   - SampleNames - cell - Size equal to pnm. Name of sample
            %       folder will correspond to once it's loaded into
            %       CytoMAP.
            
            SampleNames = cell(numel(pnm), 1);
            for pnm_i = 1:numel(pnm)
                    % If there is only one file associated with this folder
                    if numel(fnm{pnm_i})==1
                        DName = strrep(fnm{pnm_i}{1}, '.csv', '');
                    else
                        % Get name of Sample for this dataset
                        DName = split(pnm{pnm_i}, filesep);
                        DName = DName{end};
                    end

                % If there isn't a similar string in the filenames call this data, data1
                if isempty(DName)
                    if strcmp(app.DataN.Value, 'No Data Loaded')
                        DName = 'Sample_1';
                    else
                        DName = ['Sample_' num2str((numel(app.DataN.Items)+1))];
                    end
                else
                    DName = strcat('Sample_', DName);
                end

                DName = strrep(DName, 'export_', '');
                DName = Helper.valid_var(DName);

                SampleNames{pnm_i} = DName;
            end
        end

        function dat = get_dat(pnm, fnm, SampleNames)
            % GET_DAT Returns an early stage app.data which contains all
            % samples to be loaded, but not yet fully formed.
            %
            % Input:
            %   - pnm - cell - Paths to folders containing .csv files.
            %   - fnm - cell - Names of files of .csv files within these
            %       folders.
            %   - SampleNames - cell - Names of samples corresponding to
            %       given folders.
            %
            % Output:
            %   - dat - struct - Intermediate app.data struct. It will
            %       contain all samples as fields, but with only GatesTMP,
            %       which later will be extracted for more efficient
            %       program manipulation.
            
            n=0;
            dat = struct;
            for pnm_idx = 1:numel(pnm)
                pnm_name = SampleNames{pnm_idx};
                % Determine a similar pattern at start of each filename for given sample
                % It is considered redundant and will be removed from each gate.
                if numel(fnm{pnm_idx}) == 1
                    % Don't rename if one .csv is loaded in.
                    same_start_str = 1;
                else
                    same_start_str = fnm{pnm_idx}{1};
                    while ~all(startsWith(fnm{pnm_idx}, same_start_str, 'IgnoreCase', true))
                        if isempty(same_start_str)
                            break;
                        else
                            same_start_str = same_start_str(1:end-1);
                        end
                    end
                    same_start_str = numel(same_start_str) + 1;  % MatLab has inclusive indexing.
                end
                if pnm_idx==1
                    redefineX = 0;
                    redefineY = 0;
                    redefineZ = 0;
                end
                    
                for j = 1:size(fnm{pnm_idx}, 2)
                    n=n+1;

                    %open the file
                    dat00=importdata([pnm{pnm_idx} fnm{pnm_idx}{j}]);
                    % if there are no cells in the csv do something 
                    if iscell(dat00)
                        dtemp = struct;
                        dtemp.textdata = strsplit(dat00{1}, ',');
                        dtemp.colheaders = strsplit(dat00{1}, ',');
% % %                         dtemp.data = zeros(1, numel(dtemp.textdata));
                        dtemp.data = zeros(size(dat00, 1)-1, numel(dtemp.textdata));
                        % Find the number of headers
                        for ind_i = 1: size(dat00, 1)-1
                            try
                                dtemp.data(ind_i, :) = strsplit(dat00{ind_i+1}, ',');
                            catch
                                splitdat = strsplit(dat00{ind_i+1}, ',');                                
                                for col_i = 1:size(numel(dtemp.textdata))
                                    try
                                        dtemp.data(ind_i, col_i) = splitdat(col_i);
                                    catch
                                        
                                    end
                                end
                            end
                            
                        end
                        dat00 = dtemp;
                    elseif size(dat00.textdata, 1) > 1
                        % the first column is text so stuff loads weird
                        dtemp = struct;
                        dtemp.textdata = dat00.textdata(1, :);
% % %                         dtemp.data = zeros(1, numel(dtemp.textdata));
                        dtemp.data = cell(size(dat00.textdata, 1)-1, numel(dtemp.textdata));
                        % find columns of text data
                        colsIND = sum(cellfun(@any,dat00.textdata(:, :)))==1;
                        % Find the number of headers
                        for ind_i = 1: size(dat00.textdata, 1)-1
% % %                             colsIND = cellfun(@isempty,dat00.textdata(ind_i+1, :))
                            dtemp.data(ind_i, ~colsIND) = dat00.textdata(ind_i+1, ~colsIND);
                            dtemp.data(ind_i, colsIND) = num2cell(dat00.data(ind_i, :));
                        end
                        dat00 = dtemp;
                    end
                                        
                    %Define the gate name
                    if numel(fnm{pnm_idx}) == 1
                        gate = 'All';
                    else
                        %take out any spaces in the gate name
                        gate = fnm{pnm_idx}{j}(same_start_str:end);
                        gate = split(gate, '.');
                        gate = gate{1};
                        if isempty(gate)
                            gate = 'All';                        
                        end
                    end
                    gate = Helper.valid_gate(gate);
                    if iscell(gate)
                        gate = gate{1};
                    end

                    %Clean up the variable names
                    if isfield(dat00, 'colheaders')
                        VNames=strrep(dat00.colheaders, '-','');
                    else
                        VNames=strrep(dat00.textdata, '-','');
                    end

                    VNames = Helper.valid_other(VNames);
                    VNames = Helper.valid_channel(VNames);
                    
                    % if it is not the first sample loaded in make sure
                    % channel index didn't move around
                    if pnm_idx~=1 
                        if redefineX ~= 0 && ~ismember('X', VNames) && ismember(redefinedX, VNames)
                            % Do nothing
                        else % For all other cases ask to redefine 
                            redefineX = 0;
                        end
                        if redefineY ~= 0 && ~ismember('Y', VNames) && ismember(redefinedY, VNames)
                            % Do nothing
                        else % For all other cases ask to redefine
                            redefineY = 0;   
                        end
                        if redefineZ == 1 && ~ismember('Z', VNames) && ismember(redefinedZ, VNames)
                            % Do nothing
                        elseif redefineZ == 2
                            % There is no Z
                        else % For all other cases ask to redefine Z
                            redefineZ = 0;  
                        end
                    end
                    
                    % See if you need to rename the positional channels
                    if redefineX==1
                        VNames = strrep(VNames, redefinedX, 'X');
                    end
                    if redefineY==1
                        VNames = strrep(VNames, redefinedY, 'Y');
                    end
                    if redefineZ==1
                        VNames = strrep(VNames, redefinedZ, 'Z');
                    end
                    % remove duplicate channel names
                    valid_cols = Helper.remove_duplicates(VNames);
                    
                    % Convert to a table
                    if iscell(dat00.data)
                        dtemp = cell2mat(cellfun(@str2double,dat00.data,'un',0));
                        INDNAN = any(isnan(dtemp));
                        if all(INDNAN) % nothing converted from string, assume it imported weird
                            dtemp = array2table(dat00.data);
                            nmstmp = dtemp.Properties.VariableNames;
                            for col_i = find(colsIND)
                                dtemp.(nmstmp{col_i}) = cell2mat(dtemp.(nmstmp{col_i}));
                            end
                        else % the conversion from string worked
                            INDNAN = find(INDNAN);
                            dtemp = array2table(dtemp);
                            nmstmp = dtemp.Properties.VariableNames;
                            for col_i = INDNAN
                                dtemp.(nmstmp{col_i}) = dat00.data(:, col_i);
                            end
                        end
                    else
                        dtemp = array2table(dat00.data);
                    end
           
                    % Force numerical values in table
                    % ignore string columns in cell dat
                    IND_numeric = any(cellfun(@isnumeric, table2cell(dtemp)),1);
                    SubChNames = dtemp(:,IND_numeric).Properties.VariableNames;
                    % Force Numerical columns to be the same format
                    dtemp_fmt = array2table(cell2mat(table2cell(dtemp(:,IND_numeric))));
                    dtemp_fmt.Properties.VariableNames = SubChNames;
                    for chnm_i = 1:numel(SubChNames)
                        dtemp.(SubChNames{chnm_i}) = dtemp_fmt.(SubChNames{chnm_i});
                    end
          
                    % make a table for each cell/gate type
                    dat.(pnm_name).GatesTMP.(gate) = dtemp;
                    dat.(pnm_name).GatesTMP.(gate) = dat.(pnm_name).GatesTMP.(gate)(:, valid_cols);
                    dat.(pnm_name).GatesTMP.(gate).Properties.VariableNames = VNames(valid_cols);

                    % If I can't find X have user define them
                    if ~ismember('X', dat.(pnm_name).GatesTMP.(gate).Properties.VariableNames)
                        VNamestmp = dat.(pnm_name).GatesTMP.(gate).Properties.VariableNames;
                        indx = listdlg('PromptString', 'Select the X axis:','SelectionMode','multi', 'ListString', VNamestmp);
                        if size(indx,2) ~=1                            
                            dat.(pnm_name).GatesTMP.(gate).X = mean(dat.(pnm_name).GatesTMP.(gate){:, VNamestmp(indx)}, 2);
                        else
                            redefineX = 1;
                            redefinedX = VNamestmp{indx};
                            VNamestmp{indx} = 'X';
                            dat.(pnm_name).GatesTMP.(gate).Properties.VariableNames = VNamestmp;
                        end
                    end

                    % If I can't find Y have user define them
                    if ~ismember('Y', dat.(pnm_name).GatesTMP.(gate).Properties.VariableNames)
                        redefineY = 1;
                        VNamestmp = dat.(pnm_name).GatesTMP.(gate).Properties.VariableNames;
                        indy = listdlg('PromptString', 'Select the Y axis:','SelectionMode','multi', 'ListString', VNamestmp);
                        if size(indy,2) ~=1                            
                            dat.(pnm_name).GatesTMP.(gate).Y = mean(dat.(pnm_name).GatesTMP.(gate){:, VNamestmp(indy)}, 2);
                        else
                            redefineY = 1;
                            redefinedY = VNamestmp{indy};
                            VNamestmp{indy} = 'Y';
                            dat.(pnm_name).GatesTMP.(gate).Properties.VariableNames = VNamestmp;
                        end
                    end

                    % If I can't find Z eiher have user define them or create an empty channel
                    if redefineZ==2
                        dat.(pnm_name).GatesTMP.(gate).Z = 0.*dat.(pnm_name).GatesTMP.(gate).X;
                    end
                                        
                    if ~ismember('Z', dat.(pnm_name).GatesTMP.(gate).Properties.VariableNames)
                        VNamestmp = dat.(pnm_name).GatesTMP.(gate).Properties.VariableNames;
                        indz = listdlg('PromptString', 'Select the Z axis:','SelectionMode','multi', 'ListString', [{'There is no Z (make a fake one, otherwise everything breaks)'},VNamestmp(:)']);
                        %if Z data does not exist add it in as a column of
                        %zeros so stuff doesn't break later
                        if size(indz,2) ~=1
                            indz = indz-1;
                            dat.(pnm_name).GatesTMP.(gate).Z = mean(dat.(pnm_name).GatesTMP.(gate){:, VNamestmp(indz)}, 2);
                        else
                            if indz==1
                                redefineZ = 2;
                                dat.(pnm_name).GatesTMP.(gate).Z = 0.*dat.(pnm_name).GatesTMP.(gate).X;
                            else
                                indz = indz-1;
                                redefineZ = 1;
                                redefinedZ = VNamestmp{indz};
                                VNamestmp{indz} = 'Z';
                                dat.(pnm_name).GatesTMP.(gate).Properties.VariableNames = VNamestmp;
                            end
                        end
                    end
                end % end loop through files
            end % end loop through folders
        end

        function func_load_execute(UIfig, app, dat, fnm, pnm, DName, t_dat, swch, numrep)
            % FUNC_LOAD_EXECUTE it is main backend function which actually
            % performs modification of app to include loaded samples.
            %
             
            %   - UIfig - uifigure - Figure used for user interface in
            %       loading data.
            %   - app - Instance of CytoMAP.
            %   - dat - struct - Intermediate result of app.data. Should
            %       have all samples wth GatesTMP field.
            %   - pnm - cell - Paths to folders containing .csv/.mat files.
            %   - fnm - cell - Names of files of .csv/.mat files within
            %       these folders.
            %   - DName - cell - Name of dataset to be inserted into
            %       app.data. Ignored if size of pnm >= 1.
            %   - t_dat - cell - Data from user interface table, which was
            %       used to redefine some names etc.
            %
            % Modifies:
            %   - app - Specifically app.data, which after this function
            %       call will include new datasets
            
            if nargin == 7
                swch = 'NaN';
                numrep = 0;
            else
                swch = swch.Value;
            end
            
            
            SNsOLD = fieldnames(dat);
            vPD2 = waitbar(0, 'Doing some stuff');

            for pnm_j = 1:numel(pnm)
                DName = t_dat{pnm_j, 3};
                DName = strrep(DName, ' ', '');

                CeNames = t_dat(:, 9);
                ChNames = t_dat(:, 13);

                % Put the sample metadata in
                dat.(SNsOLD{pnm_j}).MetaData = struct;
                dat.(SNsOLD{pnm_j}).MetaData.Group = t_dat{pnm_j, 1};
                dat.(SNsOLD{pnm_j}).MetaData.Annotation = t_dat{pnm_j, 2};

                % remove empty elements
                ChNames = ChNames(~cellfun('isempty',ChNames));
                CeNames = CeNames(~cellfun('isempty',CeNames));
                for i=1:numel(ChNames)
                    ChNames{i} = strrep(ChNames{i}, ' ', '');
                end
                % re-name the gates/cells
                CeNamesOld = fieldnames(dat.(SNsOLD{pnm_j}).GatesTMP);
                ChNamesOld = dat.(SNsOLD{pnm_j}).GatesTMP.(CeNamesOld{1}).Properties.VariableNames;
                % Loop through all of the cell names
                for i=1:numel(CeNamesOld)
                    % reorder channels, so they are in the same order
                    % in all cell phenotypes and exclude channels not in
                    % all phenotypes
                    
                    % Find the index of the selected cell name in the table
                    cell_i = find(strcmp(t_dat(:, 6), CeNamesOld{i}));
                    % TODO: find a way to keep these not in all channels
                    logic = ismember(ChNamesOld, dat.(SNsOLD{pnm_j}).GatesTMP.(CeNamesOld{i}).Properties.VariableNames);
                    tmpChNames = ChNamesOld(logic);
                    dat.(SNsOLD{pnm_j}).GatesTMP.(CeNamesOld{i}) = dat.(SNsOLD{pnm_j}).GatesTMP.(CeNamesOld{i})(:,tmpChNames);
                    % If the cell is selected 
                    if t_dat{cell_i, 5}
                        try
                            CeNames{cell_i} = strrep(CeNames{cell_i}, ' ', '');
                            dat.(SNsOLD{pnm_j}).Gates.(CeNames{cell_i}) = dat.(SNsOLD{pnm_j}).GatesTMP.(CeNamesOld{i});
                        catch
                            errordlg('Something went wrong. Usually it is channel names')
                            return
                        end
                        try
                            if size(ChNames, 1) ~= 1
                                ChNames = ChNames';
                            end
                            dat.(SNsOLD{pnm_j}).Gates.(CeNames{cell_i}).Properties.VariableNames = ChNames(logic);
                        catch
                            errordlg('Something went wrong. Usually it is channel names');
                            return
                        end
                    end
                end
                dat.(SNsOLD{pnm_j}) = rmfield(dat.(SNsOLD{pnm_j}), 'GatesTMP');

                % Pull the names of the different types of cells
                dat.(SNsOLD{pnm_j}).CNames = fieldnames(dat.(SNsOLD{pnm_j}).Gates);

                %% Make dat.AllCells
                %%% This list contains the x-y-z of all unique cells across
                %%% all of your gates. This is useful for knowing your
                %%% total number of cells and what percentage of total does
                %%% each individual sub-population make up.
                waitbar(0.25, vPD2, 'Checking data format');
                clear UnsortedDat DatSub Logic
                % Initialize the table of all cells
                dat.(SNsOLD{pnm_j}).AllCells = dat.(SNsOLD{pnm_j}).Gates.(dat.(SNsOLD{pnm_j}).CNames{1});
                % Alternativly, initialize a cell data structure
                for i=2:numel(dat.(SNsOLD{pnm_j}).CNames)
                    DatSub = dat.(SNsOLD{pnm_j}).Gates.(dat.(SNsOLD{pnm_j}).CNames{i});
                    SubNames = DatSub.Properties.VariableNames;
                    % Only use channels that are shared between the first
                    % gate and other gates in all cells
                    % Find the index of each channel shared with AllCells
                    IND = zeros(numel(ChNames), 1);
                    INDAll = zeros(numel(ChNames), 1);
                    for ch_i=1:numel(ChNames)
                        % If the channel is in both AllCells and the
                        % selected Gate, find the index
                        if sum(strcmp(ChNames{ch_i}, SubNames))==1
                            IND(ch_i) = find(strcmp(SubNames, ChNames{ch_i}));
                            INDAll(ch_i) = ch_i;
                        else %If the channel is not in the Gate delete that element of IND
                            IND(ch_i) = 0;
                            INDAll(ch_i) = 0;
                        end
                    end
                    IND(IND==0) = [];
                    INDAll(INDAll==0) = [];
                    % rearange your gate to be in the same order as the first channel
                    DatSub = DatSub(:, IND);
                    % Take out any channels that are not present in both
                    % AllCells and the gate
                    dat.(SNsOLD{pnm_j}).AllCells = dat.(SNsOLD{pnm_j}).AllCells(:,INDAll);
                    
                    % ignore string columns in cell dat
                    IND_numeric = any(cellfun(@isnumeric, table2cell(DatSub)),1);
                    
                    Logic = ismember(DatSub(:,IND_numeric), dat.(SNsOLD{pnm_j}).AllCells(:,IND_numeric));
                    if max(~Logic)==1 %if any cells are not already in AllCells, add them
                        dat.(SNsOLD{pnm_j}).AllCells = [dat.(SNsOLD{pnm_j}).AllCells; DatSub(~Logic,:)];
                    end
                    ChNames = dat.(SNsOLD{pnm_j}).AllCells.Properties.VariableNames;
                end
                % If you took out elements, redefine the channel names

                % Now AllCells, should contain an element for every cell in
                % your data, and only have channels that are present in all
                % of your gates. This means if one of your gates does not
                % have a channel, all cells will not have that channel
                % either
                waitbar(1, vPD2, 'Checking data format');

                %% Build a binary gate structure
                %%%% Building a table with 1/0 elements for which gates
                %%%% each cells belong to make it easier to count cells
                %%%% when I find the number of cells in each neighborhood.
                %%%% I realize this is also a little redundant.
                waitbar(0.5, vPD2, 'Checking your gates');
                clear UnsortedDat DatSub Logic
                dat.(SNsOLD{pnm_j}).SortedDat = dat.(SNsOLD{pnm_j}).AllCells(:, {'X', 'Y', 'Z'});
                
                % loop through the cell names
                for ChellName_i=1:numel(fieldnames(dat.(SNsOLD{pnm_j}).Gates))
                    CellNameTMP = dat.(SNsOLD{pnm_j}).CNames{ChellName_i};
                    
                    DatSub = dat.(SNsOLD{pnm_j}).Gates.(CellNameTMP);
                    SubChNames = DatSub.Properties.VariableNames;
                    
                    % re-arange the channels in the gate to be in the same
                    % order as AllCells, also ignore any extra channels in
                    % the gate
                    IND = zeros(numel(ChNames), 1);
                    for ch_n=1:numel(ChNames)
                        IND(ch_n) = find(strcmp(SubChNames, ChNames{ch_n}));
                    end
                    
                    % rearange your gate to be in the same order as AllCells
                    DatSub = DatSub(:, IND);
                    
                    % ignore string columns in cell dat
                    IND_numeric = any(cellfun(@isnumeric, table2cell(DatSub)),1);
                    
                    % Find which cells in AllCells belong to this specific gate
                    Logic = ismember(dat.(SNsOLD{pnm_j}).AllCells(:,IND_numeric), DatSub(:,IND_numeric));
                    
                    % Put a one in the column corresponding to the gate for
                    % each cell that belongs to that group
                    dat.(SNsOLD{pnm_j}).SortedDat.(CellNameTMP)  = Logic;
                end

                %% Update app properties
                dat.(SNsOLD{pnm_j}).FNames = fnm';
                dat.(SNsOLD{pnm_j}).Path = pnm;
                dat.(SNsOLD{pnm_j}).Surfaces = struct;

                % Load the data into the main data structure
                app.data.(DName) = dat.(SNsOLD{pnm_j});
                app.data.(DName) = csv2tree(app.data.(DName));
                % remove the Gate_0 tag
                if ismember({strcat(Constants.gate_tag, '0')}, app.data.(DName).AllCells.Properties.VariableNames)
                    app.data.(DName).AllCells = removevars(app.data.(DName).AllCells, ...
                                                    {strcat(Constants.gate_tag, '0')});
                end
                
                switch swch
                    case 'NaN'
                        
                    case 'Value:'
                        INDNUM = cellfun(@isnumeric, table2cell(app.data.(DName).AllCells));
                        INDNAN = isnan(app.data.(DName).AllCells{:, any(INDNUM)});
                        app.data.(DName).AllCells{:, any(INDNUM)}(INDNAN) = numrep;
                end
            
            
            end
            app.DataN.Items = fieldnames(app.data);
            app.DataN.Value = DName;
            app.PhList.Name = app.DataN.Value;

            if ~isempty(UIfig) && isvalid(UIfig)
                close(UIfig)
            end

            waitbar(1, vPD2, 'Done!');
            close(vPD2);
            function treed_struct = csv2tree(csv_struct)
                % CSV2TREE Takes a app.data.X struct from just loaded .csv file, and converts it to style that .wsp files are loaded in.
                % Additionally checks for version, so that if .csv file was already new and in style of .wsp, then no redundant jobs are done.

                % If it's not legacy version, then use it's data format.
                % This assumes all versions after legacy have structure similar to the .wsp/.mat files.
                
                if isfield(csv_struct.MetaData, 'Ver') && ~strcmp(csv_struct.MetaData.Ver, Constans.LEGACY_VER)
                    treed_struct = csv_struct;
                else
                    treed_struct = struct;
                    treed_struct.FNames = csv_struct.FNames;
                    treed_struct.Path = csv_struct.Path{1, 1};
                    treed_struct.AllCells = csv_struct.AllCells;

                    if isfield(csv_struct, 'MetaData')
                        treed_struct.MetaData = csv_struct.MetaData;
                    else
                        treed_struct.MetaData = struct;
                        treed_struct.MetaData.Group = 1;
                        treed_struct.MetaData.Annotation = 0;
                    end
                    treed_struct.MetaData.Ver = Constants.CURR_VER;
                end

                % Check for non-basic fields. If they do not exist, either create stubs, or ignore.
                if isfield(csv_struct, 'Surfaces')
                    treed_struct.Surfaces = csv_struct.Surfaces;
                else
                    treed_struct.Surfaces = struct;
                    treed_struct.Surfaces.UDS = struct;
                    treed_struct.Surfaces.CCN = struct;
                    treed_struct.Surfaces.RSN = struct;
                end

                % If tree and GateTags already exist, there is no need to repeat the process.
                if isfield(treed_struct, 'tree') && isfield(treed_struct, 'GateTags')
                    return
                end

                treed_struct.GateTags = table;
                gate_names = csv_struct.SortedDat.Properties.VariableNames(~ismember(csv_struct.SortedDat.Properties.VariableNames, {'X' 'Y' 'Z'}));
                main_tree = tree('All', 'tag', [Constants.gate_tag '0']);
                tag = 1;
                for gate_kid=gate_names
                    kid = tree(gate_kid{1, 1}, 'gate_type', 'logic', 'gate_points', csv_struct.SortedDat(:, gate_kid), 'tag', sprintf([Constants.gate_tag '%d'], tag));
                    main_tree = main_tree.add_kid(kid);
% %                     treed_struct.GateTags.(sprintf([Constants.gate_tag '%d'], tag)) = strcat(main_tree.name, '/', kid.name);

                    % Add the full name
                    FullName = strcat(main_tree.name, '/', kid.name);
                    % Add only the immediate parent and kid name
                    ShortName = strsplit(main_tree.name, '/');
                    ShortName = [ShortName{end} '/' kid.name];
                    ShortName = strrep(strrep(strrep(strrep(ShortName, 'Gate_', ''), '_', ' '), 'NEG', '-'), 'POS', '+');

                    treed_struct.GateTags.(sprintf([Constants.gate_tag '%d'], tag)) = {FullName; ShortName};
                    tag = tag + 1;
                end
                treed_struct.AllCells = Helper.reorder_cols(treed_struct.AllCells);
                treed_struct.AllCells.Properties.VariableNames = Helper.valid_channel(treed_struct.AllCells.Properties.VariableNames);
                treed_struct.AllCells.Properties.VariableNames = Helper.valid_other(treed_struct.AllCells.Properties.VariableNames);
                treed_struct.AllCells = main_tree.gate_cells(treed_struct.AllCells);
                treed_struct.tree = main_tree;
            end

        end

        function [pathname] = uigetdir2(start_path, dialog_title)
            % UIGETDIR2 Allows user to load multiple folders/files in one
            % window. This overrides default MatLAB behaviour which allows
            % for only one folder.
            %
            % Input:
            %   - start_path - char, string - path at which to start dialog
            %       on.
            %   - dialog_title - char, string - Title of the window.
            %
            % Output:
            %   - pathname - cell - Paths to all the currently chosen
            %       files, and folders.
            
            import javax.swing.JFileChooser; % I am not sure exactly how this works
            if nargin == 0 || isempty(start_path) || start_path == 0 % Allow a null argument.
                start_path = pwd;
            end
            jchooser = javaObjectEDT('javax.swing.JFileChooser', start_path);
            jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            if nargin > 1
                jchooser.setDialogTitle(dialog_title);
            end
            jchooser.setMultiSelectionEnabled(true);
            status = jchooser.showOpenDialog([]);
            if status == JFileChooser.APPROVE_OPTION
                jFile = jchooser.getSelectedFiles();
                pathname{size(jFile, 1)}=[];
                for i=1:size(jFile, 1)
                    pathname{i} = char(jFile(i).getAbsolutePath);
                end
            elseif status == JFileChooser.CANCEL_OPTION
                pathname = [];
            else
                error('Error occured while picking file.');
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save/Export data functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function func_save(app, file)
            % FUNC_SAVE Exports data from current CytoMAP to either .csv or
            % .mat file.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - file - string - Name of file to which data should be
            %       saved. If not provided user will be asked to choose a
            %       file.
            %
            % Writes:
            %   - file - Creates new file, which depending on filetype
            %       (either .mat or .csv) will contain either all sensible
            %       information about this instance of CytoMAP (.mat) or
            %       only AllCells from default sample (i.e.
            %       app.DataN.Value) (.csv).
            
            if ~exist('file', 'var')  % Use most recent dataset
                [file,path] = uiputfile({'*.mat'; '*.csv'; '*.*'}, 'Save file name');
                if isempty(file)
                    return;
                end
                cd(path);
            end
            vPD = waitbar(0.25, 'Saving file');
            vPD.Position(4) = 85;
            if strcmp(file((end-3):end), '.mat')
                dat = struct;
                dat.data = app.data;
                dat.net = app.net;
                dat.map = app.map;
                dat.points = app.points;
                dat.polygons = app.polygons;
                dat.DataNItems = app.DataN.Items;
                dat.CellInfo = app.CellInfo;
                if ismember('NewGate', fields(app))
                    % "Gates Detected"
                    dat.NewGate = app.NewGate;
                else
                    % "No Gates Detected"
                end
                waitbar(0.5, vPD, 'Writing File');
                save(file, 'dat', '-v7.3', '-nocompression');
            elseif strcmp(file((end-3):end), '.csv')
                dat = app.data.(app.DataN.Value).AllCells;
                writetable(dat, file);
            end
            if isvalid(vPD)
                close(vPD);
            end
        end

        function func_ExportCSV(app, web)
            % FUNC_EXPORTCSV Allows user to export current state of
            % instance of CytoMAP to csv files. While similar to func_save,
            % it gives user more control over what specifically to save.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %
            % Writes:
            %   - path - Creates a new file in given path, with user chosen
            %       data (processed, and filtered in also way chosen by
            %       user) in it.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            tbldat = Helper.populate_table(app, 'smpls', {app.DataN.Value}, 'MFI', 'AllCells');
            tbldat = [tbldat; tbldat(end,:)];
            tbldat(:,1) = [];
            tbldat(:,1) = {false};
            tbldat(:,3) = [];
            tbldat(end,4) = {'Select All'};
            tbldat(end,3) = {false};
            
            tbldat(end,1) = {false};
            tbldat(end,2) = {'Select All'};
            tbldat(end,5) = {false};
            tbldat(end,6) = {'Select All'};

            UIfig = uifigure('Name', 'Export Data');
            UIfig.Resize = 'on';   
            UIfig.Scrollable = 'on';
            if web==1
                UIfig.Visible='OFF';
            end
            UIfig.Position = alpha*[10 10 750 800];
            t = uitable(UIfig);
            t.Data = tbldat;
            t.ColumnWidth = {alpha*50 alpha*200 alpha*50 alpha*200 alpha*50,alpha*200};
            t.Position = alpha*[0 50 750 750];
            t.ColumnName = {'', 'Include Cell Types', '', 'Include Channels', '', 'Export csv for Samples'};
            t.ColumnEditable = true;
            t.RowName = ([]);

            % Select a specific region to export
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Export only Regions: (e.g.: 1,2,6):');
            lbl.FontColor =app.GUIOPTS.txtclr;
            lbl.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
            lbl.Position = alpha*[10+245 35 400 15];
            
            RegSel = uieditfield(UIfig,'text');
            RegSel.Value = '';
            RegSel.Position = alpha*[10+240+100, 5, 100, 30];
            RegSel.BackgroundColor = app.GUIOPTS.bgclr;

            Model = uidropdown(UIfig);
            Model.Position = alpha*[10+240 5 100 30];
            Model.BackgroundColor = app.GUIOPTS.bgclr;

            if ~isempty(fieldnames(app.net))
                Model.Items = fieldnames(app.net);
                Model.Value = Model.Items{1};
            end

            % Select Data Type to Export
            lbl = uilabel(UIfig);
            lbl.Text = sprintf('Select Data Type to Export:');
            lbl.FontColor = app.GUIOPTS.txtclr;
            lbl.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
            lbl.Position = alpha*[10 35 400 15];
            DataType = uidropdown(UIfig);
            DataType.Items = {'Individual Cells', ...
                              'Raster Scanned Neighborhoods', ...
                              'Cell Centered Neighborhoods',};
            DataType.Value = {'Individual Cells'};
            DataType.Position = alpha*[10 5 235 30];
            DataType.BackgroundColor = app.GUIOPTS.bgclr;

            t.CellEditCallback = @(dd, p) switched_sample(app, DataType.Value, dd, p);
            DataType.ValueChangedFcn = @(~, ~) ChangeDataType(app, t, DataType);

            
% % %             % Create an include gate binaries option
% % %             lbl = uilabel(UIfig);
% % %             lbl.Text = sprintf('Include Gate Logicals');
% % %             lbl.FontColor = app.GUIOPTS.txtclr;
% % %             lbl.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
% % %             lbl.Position = alpha*[750-145-135 35 400 15];
% % %             gateL = uiswitch(UIfig, 'Orientation', 'horizontal');
% % %             gateL.Items = {'on','off'};
% % %             gateL.Value = gateL.Items{2};
% % %             gateL.Position = alpha*[750-145-115 5 150 30];
            
                        % Create a plot types button group
            gateL = uibuttongroup(UIfig, 'Visible','off');
            gateL.Position = [750-145-150, 5, 165, 45];
            Helper.func_SetCLR(app, gateL, 'UICpopup')
            % Create two radio buttons in the button group.
            r1 = uicheckbox(gateL, 'Text','Individual .csv for each cell',...
                'Position',[5, 2.5, 150, 20], 'Value', true);
            Helper.func_SetCLR(app, r1, 'table')
            r2 = uicheckbox(gateL, 'Text','Include Gate Logicals',...
                'Position',[5, 22.5, 150, 20], 'Value', false);
            Helper.func_SetCLR(app, r2, 'table')
            gateL.Visible = 'on';
            
            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) wrap_backend(t, app, DataType, RegSel, Model, gateL));
            btn.Position = alpha*[750-120, 5, 100, 40];
            btn.Text = 'Export';
            btn.BackgroundColor = app.GUIOPTS.bgclr;
            
            function ChangeDataType(app, t, DataType)
                % changes visiblity of manual choice numerical input.
                p.Indices = [5 5];
                p.EditData = 0;
                switched_sample(app, DataType.Value, t, p);
            end
            
            function switched_sample(app, DataType, dd, p)
                if p.Indices(2) == 5
                    % Make sure that at least one thing is selected
                    if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 5)), 5)))
                        dd.Data{p.Indices(1), 5} = true;
                        return;
                    end
                    % Process Select All button
                    fill_select_all = false;
                    if p.Indices(1) == find(strcmp(dd.Data(:, 6), 'Select All')) && p.Indices(2) == 5
                        ind = ~cellfun('isempty', dd.Data(:, 5));
                        dd.Data(ind, 5) = {logical(p.NewData)};
                        if ~logical(p.NewData)
                            dd.Data{strcmp(dd.Data(:, 6), app.DataN.Value), 5} = true;
                        end
                        fill_select_all = logical(p.NewData);
                    end
                    if strcmp(DataType, 'Individual Cells')
                        scan_type = 'AllCells';
                    elseif strcmp(DataType, 'Cell Centered Neighborhoods')
                        scan_type = 'MFICCN';
                    elseif strcmp(DataType, 'Raster Scanned Neighborhoods')
                        scan_type = 'MFIRSN';
                    end

                    ind = ~cellfun('isempty', dd.Data(:, 5));
                    ind(end) = false;
                    ind(ind) = logical(cell2mat(dd.Data(ind, 5)));
                    smpls = dd.Data(ind, 6);
                    if isempty(smpls)
                        smpls = {app.DataN.Value};
                    end
                    
                    tmpDat = cell(size(dd.Data, 1)-1, 8);

                    tmpDat(1:end, 1) = {1}; % Stub. Needed for populate table
                    tmpDat(1:end, 2) = dd.Data(1:end-1, 1);
                    tmpDat(1:end, 3) = dd.Data(1:end-1, 2);
                    tmpDat(1:end, 4) = {1}; % Stub. Needed for populate table
                    tmpDat(1:end, 5) = dd.Data(1:end-1, 3);
                    tmpDat(1:end, 6) = dd.Data(1:end-1, 4);
                    tmpDat(1:end, 7) = dd.Data(1:end-1, 5);
                    tmpDat(1:end, 8) = dd.Data(1:end-1, 6);
                    
                    tdataTMP = Helper.populate_table(app, ...
                        'smpls', smpls, ...
                        'MFI', scan_type, ...
                        'prev_table', tmpDat ...
                    );

                    tdataTMP(:,1) = [];
                    tdataTMP(:,3) = [];
                    tdataTMP(end + 1, :) = {false, 'Select All', false, 'Select All', false, 'Select All'}; 
                    dd.Data = tdataTMP;

                elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
                    ind = ~cellfun('isempty', dd.Data(:, 2));
                    dd.Data(ind, 1) = {logical(p.NewData)};
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 4), 'Select All')) && p.Indices(2) == 3
                    ind = ~cellfun('isempty', dd.Data(:, 4));
                    dd.Data(ind, 3) = {logical(p.NewData)};
                end
            end
            
            function wrap_backend(t, app, DataType, RegSel, Model, gateL)
                path = uigetdir(cd, 'Save file Path');
                if path==0
                    return
                end
                tempdat = t.Data(1:end-1, :);
                IO.func_ExportCSV_backend( ...
                    app, ...
                    path, ...
                    tempdat([tempdat{:,5}]==1, 6), ... Samples
                    tempdat([tempdat{:,1}]==1, 2), ... Phenotypes
                    tempdat([tempdat{:,3}]==1, 4), ... Channels
                    DataType.Value, ...
                    RegSel.Value, ...
                    Model.Value, ...
                    [gateL.Children.Value] ...
                )
            end
        end
        
        function func_ExportCSV_backend(app, path, smpl, CellNames, ChNames, DataType, RegSel, Model, gateL)
            % FUNC_EXPORTCSV_BACKEND it is back-end of function which allows
            % user to export current state of instance of CytoMAP to csv
            % files. While similar to func_save, it gives user more control
            % over what specifically to save.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - path - char, string - Path to which save file.
            %   - smpl - cell - Samples from which to save data
            %   - CellNames - cell - Names of Phenotypes (in full version),
            %       on which to gate data for saving.
            %   - ChNames - cell - Names of Channels which to save.
            %   - DataType - string - Wherther to save MFIRSN or MFICCN.
            %   - RegSel - char, int - Either '' or integer, which
            %       correspond to include which region to filter on from
            %       chosen model. Ignored if EXPcells is true.
            %   - Model - char, string - Name of model on which to look for
            %       regions to filter on. Ignored if EXPcells is true.
            %   - EXPcells - boolean - Whether to save AllCells or MFIs.
            %
            % Writes:
            %   - path - Creates a new file in given path, with given data
            %       in it.
            if isempty(ChNames) || isempty(smpl) || isempty(CellNames)
                return;
            end
            if ~isempty(RegSel)
                RegSel = str2double(split(RegSel, ','));
            end
            
            current_path = pwd;
            cd(path);

            for smpl_i = 1:numel(smpl)
                Ch = Helper.valid_channel(ChNames);
                GateNames = Helper.gate_full2tag(app, CellNames, smpl{smpl_i});
                switch DataType
                    case 'Individual Cells'
                        dat = app.data.(smpl{smpl_i}).AllCells;
                        Gates = app.data.(smpl{smpl_i}).AllCells(:,GateNames);
                        % Include the gate logicals?
                        if gateL(1)==1
                            Ch = [Ch; GateNames];
                        end
                    case 'Raster Scanned Neighborhoods'
                        dat = app.data.(smpl{smpl_i}).MFIRSN;
                        Gates = app.data.(smpl{smpl_i}).MFIRSN(:,GateNames);
                        % remove empty neighborhoods
                        Gates(dat.NCells==0, :) = [];
                        dat(dat.NCells==0, :) = [];
                        
                        Ch = [Ch; GateNames];
                        nms = [Helper.valid_var(ChNames); Helper.valid_var(CellNames)];
                    case 'Cell Centered Neighborhoods'
                        dat = app.data.(smpl{smpl_i}).MFICCN;
                        Gates = app.data.(smpl{smpl_i}).MFICCN(:,GateNames);
                        
                        Gates(dat.NCells==0, :) = [];
                        dat(dat.NCells==0, :) = [];
                        Ch = [Ch; GateNames];
                        nms = [Helper.valid_var(ChNames); Helper.valid_var(CellNames)];
                end
                
                % Reorder columns, so X, Y, Z are always first.
                if ismember(strcat(Constants.channel_tag, 'chID'), dat.Properties.VariableNames)
                   dat = movevars(dat, 'X', 'After', strcat(Constants.channel_tag, 'chID'));
                elseif ismember(strcat(Constants.channel_tag, 'ID'), dat.Properties.VariableNames)
                   dat = movevars(dat, 'X', 'After', strcat(Constants.channel_tag, 'ID'));
                end
                dat = movevars(dat, 'Y', 'After', 'X');
                dat = movevars(dat, 'Z', 'After', 'Y');
                
                % Pull only the selected channels
                dat = dat(:,Ch);

                %% Write Table
                    
                % If the user selected specific cells
                if strcmp(DataType, 'Individual Cells')
                    if gateL(2)==1
                        for cell_i=1:numel(CellNames)
                            % Pull only the selected cell types
                            datSUB = dat(Gates.(GateNames{cell_i})==1, :);
                            if ~isempty(RegSel)
                                if numel(RegSel)==1
                                    datSUB = datSUB(datSUB.([Constants.other_tag, Model])==RegSel, :);
                                else
                                    datSUB = datSUB(sum(datSUB.([Constants.other_tag, Model])==RegSel', 2)==1, :);
                                end
                            end
                            % Build the file name
                            fname = ['CytoMAP_' smpl{smpl_i} 'Gate__' strrep(CellNames{cell_i}, '/', '_') '.csv'];
                            % Make the channel names better
                            names = datSUB.Properties.VariableNames;
                            names = strrep(strrep(names, Constants.other_tag, ''), Constants.channel_tag, '');
                            if gateL(1)==1
                                names((end-numel(CellNames)+1):end) = CellNames; 
                            end
                            datSUB.Properties.VariableNames = Helper.valid_var(names);
                            writetable(datSUB,fname);
                        end
                    else
                        % Pull only the selected cell types
                        datSUB = dat;
                        % build logical
                        cell_logic = zeros(size(dat,1),1);
                        for cell_i=1:numel(CellNames)
                            cell_logic = (Gates.(GateNames{cell_i}) + cell_logic) ~=0;
                        end
                        datSUB = datSUB(cell_logic==1, :);
                        if ~isempty(RegSel)
                            if numel(RegSel)==1
                                datSUB = datSUB(datSUB.([Constants.other_tag, Model])==RegSel, :);
                            else
                                datSUB = datSUB(sum(datSUB.([Constants.other_tag, Model])==RegSel', 2)==1, :);
                            end
                        end
                        % Build the file name
                        fname = ['CytoMAP_' smpl{smpl_i} '.csv'];
                        % Make the channel names better
                        names = datSUB.Properties.VariableNames;
                        names = strrep(strrep(names, Constants.other_tag, ''), Constants.channel_tag, '');
                        if gateL(1)==1
                            names((end-numel(CellNames)+1):end) = CellNames; 
                        end
                        datSUB.Properties.VariableNames = Helper.valid_var(names);
                        writetable(datSUB,fname);
                    end
                    
                else
                    % Pull all neighborhoods
                    datSUB = dat;

                    if ~isempty(RegSel)
                        if numel(RegSel)==1
                            datSUB = datSUB(datSUB.([Constants.other_tag, Model])==RegSel, :);
                        else
                            datSUB = datSUB(sum(datSUB.([Constants.other_tag, Model])==RegSel', 2)==1, :);
                        end
                    end

                    fname = ['CytoMAP_' smpl{smpl_i} '_Neighborhoods' '.csv'];

                    % Make the channel names better
                    names = nms;
                    names = strrep(strrep(names, Constants.other_tag, ''), Constants.channel_tag, '');
                    datSUB.Properties.VariableNames = Helper.valid_var(names);

                    writetable(datSUB,fname);
                            
                end
            end
            cd(current_path)
        end

        function func_ExportPrism(app, web)
            % FUNC_EXPORTPRISM Allows user to export data in the csv format
            % which is compatible with prism.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %
            % Writes:
            %   - Creates files in the directory chosen by user, and with
            %       contents determined by user. One file per chosen
            %       sample.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            if ~Helper.any_sample(app)
                return;
            end
            
            [~, sample] = Helper.find_MFI(app);
            if ~iscell(sample)
                sample = {sample};
            end
            
            DataType = 'Individual Cells';
            tbldat = Helper.populate_table(app, 'smpls', sample, 'MFI', 'AllCells', 'fill_checkbox', false);
            tbldat = [tbldat; tbldat(end,:)];
            tbldat(:,1) = []; tbldat(:,3) = [];
            tbldat(end + 1, :) = {false, 'Select All', false, 'Select All', false, 'Select All'}; 

            UIfig = uifigure('Name', 'Export Data');
            UIfig.Resize = 'on';   
            UIfig.Scrollable = 'on';
            if web==1
                UIfig.Visible='OFF';
            end
            UIfig.Position = alpha*[10 10 750 800];
            
            t = uitable(UIfig);
% % %             s1 = uistyle('BackgroundColor',app.GUIOPTS.bgclr, 'FontColor',app.GUIOPTS.txtclr,'FontName',app.GUIOPTS.FontName);
% % %             addStyle(t,s1);
            t.Data = tbldat;
            t.ColumnWidth = {alpha*50 alpha*200 alpha*50 alpha*200 alpha*50,alpha*200};
            t.Position = alpha*[0 50 750 750];
            t.ColumnName = {'Export','Include Cell Types','', 'Include Channels','','Export csv for Samples'};
            t.ColumnEditable = true;
            t.RowName = ([]);
            t.CellEditCallback = @(dd, p) switched_sample(app, DataType, dd, p);

            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) wrap_backend(t, app));
            btn.Position = alpha*[10, 5, 100, 40];
            btn.Text = 'Ok';
            
            function switched_sample(app, DataType, dd, p)
                if p.Indices(2) == 5
                    % Make sure that at least one thing is selected
                    if p.EditData == 0 && ~any(cell2mat(dd.Data(~cellfun('isempty', dd.Data(:, 5)), 5)))
                        dd.Data{p.Indices(1), 5} = true;
                        return;
                    end
                    % Process Select All button
                    fill_select_all = false;
                    if p.Indices(1) == find(strcmp(dd.Data(:, 6), 'Select All')) && p.Indices(2) == 5
                        ind = ~cellfun('isempty', dd.Data(:, 5));
                        dd.Data(ind, 5) = {logical(p.NewData)};
                        if ~logical(p.NewData)
                            dd.Data{strcmp(dd.Data(:, 6), app.DataN.Value), 5} = true;
                        end
                        fill_select_all = logical(p.NewData);
                    end
                    if strcmp(DataType, 'Individual Cells')
                        scan_type = 'AllCells';
                    elseif strcmp(DataType, 'Cell Centered Neighborhoods')
                        scan_type = 'MFICCN';
                    elseif strcmp(DataType, 'Raster Scanned Neighborhoods')
                        scan_type = 'MFIRSN';
                    end

                    ind = ~cellfun('isempty', dd.Data(:, 5));
                    ind(end) = false;
                    ind(ind) = logical(cell2mat(dd.Data(ind, 5)));
                    smpls = dd.Data(ind, 6);
                    if isempty(smpls)
                        smpls = {app.DataN.Value};
                    end
                    
                    tmpDat = cell(size(dd.Data, 1)-1, 8);

                    tmpDat(1:end, 1) = {1}; % Stub. Needed for populate table
                    tmpDat(1:end, 2) = dd.Data(1:end-1, 1);
                    tmpDat(1:end, 3) = dd.Data(1:end-1, 2);
                    tmpDat(1:end, 4) = {1}; % Stub. Needed for populate table
                    tmpDat(1:end, 5) = dd.Data(1:end-1, 3);
                    tmpDat(1:end, 6) = dd.Data(1:end-1, 4);
                    tmpDat(1:end, 7) = dd.Data(1:end-1, 5);
                    tmpDat(1:end, 8) = dd.Data(1:end-1, 6);
                    
                    tdataTMP = Helper.populate_table(app, ...
                        'smpls', smpls, ...
                        'MFI', scan_type, ...
                        'prev_table', tmpDat ...
                    );

                    tdataTMP(:,1) = [];
                    tdataTMP(:,3) = [];
                    tdataTMP(end + 1, :) = {false, 'Select All', false, 'Select All', false, 'Select All'}; 
                    dd.Data = tdataTMP;

                elseif p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
                    ind = ~cellfun('isempty', dd.Data(:, 2));
                    dd.Data(ind, 1) = {logical(p.NewData)};
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 4), 'Select All')) && p.Indices(2) == 3
                    ind = ~cellfun('isempty', dd.Data(:, 4));
                    dd.Data(ind, 3) = {logical(p.NewData)};
                end
            end
            
            function wrap_backend(t, app)
                path = uigetdir(cd, 'Save file Path');
                if isempty(path)
                    return
                end
                tmp_data = t.Data(1:end-1, :);
                IO.func_ExportPrism_backend( ...
                    app, ...
                    path, ...
                    tmp_data([tmp_data{:,5}]==1, 6), ... Samples
                    tmp_data([tmp_data{:,1}]==1, 2), ... Phenotypes
                    tmp_data([tmp_data{:,3}]==1, 4) ... Channel Names
                )
            end
        end

        function func_ExportPrism_backend(app, path, smpl, CellNames, ChNames)
            % FUNC_EXPORTPRISM_BACKEND Back-end of Export Data to Prism. Allows user to export data
            % in the csv format which is compatible with prism.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - path - char, string - Path to folder where to save files.
            %   - smpl - cell - Samples from which files will be created.
            %   - CellNames - cell - Phenotypes to be included in files.
            %   - ChNames - cell - Channels to be included in files.
            %
            % Writes:
            %   - Creates files in the directory chosen by user, and with
            %       contents determined by user. One file per given sample.
            
            if isempty(smpl) || isempty(CellNames) || isempty(ChNames)
                return;
            end
            
            Ch = Helper.valid_channel(ChNames);
            
            current_path = pwd;
            cd(path);

            for smpl_i = 1:numel(smpl)
                % Define the destination file name
                % Pull the gate Names for the selected cell types for
                % each sample
                GateNames = Helper.gate_full2tag(app,CellNames, smpl{smpl_i});
                dat = app.data.(smpl{smpl_i}).AllCells(:,Ch);
                Gates = app.data.(smpl{smpl_i}).AllCells(:,GateNames);
                datMAT = zeros(max(sum(table2array(Gates))),numel(GateNames)*numel(Ch));
                datMAT(datMAT==0) = NaN;
                Header = cell(2,numel(GateNames)*numel(ChNames));
                for cell_i=1:numel(GateNames)
                    datSUB = dat(Gates.(GateNames{cell_i})==1, :);
                    datSUB = table2array(datSUB);

                    st = (1+(cell_i-1)*numel(ChNames));
                    en = ((cell_i)*numel(ChNames));

                    Header{1,st} = CellNames{cell_i};
                    Header(2,st:en) = ChNames;
                    datMAT(1:size(datSUB,1), st:en) = datSUB;
% %                             if ~isempty(RegSel.Value)
% %                                 if strcmp(DataType.Value, 'Raster Scanned Neighborhoods')
% %                                     datSUB = datSUB(datSUB.RegionRSN==str2double(RegSel.Value), :);
% %                                 elseif strcmp(DataType.Value, 'Cell Centered Neighborhoods')
% %                                     datSUB = datSUB(datSUB.RegionCCN==str2double(RegSel.Value), :);
% %                                 end
% %                             end
                end % end of phenotype loop
                fname = ['CytoMAP_' smpl{smpl_i} '.csv'];
                fid = fopen(fname,'w');
                headerformat = '%s\n';
                dataformat = '%f\n';
                for i=1:(numel(GateNames)*numel(ChNames)-1)
                    headerformat = ['%s, ' headerformat];
                    dataformat = ['%f, ' dataformat];
                end
                fprintf(fid,headerformat,Header{1,:});
                fprintf(fid,headerformat,Header{2,:});
                fclose(fid);
                dlmwrite(fname, datMAT, '-append', 'precision', 5)
            end
            cd(current_path);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Phenotype Table Manipulations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function func_RemovePH(app, smpls, phns)
            % FUNC_REMOVEPH Removes phenotypes, and all their children from
            % all the samples given.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smpls - char, string, cell - Names of samples in which to
            %       remove given phenotypes.
            %   - phns - char, string, cell - Names of phenotypes to be
            %       removed.
            %
            % Modifies:
            %   - app - Specifically app.data.(smpl) for each smpl in given
            %       smpls. Almost all traces of phns will be removed (some
            %       might still be required for other parts of program to
            %       work).
            
            if ~iscell(phns)
                phns = {phns};
            end
            if ~iscell(smpls)
                smpls = {smpls};
            end

            for p=1:numel(phns)
                phn = phns{p};
                if strcmp(phn, 'Select All')
                    break
                end
                for s=1:numel(smpls)
                    spl = smpls{s};
                    if strcmp(spl, 'Select All')
                        break
                    end
                    % Find all the useful info
                    idx_in_tags = strcmp(phn, table2cell(app.data.(spl).GateTags(2, :)));
                    % If we were removing a parent of that phenotype first it already is removed.
                    if sum(idx_in_tags) == 0
                        continue;
                    end
                    tag = app.data.(spl).GateTags.Properties.VariableNames{idx_in_tags};
                    path = table2cell(app.data.(spl).GateTags(1, tag));

                    % Remove all of the children too, as they would not be in tree.
                    child_idxs = startsWith(table2cell(app.data.(spl).GateTags(1, :)), strcat(path, "/"));
                    tags = app.data.(spl).GateTags.Properties.VariableNames(child_idxs);
                    if ~ismember(tag, tags)
                        tags = [{tag}, tags]; %#ok<AGROW>
                    end

                    if ~iscell(tags)
                        tags = {tags};
                    end

                    % Remove from tables
                    app.data.(spl).AllCells = removevars(app.data.(spl).AllCells, tags);
                    app.data.(spl).GateTags = removevars(app.data.(spl).GateTags, tags);

                    % Remove any cells from AllCells who no longer have any
                    % phenotype annotation/gate tags
                    [~, short] = Helper.get_gates(app, spl);
                    AllTags = Helper.get_tag(app, short, spl);
                    INDz = any(table2array(app.data.(spl).AllCells(:, AllTags)), 2);
                    app.data.(spl).AllCells = app.data.(spl).AllCells(INDz,:);
                    try
                        % Remove from tree (it's sufficient to remove only the path,
                        % since it's kids will be also removed)
                        app.data.(spl).tree = app.data.(spl).tree.remove_kid(path);
                    catch
                        disp('Failed to remove elements from tree')
                    end
                end
            end
        end

        function func_RemoveSmpl(app, smpls)
            % FUNC_REMOVESMPL Removes samples from current instance of
            % CytoMAP.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smpls - char, string, cell - Names of samples to be
            %       removed.
            %
            % Modifies:
            %   - app - Specifically app.data, with every smpl in given
            %       smpls removed from struct. Additionally other fields in
            %       app might be affected too.
            
            for s=1:numel(smpls)
                smpl = smpls{s};
                if strcmp(smpl, 'Select All')
                    break
                end
                app.data = rmfield(app.data, smpl);
                if numel(app.DataN.Items) == 1
                    app.DataN.Items = {'No Data Loaded'}; %Initialize with no data loaded
                    app.DataN.Value = app.DataN.Items{1};
                else
                    app.DataN.Items = fieldnames(app.data);
                    app.DataN.Value = app.DataN.Items{1};
                end
            end
        end

        function func_AddPH(app, smpls)
            % FUNC_ADDPH Allows user to choose what phenotypes to add for
            % each sample.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smpls - cell - Samples to which add the phenotypes
            %       (phenotypes are chosen by user on per sample basis).
            %
            % Modifies:
            %   - app - Specifically app.data.(smpl) for each smpl in
            %       smpls, with each one added phenotypes. 
            
            for s=1:numel(smpls)
                smpl = smpls{s};
                % Pull in the data
                TEMPdat = app.data.(smpl);
                % Define the filepaths
                clear fnm pnm datt choice n
                %Input the Folder paths and File paths
                [fnm , pnm]=uigetfile({'*.csv', 'FlowJo files (*.csv)'}, ...
                    'Select either FLOWJO csv files that you want to add','Multiselect', 'on');
                %Convert character array to cell array
                if ischar(fnm)
                    fnm = {fnm};
                end
                if isempty(fnm) || isempty(pnm)
                    return;
                end
                % Import the data
                ch_names = Helper.get_channels(app, smpl);
                %% Load FlowJo .csv files
                for j = 1:size(fnm, 2)
                    %open the file
                    dat00=importdata([pnm fnm{j}]);
                    DName = app.DataN.Value;
                    % If there isn't a similar string in the filenames call this data, data1
                    gate = fnm{j}(1:end-4);

                    % Define the gate name to valid variable
                    gate = strrep(gate, 'export_', '');
                    gate = strrep(gate, 'AllCells_', '');
                    gate = strrep(gate, 'Sample_', '');
                    gate = Helper.valid_gate(gate);

                    % Figure out full name of a gate, it's tags and path
                    full_name = Helper.full_gate(gate);
                    ShortName = strcat('All/', full_name);
                    tag = Helper.get_tag(app, ShortName, smpl);
                    path = strcat('Gate_All', '/', gate);

                    %Clean up the variable names
                    VNames = strrep(dat00.colheaders, '-','');
                    VNames = Helper.valid_channel(VNames);

                    % Make sure new variable names are compatible with current sample
                    if any(~ismember(ch_names, VNames))
                        warndlg(strcat('File ', fnm{j}, ' is heavily not compatible with sample ', smpl, '.', newline, ...
                            'It will be skipped.'));
                        continue;
                    elseif any(~ismember(VNames, ch_names))
                        answer = questdlg(...
                                strcat('File ', fnm{j}, ' is somewhat not compatible with sample ', smpl, '.', newline, ...
                                'It can lead to some loss of channel information inside of current instance of MATLAB.', newline, ...
                                'Do you want to add it, skip this file/sample combination, or cancel adding any of chosen files?'), ...
                                'Incompatible Sample', 'Add', 'Skip', 'Cancel', 'Add');
                        if isempty(answer)
                            return;
                        end
                        switch answer
                            case 'Add'
                                % Keep going. This makes code less redundant
                            case 'Skip'
                                continue;
                            case 'Cancel'
                                return;
                            case ''
                                return;
                        end
                    end

                    % Make a table of new variables, compatible with current sample
                    VName_idxs = ismember(VNames, ch_names);
                    VNames = VNames(VName_idxs);
                    new_data = array2table(dat00.data(:, VName_idxs));
                    new_data.Properties.VariableNames = VNames;
                    new_data.Z = 0 .* new_data.X;

                    % Add rows to All Cells
                    chanels_only = TEMPdat.AllCells(VNames);
                    new_rows = ~ismember(new_data, chanels_only);
                    if ~isempty(new_rows) % Not commented behaviour in docs
                        new_data = new_data(new_rows, :);

                        new_data = TEMPdat.tree.gate_cells(new_data); % Gates new data points.
                        % Fill in new columns with zeros, since they can be distances/regions/new channels etc.
                        for c=1:numel(TEMPdat.AllCells.Properties.VariableNames)
                            if ismember({col}, VNames)
                                new_data.(col) = 0. * new_data.X;
                            end
                        end

                        % Add new tag to each table
                        TEMPdat.AllCells.(tag) = 0 .* TEMPdat.AllCells.X;
                        new_data.(tag) = 1 .* new_data.X;

                        % Merge two into 1 table
                        TEMPdat.AllCells = [TEMPdat.AllCells; ...
                                            new_data];
                    end

                    % Update GateTags
                    TEMPDat.GateTags.(tag) = {path; ShortName};

                    % Update tree
                    kid = tree(full_name, 'gate_type', 'logic', 'tag', tag);
                    TEMPDat.tree = TEMPDat.tree.add_kid(obj, kid);
                end

                % Update sample pointer
                app.data.(smpl) = TEMPdat;
            end
        end

        function func_MergePH(app, smpls, phns)
            % FUNC_MERGEPH Front-End of Merge Phenotypes function. Allows user to merge
            % multiple phenotypes together into one. It's useful for
            % creating gates which are union of two disjointed parts in
            % some dimensions. Old Phenotypes are not overwritten, but a
            % new phenotype is created with different name.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smpls - char, string, cell - Samples in which to merge
            %       phenotypes.
            %   - phns - cell (size > 1) - Phenotypes which should be
            %       merged.
            %
            % Modifies:
            %   - app - Specifically app.data.(smpl) for smpl in smpls,
            %       with added new phenotype, which is combination of all
            %       chosen ones.
            
            if ~iscell(phns) || numel(phns) <= 1
                warndlg(strcat('Number of phenotypes given needs to be more than 1.', ...
                        'If you want to rename your phenotype choose "Rename Phenotype" option.'))
                return;
            end
            prompt = {'Enter name of merged phenotype:'};
            dlgtitle = 'Choose name of the merged phenotype';
            dims = [1 35];
            definput = {''};
            new_name = inputdlg(prompt,dlgtitle,dims,definput);
            if isempty(new_name) || isempty(new_name{1})
                return;
            end

            [~, parents] = Helper.get_gates(app, smpls);

            if ~ismember('All Cells', parents)
                parents = [{'All Cells'}, parents];
            end

            done = false;
            while ~done
                [indx, tf] = listdlg(...
                    'PromptString', "Select a parent phenotype:", ...
                    'SelectionMode', 'single', ...
                    'InitialValue', find(strcmp(parents, "All Cells")), ...
                    'ListString', parents, ...
                    'ListSize', [350, 200], ...
                    'Name', 'Parent Phenotype selection' ...
                );
                if ~tf
                    return;
                end
                if indx ~= find(strcmp(parents, "All Cells"))
                    isdone = questdlg( ...
                            "Note that the merged phenotype has to be a subset of a parent you chose." + newline + ...
                                "Do you want to continue?", ...
                            "Confirm Parent", ...
                            "Yes", ...
                            "Go Back", ...
                            "Cancel", ...
                            "Yes" ...
                    );
                    switch isdone
                        case "Go Back"
                            % Done stays false
                        case "Yes"
                            done = true;
                        otherwise
                            return;
                    end
                else
                    done = true;
                end
            end
            IO.func_MergePH_backend(app, smpls, phns, new_name, parents{indx});
        end

        function func_MergePH_backend(app, smpls, phns, name, parent)
            % FUNC_MERGEPH_BACKEND Back-End of Merge Phenotypes function. Allows user to merge
            % multiple phenotypes together into one. It's useful for
            % creating gates which are union of two disjointed parts in
            % some dimensions. Old Phenotypes are not overwritten, but a
            % new phenotype is created with different name.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smpls - char, string, cell - Samples in which to merge
            %       phenotypes.
            %   - phns - cell (size > 1) - Phenotypes which should be
            %       merged.
            %   - name - char, string - Name of the new phenotype to be
            %       created.
            %   - parents - char, string - Name of parent under which to
            %       make this new phenotype (new phenotype has to be subset
            %       of the parents).
            %
            % Modifies:
            %   - app - Specifically app.data.(smpl) for smpl in smpls,
            %       with added new phenotype under given parent, which is
            %       combination of all chosen ones.
            
            if ~iscell(smpls)
                smpls = {smpls};
            end

            valid_name = Helper.valid_gate(name);

            for smpl_idx=1:numel(smpls)
                smpl = smpls{smpl_idx};
                combined_tags = Helper.gate_full2tag(app, phns, smpl);

                % Make path/short names
                if strcmp(parent, "All Cells")
                    path = strcat('Gate_All/', valid_name);
                    short = strcat('All/', name);
                else
                    parent_tag = Helper.gate_full2tag(app, parent, smpl);
                    path = Helper.gate_full2path(app, parent, smpl);
                    path = strcat(path{1}, '/', valid_name);
                    short = split(parent, '/');
                    short = strcat(short{2}, '/', name);
                end

                % Get a tag for new PH in this sample
                tag = char(Helper.get_tag(app, short, smpl));

                % Create logic for AllCells
                logic = zeros(size(app.data.(smpl).AllCells, 1), 1);
                for cmbd_tag_idx=1:numel(combined_tags)
                    logic = logic | app.data.(smpl).AllCells{:, combined_tags{cmbd_tag_idx}};
                end
                if ~strcmp(parent, "All Cells")
                    logic = logic & app.data.(smpl).AllCells{:, Helper.gate_full2tag(app, parent, smpl)};
                end

                if isfield(app.data.(smpl), 'MFIRSN')
                    % Create estimates of the each included tag in each region.
                    % Account for overlap etc.
                end

                if isfield(app.data.(smpl), 'MFICCN')
                end

                app.data.(smpl).GateTags{:, tag} = {char(path); char(short)};
                app.data.(smpl).AllCells{:, tag} = logic;
                app.data.(smpl).tree = app.data.(smpl).tree.add_kid(tree(name, 'tag', tag, 'gate_type', 'logic'));
            end
        end

        function func_MergeSmpls(app, smpls, new_smpl_name, rename)
            % FUNC_MERGESMPLS Allows user to merge multiple samples into
            % one. This is useful in case in which samples had to be
            % splitted due to memory constraints. It also can be used for
            % renaming already loaded sample.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - smpls - cell - Samples to be merged
            %   - new_smpl_name - char, string, cell (only 1st element
            %       used) - Name of new sample to be created from merged
            %       ones.
            %
            % Modifies:
            %   - app - Specifically adds app.data.(new_smpl_name),
            %       representing a new sample, which is a union of the ones
            %       given.
            
            if iscell(new_smpl_name)
                new_smpl_name = new_smpl_name{1};
            end
            if numel(smpls) <= 0
                errordlg('You need at least 1 sample to be merged.');
                return;
            elseif numel(smpls) == 1 && rename==1
                % If there is only one sample, just rename, since there is nothing to merge.
                app.data.(new_smpl_name) = app.data.(smpls{1});
                app.DataN.Items(strcmp(smpls, app.DataN.Items)) = {new_smpl_name};
                app.data = rmfield(app.data, smpls{1});
                if strcmp(app.DataN.Value, smpls)
                    app.DataN.Value = {new_smpl_name};
                end
                return;
            end

            vPD = waitbar(0, 'Merging Samples');
            vPD.Position(4) = 85;
            % Find the intersection of the two samples you are merging
            new_smpl_struct = Helper.intersect_smpls(app, smpls, vPD);
            % Check to make sure there is something to merge
            if isempty(new_smpl_struct)
                waitdlg('Combined Samples create an empty sample. Aborting');
                close(vPD);
                return;
            end
            % add the new merged sample to the data structure
            app.data.(new_smpl_name) = new_smpl_struct;
            if ~any(strcmp(new_smpl_name, app.DataN.Items))
                app.DataN.Items(end + 1) = {new_smpl_name};
            end
            close(vPD);
        end

        function func_NewSmpls(app, smpls, phns, new_smpl_name, phnsKeep)
            % FUNC_MERGESMPLS Allows user to merge multiple samples into
            % one. This is useful in case in which samples had to be
            % splitted due to memory constraints. It also can be used for
            % renaming already loaded sample.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - smpls - cell - Samples to be merged
            %   - new_smpl_name - char, string, cell (only 1st element
            %       used) - Name of new sample to be created from merged
            %       ones.
            %
            % Modifies:
            %   - app - Specifically adds app.data.(new_smpl_name),
            %       representing a new sample, which is a union of the ones
            %       given.
            
            IO.func_MergeSmpls(app, smpls, new_smpl_name, 0);
            IO.func_RemovePH(app, new_smpl_name, phns);
            
            % Clean Up All Cells
            gatesKeep = Helper.gate_full2tag(app, phnsKeep, new_smpl_name);
            INDkeep = table2array(app.data.(new_smpl_name{1}).AllCells(:, gatesKeep));
            if size(INDkeep, 2)>1
                INDkeep = any(INDkeep, 2);
            end
            app.data.(new_smpl_name{1}).AllCells = app.data.(new_smpl_name{1}).AllCells(INDkeep==1, :);
            
            
        end
        
        function func_Rename(app, smpl, old, new)
            % FUNC_RENAME Allows user to rename phenotypes within given
            % samples.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - smpl - char, string, cell - Sample in which to rename
            %       phenotype.
            %   - old - char, string - Phenotype to be renamed.
            %   - new - char, string - New name of that phenotype.
            %
            % Modifies:
            %   - app - More specifically app.data.(s) for s in smpl.
            %       Creates a copy of given phenotype with a new name. Note
            %       that the old phenotype is not removed.
            
            if ~iscell(smpl)
                smpl = {smpl};
            end
            
            parent = split(old, '/');
            parent = parent(1);

            new_valid = split(new, '/');
            new = new_valid(2);
            new_valid = Helper.valid_gate(new_valid(2));

            for s=1:numel(smpl)
                spl = smpl{s};
                idx_in_tags = strcmp(old, table2cell(app.data.(spl).GateTags(2, :)));
                old = split(old, '/');
                old = old(2);
                tag = app.data.(spl).GateTags.Properties.VariableNames{idx_in_tags};
                parent_tag = Helper.gate_tag2parent_tag(app, tag, spl);
                old_path = table2cell(app.data.(spl).GateTags(1, tag));
                if strcmp(parent_tag, strcat('Gate_0'))
                    parent_path = 'Gate_All';
                else
                    parent_path = app.data.(spl).GateTags{1, parent_tag};
                end

                % Define new path
                new_path = strcat(parent_path, '/', new_valid);

                % Find children to rename (but not the current one)
                child_idxs = startsWith(table2cell(app.data.(spl).GateTags(1, :)), path) & ...
                    ~strcmp(table2cell(app.data.(spl).GateTags(1, :)), path);
                tags = app.data.(spl).GateTags.Properties.VariableNames(child_idxs);
                for t=1:numel(tags)
                    % For each child, switch path, and name, if necessary
                    tg = tags{t};
                    kid_path = app.data.(spl).GateTags{1, tg};
                    kid_name = app.data.(spl).GateTags{2, tg};

                    kid_path = strrep(kid_path, path, new_path);
                    if startsWith(kid_name, strcat(old, '/'))
                        kid_name = strrep(kid_name, old, new);
                    end
                    if iscell(kid_path)
                        kid_path = kid_path{1};
                    end
                    if iscell(kid_name)
                        kid_name = kid_name{1};
                    end
                    app.data.(spl).GateTags{1, tg} = kid_path;
                    app.data.(spl).GateTags{2, tg} = kid_name;
                end

                % Update Gate Tags
                if iscell(new_path)
                    new_path = new_path{1};
                end
                if iscell(parent)
                    parent = parent{1};
                end
                if iscell(new)
                    new = new{1};
                end
                
                app.data.(spl).GateTags.(tag){1} = new_path;
                app.data.(spl).GateTags.(tag){2} = strcat(parent, '/', new);
                % why aren't these equivelent?
% % %                 app.data.(spl).GateTags{1, idx_in_tags} = new_path;
% % %                 app.data.(spl).GateTags{2, idx_in_tags} = strcat(parent, '/', new);

                % Update tree
                app.data.(spl).tree = app.data.(spl).tree.rename_kid(new, old_path);
            end
        end
    end
end