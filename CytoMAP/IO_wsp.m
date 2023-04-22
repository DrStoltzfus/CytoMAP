classdef IO_wsp
    methods (Static)
        function Import_Definitions_Func
% % %             import tree.*;
% % %             import Constants.*;
% % %             import Helper.*;
        end

        function func_WSPload(app)
            % FUNC_WSPLOAD Front-End of loading .wsp files function. Allows user to load
            % output of FlowJo along with all of the .fcs files mentioned
            % in that .wsp file.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %
            % Modifies:
            %   - app - Specifically creates app.data.(smpl) where smpl is
            %       name of sample corresponding to this .wsp file.
            
            %% Load .wsp file and data from all .fsj files connected to it.
            % Get filename and path from user
            [fnm, pnm]=uigetfile({'*.wsp', 'FlowJo files (*.wsp)'}, ...
                            'Select FLOWJO workspace file','Multiselect', 'off');
            % Window closed/File not chosen case.
            if (fnm == 0) % || pnm == 0)
                return;
            end
            IO_wsp.func_WSPload_backend(app, pnm, fnm);
        end

        function func_WSPload_backend(app, pnm, fnm)
            % FUNC_WSPLOAD_BACKEND Back-End of loading .wsp files function. Allows user to load
            % output of FlowJo along with all of the .fcs files mentioned
            % in that .wsp file.
            %
            % Input:
            %   - app - Instance of CytoMAP.
            %   - fnm - char, string - Filesname of .wsp which to load into
            %       CytoMAP.
            %   - pnm - char, string - Folder in which this file is.
            %
            % Modifies:
            %   - app - Specifically creates app.data.(smpl) where smpl is
            %       name of sample corresponding to this .wsp file.
            
            vPD = waitbar(0, 'Initializing');
            waitbar(0.5, vPD, 'Loading Data');
            % Get trees, tags, and datasets from .wsp.
            [trees, tables, datas, smpl2group_tab, sample_names] = load_wsp(strcat(pnm, fnm));
            if isempty(datas{end})
                close(vPD)
                return
            end

            waitbar(0.8, vPD, 'Building Trees');
            % Arrange each tuple of tree, tags, and datasets into a CytoMAP
            % friendly data structure.
            for t=1:numel(trees)
                dat_t = struct;

                dat_t.Surfaces = struct;
                dat_t.Surfaces.UDS = struct;
                dat_t.Surfaces.CCN = struct;
                dat_t.Surfaces.RSN = struct;

                dat_t.MetaData = struct;
                dat_t.MetaData.Ver = Constants.CURR_VER;
                dat_t.MetaData.Annotation = 0;
                if ~isempty(smpl2group_tab) && ismember(trees{t}.name, smpl2group_tab.Properties.VariableNames)
                    dat_t.MetaData.Group = smpl2group_tab.(trees{t}.name);
                else
                    dat_t.MetaData.Group = {1};
                end

                dat_t.AllCells = datas{t};
                dat_t.GateTags = tables{t};
                spl_name = sample_names{t};

                % Add a short name row to the GateTags folder
                dat_t.GateTags = [dat_t.GateTags; dat_t.GateTags];
                for gaten=1:numel(dat_t.GateTags(1,:))
                    % Pull the full name
                    FullName = dat_t.GateTags{1,gaten};
                    % Add only the immediate parent and kid name
                    ShortName = strsplit(FullName, '/');
                    if numel(ShortName)==1
                        ShortName = ShortName{1};
                    else
                        ShortName = [ShortName{end-1} '/' ShortName{end}];
                    end
                    ShortName = strrep(strrep(strrep(strrep(ShortName, 'Gate_', ''), '_', ' '), 'NEG', '-'), 'POS', '+');
                    % to asign a string to an element of a table that has a
                    % different length than the string currently occupying
                    % that element you MUST use dot notation?
                    % Pull element designation
                    dat_t.GateTags.(dat_t.GateTags(2,gaten).Properties.VariableNames{1}) = [{FullName}; {ShortName}];
                end
                % Remove Gate_0 if necessary
                if ismember(strcat(Constants.gate_tag, '0'), dat_t.GateTags.Properties.VariableNames)
                    dat_t.GateTags = removevars(dat_t.GateTags, strcat(Constants.gate_tag, '0'));
                end


                dat_t.FNames = fnm;
                dat_t.Path = pnm;
                dat_t.tree = trees{t};
                % Add new data to app.data.(name of that sample)
                app.data.(spl_name) = dat_t;
                if strcmp(app.DataN.Items{end}, 'No Data Loaded')
                    app.DataN.Items{end} = spl_name;
                else
                    app.DataN.Items{end + 1} = spl_name;
                end
            end

            waitbar(1, vPD, 'Done!');
            close(vPD)

            function [experiment_trees, tag_table, data, smpl2group_tab, sample_names] = load_wsp(file_name)
                %{
                    Description:
                        process_wsp takes a path to a .wsp file and returns:
                            - a tree describing it's populations,
                            - a table from tags to paths, relating a tag to a tree node inside of previously returned tree
                            - a data table which reads information from related .fcs files to a given .wsp file.

                    Input Args:
                        - file_name - an (absolute?) path to a .wsp file.
                           *** Requirement: all .fcs files this .wsp file refers to,     ***
                           *** are either in the same directory, or in subdirectory      ***
                           *** specified within the .wsp file                            ***
                           *** i.e. structure of directories out of flowjo has to remain ***

                    Output Args:
                        - experiment_trees - cell of trees containing the
                            structure of samples. Each sample has it's own tree
                            in cell with each gate being represented as a node,
                            containing a gate information, and a tag related to it.
                        - tag_table - cell of tables. Each table, given a tag, points to absolute
                            path in tree, to reach a node that tag relates to.
                        - data - data from the .fcs files this .wsp is pointing to,
                            with all the channels shared in all of .fcs files, as well
                            as gating information given the tree.
                        - Sample Names - Name of a sample. Index-wise
                            corresponding to the experiment tree/tag table
                %}

                % Get the .wsp file and access samples.
                xDoc = xml2struct(file_name);
                groups = xDoc.Children(get_all_nodes(xDoc, 'Groups'));
                sample_list = xDoc.Children(get_all_nodes(xDoc, 'SampleList'));
                clear xDoc
                sample_idxs = get_all_nodes(sample_list(1), 'Sample');

                smpl2group_tab = get_groups(groups, sample_list);

                % Initialize the variables to be returned.
                experiment_trees = cell(size(sample_idxs, 1));
                tag_table = cell(size(sample_idxs, 1));
                data = cell(size(sample_idxs, 1));

                sample_names = cell(numel(sample_idxs),1);
                for sample_index=1:numel(sample_idxs)
                    % Loop over each sample.
                    s_i = sample_idxs(sample_index);

                    % Get it's name.
                    sample = sample_list.Children(s_i);
                    sample_node = sample.Children(get_all_nodes(sample(1), 'SampleNode'));
                    for attr=1:numel(sample_node.Attributes)
                        if strcmp(sample_node.Attributes(attr).Name, 'name')
                            sample_name = strsplit(sample_node.Attributes(attr).Value, '.');
                            sample_name = sample_name{1};
                        end
                    end
                    sample_name = Helper.valid_var(sample_name);
                    if iscell(sample_name)
                        sample_name = sample_name{1};
                    end
                    sample_names{sample_index} = sample_name;

                    % Add subpopulation of the sample.
                    the_pop = sample_node.Children(get_all_nodes(sample_node(1), 'Subpopulations'));
                    param_map = param_to_names(sample);  % Map parameters to actual names: eg. Parameter 1 -> CD3
                    [minimum, maximum] = find_min_max_ranges(sample); % Find minimal and maximal values in the sample.
                    sample_tree = tree('All');  % Create an empty tree node, with the name of sample.
                    sample_tree = sample_tree.add_kid(process_subpopulation(the_pop, NaN));  % Add multiple kids under the sample
                    data{sample_index} = get_fcs_dataset(file_name, sample, sample_name); % Reads a fcs dataset in any of the subdirectories.
                    if isempty(data{sample_index})
                        return;
                    end
                    [experiment_trees{sample_index}, tag_table{sample_index}] = add_gates(sample_tree, table); % Creates a tree with a sample name as a root, and tag table from tags to tree nodes.
                    data{sample_index} = experiment_trees{sample_index}.gate_cells(data{sample_index});
                end

                function smpl_to_group_tab = get_groups(group_handle, sample_list_handle)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    % Input:
                    %   - group_handle - Handle to Group node from wsp/xml
                    %       struct.
                    %
                    %   - sample_list_handle - Handle to Sample List node 
                    %       from wsp/xml struct. 
                    %
                    % Output:
                    %   - smpl_to_group_tab - Table, where columns are
                    %       named after the samples, and rows corresponds
                    %       groups.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    group_nodes = get_all_nodes(group_handle, 'GroupNode');

                    % If there are no group just give back an empty vector.
                    if isempty(group_nodes)
                        smpl_to_group_tab = [];
                        return;
                    end

                    % First find mappings from SampleIDs to Sample Names.
                    smpl_to_group_tab = table;
                    smpl_ids_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
                    smpls = get_all_nodes(sample_list_handle, 'Sample');
                    for smpl_i=1:numel(smpls)
                        smpl = sample_list_handle.Children(smpls(smpl_i));
                        smpl_node = smpl.Children(get_all_nodes(smpl, 'SampleNode'));
                        smpl_id = 0;
                        smpl_name = '';
                        for attr_idx = 1:numel(smpl_node.Attributes)
                            if strcmp(smpl_node.Attributes(attr_idx).Name, 'sampleID')
                                smpl_id = str2double(smpl_node.Attributes(attr_idx).Value);
                            end
                            if strcmp(smpl_node.Attributes(attr_idx).Name, 'name')
                                smpl_name = strsplit(smpl_node.Attributes(attr_idx).Value, '.');
                                smpl_name = Helper.valid_var(smpl_name);
                                smpl_name = smpl_name{1};
                            end
                        end
                        if smpl_id ~= 0 && ~isempty(smpl_name)
                            smpl_ids_map(smpl_id) = smpl_name;
                        end
                    end

                    % Process Groups.
                    for gn_idx = 1:numel(group_nodes)
                        gn_i = group_handle.Children(group_nodes(gn_idx));
                        g_i = gn_i.Children(get_all_nodes(gn_i, 'Group'));
                        g_name = '';
                        for g_a_i=1:numel(g_i.Attributes)
                            if strcmp(g_i.Attributes(g_a_i).Name, 'name')
                                g_name = g_i.Attributes(g_a_i).Value;
                            end
                        end
                        if iscell(g_name)
                            g_name = g_name{1};
                        end
                        sample_refs_root = g_i.Children(get_all_nodes(g_i, 'SampleRefs'));
                        if isempty(sample_refs_root)
                            break;
                        end
                        samples_refs = get_all_nodes(sample_refs_root, 'SampleRef');
                        for s_r_idx = 1:numel(samples_refs)
                            % Extract the name of sample
                            s_r = sample_refs_root.Children(samples_refs(s_r_idx));
                            s_name = smpl_ids_map(str2double(s_r.Attributes.Value));

                            % Process the name
                            s_name = Helper.valid_gate(s_name);
                            if iscell(s_name)
                                s_name = s_name{1};
                            end

                            % If that sample is not in table yet, add it as a column
                            if size(smpl_to_group_tab, 1) == 0
                                smpl_to_group_tab.(s_name) = {g_name};
                                continue;
                            elseif ~ismember(s_name, smpl_to_group_tab.Properties.VariableNames)
                                smpl_to_group_tab.(s_name) = [{g_name}, cell(size(smpl_to_group_tab, 1) - 1)];
                                continue;
                            end
                            % Column exists, table is not fully empty.

                            add_idx = 1 + sum(~cellfun('isempty', smpl_to_group_tab.(s_name)));
                            % If that column is filled append the whole row
                            if add_idx > size(smpl_to_group_tab, 1)
                                % Make a row with all empty values, but one at index of the sample.
                                one_hot_row = cell(1, size(smpl_to_group_tab, 2));
                                one_hot_row{strcmp(smpl_to_group_tab.Properties.VariableNames, s_name)} = g_name;
                                smpl_to_group_tab = [smpl_to_group_tab; one_hot_row];
                            else
                                smpl_to_group_tab(add_idx, s_name) = {g_name};
                            end
                        end
                    end
                end

                % DEPRACATED
                function data = combine_samples(dataset) %#ok<DEFNU>
                    %------------------------------------------------------
                    % Combines a table of datasets (also tables).
                    % A new, combined table contains every cell from all
                    % datasamples (without repetition), and columns which
                    % are in every dataset.
                    %------------------------------------------------------

                    % If there is only one dataset, then just return it.
                    if size(dataset, 1) == 1
                        data = dataset{1};
                    else
                        % Keep track of columns contained in each dataset.
                        valid_vars = dataset{1}.Properties.VariableNames;

                        % For each dataset update valid columns.
                        for i=dataset{2:end}
                            valid_vars = valid_vars(ismember(valid_vars, i.Properties.VariableNames));
                        end

                        % Take only valid columns, from first dataset.
                        data = dataset{1};
                        data = data(:, valid_vars);

                        % For every other dataset, concatanate every cell
                        % that is not yet contained in new dataset.
                        for i=dataset{2:end}
                            j = i(~ismember(i(:, valid_vars), data), valid_vars);
                            data = unique([data; j]);  % Should not be neccessary. Sanity check.
                        end
                    end
                end

                function dataset = get_fcs_dataset(wsp_path, sample, sample_name)
                    % Loads a dataset given a path to local .wsp file,
                    % and a sample specifying path to .fcs file.
                    % Sample can point to path on other device,
                    % but dataset will be loaded successfully,
                    % if relative location of .wsp and .fcs files
                    % have not changed.
                    % Otherwise, it will ask user to point program to path
                    % of the .fcs file.

                    % Get directory of .wsp file.
                    gen_path = strsplit(wsp_path, filesep);
                    gen_path(1:end-1);
                    gen_path = strjoin(gen_path(1:end-1), filesep);
                    dataset = [];

                    % Get path to .fcs file (possibly on other device)/
                    fcs_name = sample.Children(get_all_nodes(sample, 'DataSet'));
                    for key=fcs_name.Attributes
                        if strcmp(key.Name, 'uri')
                            % Get name and path to .fcs file.
                            fcs_name = strsplit(key.Value, '/');
                            fcs_path = fcs_name(1:end-1);
                            fcs_name = fcs_name{end};
                            fcs_name = strrep(fcs_name, '%20', ' ');

                            total_path = strjoin({gen_path fcs_name}, filesep);
                            % While path to .fcs doesn't exist, go deeper, with name specified by
                            % combination of local path to .wsp file and end of path to .fcs file.
                            n = 0;
                            while ~isfile(total_path)
                                n=n+1;
                                try
                                    % Add last subdirectory of .fcs path to currently searched path.
                                    fcs_name = strjoin({fcs_path{end} fcs_name}, filesep);
                                    fcs_name = strrep(fcs_name, '%20', ' ');
                                    fcs_path = fcs_path(1:end-1);

                                    total_path = strjoin({gen_path fcs_name}, filesep);
                                catch
                                    % Pull just the file name from the path
                                    [~, name, ectension] = fileparts(fcs_name);
                                    fcs_name = [name ectension];
                                    if isfile(fcs_name)
                                        total_path =[pwd '\' fcs_name];
                                    else
                                        [~, fcs_path] = uigetfile(...
                                            {'*.fcs', 'FlowJo files (*.fcs)'}, ...
                                            ['Select the .FCS file; ' fcs_name],'Multiselect', 'off'...
                                        );                                           
                                        if isempty(fcs_name) || isempty(fcs_path) || all(fcs_name == 0) || all(fcs_path == 0)
                                            errordlg("Aborting. No .fcs file chosen");
                                            return;
                                        end
                                        cd(fcs_path)
                                        total_path = [fcs_path fcs_name];
                                    end
                                end
                                % Just in case
                                if n>100
                                    break
                                end
                            end
                            break;
                        end
                    end

                    % Make dataset into table, with columns specified by param_map created earlier.
                    values = param_map.values();
                    vals = cell(size(values));
                    for v=1:numel(values)
                        if iscell(values{v})
                            vals(v) = values{v};
                        else
                            vals(v) = values(v);
                        end
                    end
                    tmp = fca_readfcs(total_path);

                    dataset = array2table(tmp, ...
                        'VariableNames', vals);

                    % Make dataset 3-D. Give Z=0, if Z is none existent.
                    if ~ismember('Z', vals)
                        dataset.Z = zeros(size(dataset, 1), 1);
                    end
                    % Reorder columns, so X, Y, Z are always first.
                    dataset = movevars(dataset, 'X', 'Before', 1);
                    dataset = movevars(dataset, 'Y', 'After', 'X');
                    dataset = movevars(dataset, 'Z', 'After', 'Y');

                    % Add Compensated Parameters
                    [coeffs, rows] = get_comp_mat(sample);
                    % If empty then there ain't any
                    if ~isempty(coeffs)
                        % Get columns on which spillover matrix is based
                        spill_help = table2array(dataset(:, rows));
                        for c=1:numel(coeffs.Properties.VariableNames)
                            c_n = coeffs.Properties.VariableNames{c};
                            c_coeff = coeffs.(c_n);

                            % c_coeff = table2array(c_coeff);
                            dataset.(c_n) = spill_help * c_coeff;
                        end
                    end

                    % Add Derived Parameters
                    dataset = insert_derived_channels(sample, dataset, sample_name);
                end

                function [coeffs, spill_params] = get_comp_mat(sample)
                    %------------------------------------------------------
                    % get_comp_mat returns a spillover matrix for a given handle to sample.
                    %
                    % Input:
                    %   - sample - handle from xml2struct to a sample from the flowjo file.
                    %
                    % Output:
                    %   - coeffs - Table where each column indicates for which channel the
                    %               compensation should be applied.
                    %               In example if there is a column 'CD8', it means in
                    %               original dataset, CD8 should be overwritten with
                    %               given compensation.
                    %   - spill_params - cell array of string, which indicates which
                    %               rows correspond to which channels.
                    %               It should allow easier ordering of columns
                    %               in original dataset, when a compensation is applied.
                    %------------------------------------------------------
                    % Init returns.
                    coeffs = table;
                    spill_params = {};

                    % Get Handle to Spillover Matrix, and each individual spill
                    matrix = get_all_nodes(sample, 'transforms:spilloverMatrix');
                    if isempty(matrix)  % No Spillover Matrix
                        return;
                    else
                        matrix = sample.Children(matrix);
                    end
                    spills = get_all_nodes(matrix, 'transforms:spillover');

                    % Loop over spills, and populate table to be returned.
                    for sp_idx=1:numel(spills)
                        % Get handle and name of current column.
                        sp = matrix.Children(spills(sp_idx));
                        name = param2name(sp.Attributes.Value);
                        if iscell(name)
                            name = name{1};
                        end
                        % Make a stub in table
                        coeffs.(name) = zeros(size(coeffs, 1), 1);

                        % Loop over coefficients
                        row = get_all_nodes(sp, 'transforms:coefficient');
                        for r_idx=1:numel(row)
                            r = sp.Children(row(r_idx));

                            % Figure out name and value for each row in the matrix.
                            n_r = r.Attributes(get_attr_name_idx(r, 'data-type:parameter')).Value;
                            n_r = param2name(n_r);

                            v_r = r.Attributes(get_attr_name_idx(r, 'transforms:value')).Value;
                            v_r = str2double(v_r);

                            % Populate spill_params if it's first time looking at rows.
                            if sp_idx == 1
                                if ~iscell(n_r)
                                    n_r = {n_r};
                                end
                                spill_params(end + 1) = n_r;
                            end

                            % Update the coefficients table.
                            in_table_idx = strcmp(n_r, spill_params);
                            coeffs(in_table_idx, name) = {v_r};
                        end

                        % If the spillover is literally just all 0's and one 1, then
                        % there is no spillover, and such columns can be removed coefficients.
                        if sum(coeffs{:, name}) == 1 && coeffs{strcmp(name, spill_params), name} == 1
                            coeffs = removevars(coeffs, name);
                        end
                    end
                end

                function tab = insert_derived_channels(sample, tab, sample_name)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Input:
                    %   - sample - handle from wsp/xml struct to a Sample.
                    %
                    %   - tab - table to which derived parameters should be
                    %       added
                    %
                    % Output:
                    %   - tab - table with derived parameters added to it.
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    dps_root = get_all_nodes(sample, 'DerivedParameters');
                    if isempty(dps_root)
                        return;
                    end
                    dps_root = sample.Children(dps_root);
                    dps = get_all_nodes(dps_root, 'DerivedParameter');
                    for dp_idx=1:numel(dps)
                        dp = dps_root.Children(dps(dp_idx));
                        name = dp.Attributes(get_attr_name_idx(dp, 'name')).Value;
                        name = Helper.valid_channel(name);
                        if iscell(name)
                            name = name{1};
                        end

                        formula = dp.Attributes(get_attr_name_idx(dp, 'formula')).Value;

                        tab.(name) = process_derived_formula(formula, tab, sample_name);
                    end

                end

                function col = process_derived_formula(formula, tab, sample_name)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % process_derived_formula
                    %   takes a formula from .wsp file, which is based on
                    %   columns in tab, or parameters defined in FloJo
                    %   corresponding to such columns. From this information
                    %   returns a column vector this formula corresponds to.
                    %
                    % Input:
                    %     - formula - string which contains mathematical operations
                    %           as defined in MatLab (except for edge cases, look
                    %           at the end of function for specifics), and
                    %           either arguments if form of "tab.column", or
                    %           in form "<Param name= "ParamN">.
                    %     - tab     - Table, in which each column/parameter
                    %           mentioned in formula is present as a column.
                    %
                    % Output:
                    %     - col     - Column vector, which is a result of
                    %           formula defined calculation, with parameters
                    %           specified in table.
                    %           size(col) = [size(tab, 1), 1]
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Make a new formula in terms of actual parameters

                    % Note where those parameters are encoded
                    param_starts = strfind(formula, '<Param name=');
                    param_ends = zeros(param_starts);
                    params = cell(numel(param_starts));

                    % Figure out the ordering of the parameters
                    for p=1:numel(param_starts)
                        % See where parameter xml ends
                        param_end = strfind(formula(param_starts(p):end), "/>");
                        param_end = param_starts(p) + param_end(1);
                        param_ends(p) = param_end;

                        % Get only name, and convert it to valid name.
                        param = split(formula(param_starts(p):param_end), '"');
                        param = param2name(param{2});
                        params(p) = cellstr(param);
                    end

                    % Substitute those ranges with references to specific columns
                    % in given table tab
                    if numel(params) > 1
                        new_formula = formula(1:param_starts(1) - 1);
                        for p=1:numel(param_starts) - 1
                            new_formula = strcat(new_formula, 'tab.', params{p});
                            new_formula = strcat(new_formula, ...
                                formula(param_ends(p) + 1 : param_starts(p + 1) - 1));
                        end
                        new_formula = strcat(new_formula, 'tab.', params{end}, ...
                            formula(param_ends(end) + 1:end));
                    elseif numel(params) == 1
                        new_formula = formula(1:param_starts(1) - 1);
                        new_formula = strcat(new_formula, 'tab.', params{1});
                        new_formula = strcat(new_formula, formula(param_ends(1) + 1:end));
                    else
                        new_formula = formula;
                    end

                    % Substitute differences in notation in MatLab
                    new_formula = strrep(new_formula, "/", " ./ ");
                    new_formula = strrep(new_formula, "*", " .* ");

                    %% Evaluate to a single column
                    
                    try
                        col = eval(new_formula);
                    catch
                        str = ['Error Loading Derived Parameter: ';  ...
                         params(1:end,1); ...
                        'Formula:'; ...
                        formula; ...
                        'Sample'; ...
                        sample_name; ...
                        'PARAMETER NOT IMPORTED!'];
                        msgbox(strcat(str))
                        col = zeros(size(tab, 1), 1);
                    end
                end

                function idxs = get_attr_name_idx(handle, name)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Input:
                    %   - handle - struct which is either 
                    %       parameter or transform from the wsp/xml struct.
                    %       NOTE: That handle has to contain 'Attributes' field.
                    %   - name - name of the attribute to search for.
                    %
                    % Output:
                    %   - idxs - indexes of the attribute with given name
                    %       that occurs under attributes field of the handle.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    idxs = [];
                    if ~isstruct(handle)
                        return;
                    end
                    if ~isfield(handle, 'Attributes')
                        return;
                    end

                    for i=1:numel(handle.Attributes)
                        if ~isstruct(handle.Attributes(i))
                        end
                        if ~isfield(handle.Attributes(i), 'Name')
                        end

                        if strcmp(name, handle.Attributes(i).Name)
                            idxs(end + 1) = i;
                        end
                    end
                end

                function [minimum, maximum] = find_min_max_ranges(sample)
                    % Finds minimum and maximum ranges of all parameters in
                    % the given sample (wsp/xml struct).
                    % Useful for loading gates such as Spider Gate etc.
                    transforms = sample.Children(get_all_nodes(sample, 'Transformations'));
                    indexes = get_all_nodes(transforms, 'transforms:linear');
                    minimum = realmax;
                    maximum = -realmax;
                    for kid_ind=indexes
                        if maximum < str2double(transforms.Children(kid_ind).Attributes(2).Value)
                            maximum = str2double(transforms.Children(kid_ind).Attributes(2).Value);
                        end

                        if minimum > str2double(transforms.Children(kid_ind).Attributes(3).Value)
                            minimum = str2double(transforms.Children(kid_ind).Attributes(3).Value);
                        end
                    end
                end

                function population_list = process_subpopulation(root, parent)
                    % Given a pointer to Subpopulation of .wsp/xml file,
                    % and a parent tree node,
                    % returns a tree, which contains all of it's children.
                    pop_idxs = get_all_nodes(root, 'Population');
                    population_list = tree.empty(length(pop_idxs), 0);
                    for pop_idx=1:length(pop_idxs)
                        if ~strcmp(get_name(root.Children(pop_idxs(pop_idx))), 'Histogram Gate Range')
                            population_list(pop_idx) = process_population(root.Children(pop_idxs(pop_idx)), parent);
                        end
                    end
                end

                function population_tree = process_population(root, parent)
                    % Given a pointer to Population of .wsp/xml file, 
                    % and a parent tree node,
                    % returns a tree, which contains all of it's children.
                    population_tree = tree(get_name(root));
                    population_tree = set_gates(root, population_tree, parent);
                    sub_pop = get_all_nodes(root, 'Subpopulations');
                    if ~isempty(sub_pop)
                        for i=process_subpopulation(root.Children(sub_pop), population_tree)
                            population_tree = population_tree.add_kid(i);
                        end
                    end
                end

                function population_tree = get_populations(root, path)
                    % Given pointer to Subpopulation of .wsp/xml file struct
                    % (root), and path to the top of tree, returns a tree
                    % with all populations in that subpopulation loaded.
                    if ~exist('path', 'var')
                        path = '';
                    end

                    population_tree = tree;
                    for pop_node=get_all_nodes(root, 'Population')
                        sub_root = root.Children(pop_node);
                        sub_name = get_name(sub_root);
                        sub_pop = get_all_nodes(sub_root, 'Subpopulations');
                        if ~isempty(sub_pop)
                            sub_tree = get_populations(sub_root.Children(sub_pop), strcat(path,'_',sub_name));
                            population_tree.add_kid(sub_tree);
                        else
                            population_tree.(sub_name) = table;
                        end
                    end
                end

                function name = get_name(root)
                    % Given pointer to Population of .wsp/xml file struct
                    % (root), returns it's name.
                    if isstruct(root)
                        if isstruct(root.Attributes)
                            for attr_idx=1:numel(root.Attributes)
                                if strcmp('name', root.Attributes(attr_idx).Name)
                                    name = root.Attributes(attr_idx).Value;
                                end
                            end
                        end
                    end
                end

                function out = get_all_nodes(root, node_name)
                    % Returns indexes to all Children nodes in root,
                    % with tag specified by name.
                    out = zeros(size(root.Children, 1));

                    for i=1:numel(root.Children)
                        if strcmp(node_name, root.Children(i).Name)
                            out(i) = i;
                        end
                    end
                    out = out(out ~= 0);
                end

                function tree = set_gates(population, tree, parent)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Given a population
                    % (pointer to Population in .wsp/xml struct),
                    % sets gates in corresponding tree node (tree).
                    % Since some gates are specified in terms of parent
                    % gates, the parent (tree) is also required.
                    % Returns a modified tree, with gates set.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if isstruct(population)
                        gate = get_all_nodes(population, 'Gate');
                        if ~isempty(gate)
                            gate = population.Children(gate); % There's only one Gate field in .wsp
                            graph = population.Children(get_all_nodes(population, 'Graph'));

                            for i=1:numel(gate.Children)
                                if ~strcmp(gate.Children(i).Name, '#text')
                                    switch gate.Children(i).Name
                                        case 'gating:EllipsoidGate'
                                            % DEBUG STUFF UNDER. UNCOMMENT
                                            if Constants.ellipse_debug
                                                gate_type = 'all';
                                                pts = [];
                                                axes = [];
                                            else
                                                errordlg(strcat('Due to errors in FloJo we do not support Ellipse Gates.', newline, ...
                                                'Please.'));
                                                error('Ellipses ain''t supported');

                                                gate_type = 'ellipsoid';
                                                axes = get_axes(gate.Children(i), graph);
                                                [pts, axes] = ellipsoid_ranges(gate.Children(i), axes);
                                            end
                                        case 'gating:RectangleGate'
                                            gate_type = 'rectangle';
                                            axes = get_axes(gate.Children(i), graph);
                                            pts = rectangle_ranges(gate.Children(i), parent, tree.name, axes);
                                        case 'gating:PolygonGate'
                                            gate_type = 'polygon';
                                            axes = get_axes(gate.Children(i), graph);
                                            pts = polygon_ranges(gate.Children(i));
                                        otherwise
                                            disp("Not supported Gate Type");  % For debug purposes.
                                            disp(gate.Children(i).Name)  % Show what is name of this new type
                                            continue  % If it's not a gate, ignore
                                    end
                                    tree.gate_points = pts;
                                    tree.gate_axes = axes;
                                    tree.gate_type = gate_type;
                                    break % There is only one gate per population
                                end
                            end
                        end
                    end

                    function vertices = get_vertices(gate)
                        % Returns verticies of a gate. i.e. points around
                        % which gate is specified.
                        indexes = get_all_nodes(gate, 'gating:vertex');
                        vertices = zeros(size(indexes, 1), 2);
                        for k=1:numel(indexes)
                            vars = get_all_nodes(gate.Children(indexes(k)), 'gating:coordinate');
                            for j=1:numel(vars)
                                vertices(k, j) = str2double(gate.Children(indexes(k)).Children(vars(j)).Attributes.Value);
                            end
                        end
                    end

                    function axes = get_axes(gate, graph)
                        % Returns axes of a gate.
                        % These are dimensions/parameters for which gate is
                        % specified.
                        %
                        % If there is only 1 dimension specified in gate,
                        % the get_axes looks at corresponding graph for
                        % additional information.
                        dimensions = get_all_nodes(gate, 'gating:dimension');
                        if length(dimensions) == 1
                            dimensions = get_all_nodes(graph, 'Axis');
                            axes = string.empty([length(dimensions) 0]);
                            for j=1:numel(dimensions)
                                axes(j) = param2name(graph.Children(dimensions(j)).Attributes(4).Value);
                            end
                        else
                            axes = string.empty([length(dimensions) 0]);
                            for j =1:numel(dimensions)
                                axes(j) = param2name(gate.Children(dimensions(j)).Children(2).Attributes.Value);
                            end
                        end
                    end

                    function vertices = rectangle_ranges(gate, parent, name, axes)
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Reads vertices for rectangular region.
                        % Also processes spider gates (assumes straight
                        % arms), and histogram gates from FlowJo
                        % (ones which bisect into 2 regions).
                        %
                        %
                        % Input Args:
                        %   - gate - Handle to struct made from xml, 
                        %       pointed at the gate.
                        %   - parent - Parent node in the tree (already
                        %       processed with gate points and axes)
                        %   - name - Name of the gate. It is used to
                        %       identify type of gate.
                        %   - axes - Axes on passed in gate is defined.
                        %
                        % Output Args:
                        %   - vertices - 4x2 array with, pair of numbers in
                        %       each row. Coordinates correspond to values
                        %       of ranges populated based on types of
                        %       gates (can go to minimal value found in the
                        %       sample if unbounded type etc.)
                        %       Vertices are in order:
                        %       (left-bottom;
                        %        left-top;
                        %        right-top;
                        %        right-bottom)
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        dimensions = get_all_nodes(gate, 'gating:dimension');
                        ranges = zeros(2, 2);
                        if length(dimensions) == 1
                            if startsWith(name, 'Right_Bisect')
                                axis_name = param2name(gate.Children(dimensions).Children(2).Attributes.Value);
                                other_child_axis = axes(axes ~= axis_name);

                                new_max = str2double(gate.Children(dimensions).Attributes(1).Value);
                                defined_axis = parent.gate_axes == axis_name;
                                other_parent_axis = parent.gate_axes(~defined_axis);
                                if defined_axis(1) == 1
                                    ranges(1, 1) = min(parent.gate_points(:, 1));
                                    ranges(1, 2) = new_max;
                                    if strcmp(other_child_axis, other_parent_axis)
                                        ranges(2, 1) = min(parent.gate_points(:, 2));
                                        ranges(2, 2) = max(parent.gate_points(:, 2));
                                    else
                                        ranges(2, 1) = minimum;
                                        ranges(2, 2) = maximum;
                                    end
                                else
                                    ranges(2, 1) = min(parent.gate_points(:, 2));
                                    ranges(2, 2) = new_max;
                                    if strcmp(other_child_axis, other_parent_axis)
                                        ranges(1, 1) = min(parent.gate_points(:, 1));
                                        ranges(1, 2) = max(parent.gate_points(:, 1));
                                    else
                                        ranges(1, 1) = minimum;
                                        ranges(1, 2) = maximum;
                                    end
                                end
                            elseif startsWith(name, 'Left_Bisect')
                                axis_name = param2name(gate.Children(dimensions).Children(2).Attributes.Value);
                                other_child_axis = axes(axes ~= axis_name);

                                new_min = str2double(gate.Children(dimensions).Attributes(1).Value);
                                defined_axis = parent.gate_axes == axis_name;
                                other_parent_axis = parent.gate_axes(~defined_axis);
                                if defined_axis(1) == 1
                                    ranges(1, 1) = new_min;
                                    ranges(1, 2) = max(parent.gate_points(:, 1));
                                    if strcmp(other_child_axis, other_parent_axis)
                                        ranges(2, 1) = min(parent.gate_points(:, 2));
                                        ranges(2, 2) = max(parent.gate_points(:, 2));
                                    else
                                        ranges(2, 1) = minimum;
                                        ranges(2, 2) = maximum;
                                    end
                                else
                                    ranges(2, 1) = new_min;
                                    ranges(2, 2) = max(parent.gate_points(:, 2));
                                    if strcmp(other_child_axis, other_parent_axis)
                                        ranges(1, 1) = min(parent.gate_points(:, 1));
                                        ranges(1, 2) = max(parent.gate_points(:, 1));
                                    else
                                        ranges(1, 1) = minimum;
                                        ranges(1, 2) = maximum;
                                    end
                                end
                            end
                        else
                            if startsWith(name, ["Q1_", "Q2_", "Q3_", "Q4_"])
                                % Quandrant Spider gates
                                for j=1:numel(dimensions)
                                    if strcmp(gate.Children(dimensions(j)).Attributes.Name, 'gating:max')
                                        ranges(j, 1) = minimum;
                                        ranges(j, 2) = str2double(gate.Children(dimensions(j)).Attributes.Value);
                                    else
                                        ranges(j, 1) = str2double(gate.Children(dimensions(j)).Attributes.Value);
                                        ranges(j, 2) = maximum;
                                    end
                                end

                            else
                                % Rectangular gates
                                for j=1:numel(dimensions)
                                    if length(gate.Children(dimensions(j)).Attributes) == 1
                                        if strcmp(gate.Children(dimensions(j)).Attributes.Name, 'gating:max')
                                            ranges(j, 1) = minimum;
                                            ranges(j, 2) = str2double(gate.Children(dimensions(j)).Attributes.Value);
                                        else
                                            ranges(j, 1) = str2double(gate.Children(dimensions(j)).Attributes.Value);
                                            ranges(j, 2) = maximum;
                                        end
                                    else
                                        ranges(j, :) = [str2double(gate.Children(dimensions(j)).Attributes(2).Value) ...
                                            str2double(gate.Children(dimensions(j)).Attributes(1).Value)];
                                    end
                                end
                            end
                        end
                        % Populate vertices
                        vertices = [ranges(1, 1) ranges(2, 1); ...
                            ranges(1, 1) ranges(2, 2); ...
                            ranges(1, 2) ranges(2, 2); ...
                            ranges(1, 2) ranges(2, 1)];
                    end

                    function [params, axes] = ellipsoid_ranges(gate, axes)
                        % -------------------------------------------------
                        % WARNING BUG! DOESN'T WORK.
                        % Reason is that FlowJo processes ellipsoids
                        % incorrecly, and there is no possible way of
                        % retrieving information correctly, which is a bug
                        % on side of FlowJo.
                        %
                        % Given gate and it's axes, returns params and new axes (can be switched places), which are specifiying the gate.
                        % Tries to solve an equation for ellipse, and returns parameter based on them.
                        % Params is defined as follows:
                        %   params(1) = a, params(2) = b (Major and minor axes)
                        %   params(3) = alpha (Tilt)
                        %   params(4) = x0, params(5) = y0 (Center of ellipsoid)
                        % -------------------------------------------------
                        error("It ain't working. Switch ellipse gates to polygon gates and restart." + newline + ...
                              "You should also check if the ellipse gates are the right ones.");
                        tmp_foci = get_all_nodes(gate, 'gating:foci');
                        tmp_foci = gate.Children(tmp_foci);
                        tmp_edges = get_all_nodes(gate, 'gating:edge');
                        tmp_edges = gate.Children(tmp_edges);

                        focis = get_vertices(tmp_foci);
                        edges = get_vertices(tmp_edges);

                        % Find the center point of the ellipse
                        center = [min(edges(:, 1))+(max(edges(:, 1))-min(edges(:, 1)))/2, ...
                                  min(edges(:, 2))+(max(edges(:, 2))-min(edges(:, 2)))/2];

                        % Calculate the magnitude of the semi-major axis
                        a = max(sqrt((edges(:, 1)-center(1)).^2 + (edges(:, 2)-center(2)).^2));

                        alpha = atan((focis(2, 2) - focis(1, 2)) / ...
                            (focis(2, 1) - focis(1, 1)));

                        % If major axis is oriented in y direction, flip everything
                        % (including axes), and recalculate alpha
                        if abs(alpha) > pi/4
                            axes = fliplr(axes);
                            focis = fliplr(focis);
                            edges = fliplr(edges);
                            alpha = atan((focis(2, 2) - focis(1, 2)) / ...
                                (focis(2, 1) - focis(1, 1)));
                        end

                        options = optimset('Display', 'off');

                        center = focis(1, :) + focis(2, :);
                        center = center ./ 2;

                        cos_alpha = cos(alpha);
                        sin_alpha = sin(alpha);

                        ecc = pdist([center; focis(1, :)]);
                        f = @ellipsoid_system;
                        starting_point = [1 1];
                        params = fsolve(f, starting_point, options);

                        params(end+1) = alpha;
                        params(end+1) = center(1);
                        params(end+1) = center(2);

                        function F = ellipsoid_system(c)
                            F = zeros(size(edges, 1), 1);

                            for idx=1:size(edges, 1)
                                x_diff = center(1) - edges(idx, 1);
                                y_diff = center(2) - edges(idx, 2);
                                F(idx) = (x_diff*cos_alpha + y_diff*sin_alpha)^2/(c(1)^2) + ...
                                (x_diff*sin_alpha - y_diff*cos_alpha)^2/(c(2)^2) - 1;
                            end

                            F(end+1) = c(1)^2 - c(2)^2 + ecc^2;
                        end
                    end

                    function vertices = polygon_ranges(gate)
                        % Returns points specifying polygon gate.
                        vertices = get_vertices(gate);
                    end

                end

                function map = param_to_names(sample)
                    % Given sample maps FlowJo parameter names to their 
                    % actual (valid MatLab fieldnames) equivalents.
                    % Example:
                    %   - 'Parameter_11' is mapped to 'X',
                    %   - 'Parameter_9' is mapped to 'Mean_CD11c' etc.
                    keywords = sample.Children(get_all_nodes(sample, 'Keywords'));
                    indexes = get_all_nodes(keywords, 'Keyword');

                    map = containers.Map('KeyType', 'double', 'ValueType', 'any');

                    % initialize empty channel name
                    chnl = [];
                    % Flip this if it is in the wrong order (later versions of FlowJo flipped this)
                    if ~strcmp(keywords.Children(2).Attributes(1).Value(1:2), '$P')
                        keywords.Children = flip(keywords.Children,2);
                    end

                    for j=indexes
                        if strcmp(keywords.Children(j).Attributes(1).Value(1:2), '$P') && ...
                            ~strcmp(keywords.Children(j).Attributes(1).Value, '$PAR')
                            switch keywords.Children(j).Attributes(1).Value(end)
                                case 'N'
                                    chnl = keywords.Children(j).Attributes(2).Value;
                                case 'S'
                                    % Ignore Comp-Parameter's duplicate names
                                    if ~startsWith(chnl, 'Comp-Parameter_')
                                        key = str2double(keywords.Children(j).Attributes(1).Value(3:end-1));
                                        value = Helper.valid_channel(keywords.Children(j).Attributes(2).Value);
% % %                                         {key; value} 
                                        map(key) = value;
                                    end
                                otherwise
                                    continue
                            end
                        end
                    end
                    
                end

                function name = param2name(param)
                    % Given FlowJo parameter (like 'Parameter_11'), returns it's actual name.
                    % Takes care also of Compensated Parameters
                    % Note: Derived parameters are processed by FlowJo like
                    %   normal parameters. However if something related to 
                    %   derived parameters breaks, this function is a good
                    %   place to double-check.

                    % If empty return X by default.
                    if strcmp(param, '')
                        name = 'X';
                        return
                    end
                    if startsWith(param, 'Pos ') && strcmp(param(end-4:end), '(inv)')
                        % Something weird 'Pos X (inv)' exists
                        name = param(5);
                        return
                    % Some parameters have been scaled in FloJo, These are
                    % called Comp-Parameters
                    
                    elseif startsWith(param, 'Comp-Parameter_')
                        % Pull which parameter Number this is refering to
                        param = str2double(param(16:end));
                        %%%%%%% This parameter needs to be scaled!
                    elseif startsWith(param, 'Parameter_')
                        param = str2double(param(11:end));
                    elseif endsWith(param, 'Derived')
                        name = Helper.valid_channel(param);
                        return
                    end
                    
                    % try to pull the name from the parameter map
                    try
                        name = param_map(param); 
                    catch
                        % In some version FlowJo saves the name but doesn't
                        % add it to the list of parameter names
                        name = Helper.valid_channel(param);
                    end
                end

                function [processed_tree, tmp_app] = add_gates(root, dataset, init_path)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % add_gates creates Gate Tags table, and a tree with
                    % tags based on untagged tree, dataset and path within
                    % the tree given.
                    %
                    % Input Args:
                    %   - root - root of the tree to be tagged. This tree
                    %       should not contain tags
                    %   - dataset - Table (can be empty), to be populated
                    %       with references corresponding to tree nodes.
                    %   - init_path - optional - Path to start tagging at.
                    %       Note: THIS FEATURE IS NOT WELL TESTED/CHECKED
                    %       (since in case of loading .wsp you always start
                    %       at root).
                    %
                    % Output Args:
                    %   - processed_tree - An exact copy of root, but with
                    %       tag field filled at each node.
                    %   - tmp_app - Table with only 1 row, where name of 
                    %       each column is a tag, that was also created in
                    %       the tree, and the cell underneath is path to
                    %       node corresponding to that tag.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Adds gate tags, to the dataset, and to the tree, and creates maps between them.
                    if ~exist('init_path', 'var')
                        init_path = '';
                    end
                    tmp_app = dataset;
                    index = size(dataset, 1);

                    [processed_tree, tmp_app] = recurse_tree(tmp_app, root, Helper.valid_var(init_path));

                    function [root2, app] = recurse_tree(app, root, path)
                        root_tag = strcat([Constants.gate_tag, int2str(index)]);
                        if isempty(path)
                            new_path = root.name;
                        else
                            new_path = strcat(path, '/', root.name);
                        end

                        root.tag = root_tag;
                        app.(root_tag) = new_path;
                        index = index + 1;
                        for sub_kid=1:numel(root.kids)
                            [root.kids.(sub_kid), app] = recurse_tree(app, root.kids.(sub_kid), new_path);
                        end
                        root2 = root;
                    end
                end
            end
        end
    end
end