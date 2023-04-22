classdef Helper

% Helper defines miscelanious functions for loading, converting, preparing, and analyzing data
%
%#ok<*AGROW>

    methods (Static)
        function Import_Definitions_Func
% % %             import tree.*;
% % %             import Constants.*;
        end

        function icons = Import_Icons()
            % IMPORT_ICONS This function imports the icons in the Icons
            % folder 
            
            % pull the path for the current function
            cpath = fileparts(mfilename('fullpath'));

            % Get the path of the icons folder
            folderpath = fullfile([cpath filesep 'Icons']);
            FilePathInfo = dir(fullfile([cpath filesep 'Icons'], '*.bmp'));
            % remove any directory elements without dates
            FilePathInfo = FilePathInfo(~cellfun('isempty', {FilePathInfo.date}));
            % remove any folders 
            FilePathInfo = FilePathInfo(~[FilePathInfo.isdir]);
            % remove weird stuff
            FilePathInfo = FilePathInfo(~strcmp('.', {FilePathInfo.name}));
            FilePathInfo = FilePathInfo(~strcmp('..', {FilePathInfo.name}));
            FilePathInfo = {FilePathInfo.name}';
            icons = struct;
            for i=1:numel(FilePathInfo)
                % Read an image
                [img,~] = imread([folderpath filesep FilePathInfo{i}]);
                % Pull the icon name from the filename
                FilePathInfo{i} = strrep(FilePathInfo{i}, ' ', '_');
                FilePathInfo{i} = strrep(FilePathInfo{i}, '-', '_');
                % Convert the image to double
                img = double(img./max(max(max(img))));
                % Put the image in a structure of icons
                icons.(FilePathInfo{i}(1:end-4)) = img;
            end
            
        end

        function [to_drop, new_vars] = remove_other_duplicates(vars)
            % REMOVE_OTHER_DUPLICATES Compatibility with older version. Changes tags to other_tag
            % if needed.
            % If such change would result in duplicate column, just drops
            % the tag.
        
            to_drop_ch = ismember(vars, strcat(Constants.channel_tag, Constants.other_names));
            to_drop_gate = ismember(vars, strcat(Constants.gate_tag, Constants.other_names));

            ch_len = length(Constants.channel_tag);
            for i=1:numel(to_drop_ch)
                if to_drop_ch(i) == 1
                    val = vars(i);
                    var_wo_tag = {val{1}(ch_len + 1:end)};

                    % If the channel without Ch doesn't exist in dataset, do not remove it, and change the tag
                    if ~ismember(var_wo_tag, vars) && ...
                            ~ismember(strcat(Constants.other_tag, var_wo_tag), vars)
                        vars(i) = strcat(Constants.other_tag, var_wo_tag);
                        to_drop_ch(i) = 0;
                    end
                end
            end

            to_drop = to_drop_ch | to_drop_gate;
            new_vars = vars(~to_drop);
        end

        function to_keep = remove_duplicates(name_cell)
            % REMOVE_DUPLICATES Looks over names and returns indeces to those which are not
            % duplicates (left-most indices)
            % Example for input {'a', 'b', 'c', 'b', 'd', 'a'} it would
            % return bool array [ 1 ,  1 ,  1 ,  0 ,  1 ,  0 ]
            
            if numel(name_cell) <= 1
                to_keep = TRUE;
                return;
            end
            to_keep = true(size(name_cell));
            for n=1:(numel(name_cell)-1)
                if to_keep(n) % If was already tagged as duplicate then don't check it.
                    to_keep(1:n-1) = to_keep(1:n-1) & ~ismember(name_cell(1:n-1), name_cell(n));
                    to_keep(n+1:end) = to_keep(n+1:end) & ~ismember(name_cell(n+1:end), name_cell(n));
                end
            end
        end

        function scan = get_scan(app, dataset)
            % GET_SCAN Given dataset returns it's most recent neighborhood scan.
            % Nothing more that attempt to limit number of things to write.
            
            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end

            scan = app.data.(dataset).(app.data.(dataset).ScanType);
        end

        function [gates, short_gates] = get_gates(app, dataset)
            % GET_GATES Gives easy access to common gates across multiple
            % datasets
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - dataset - cell, char, string - single or multiple names
            %       of datasets to get gates from
            %
            % Output:
            %   - gates - Path to gates which are common for all given
            %       datasets. Can be used to access tree and compare gates
            %       with other samples/datasets.
            %   - short_gates - Short names of gates which are common for
            %       all given datasets. Can be used to display information
            %       to user.
            
            if ~exist('dataset', 'var') || isempty(dataset)  % Use most recent dataset
                dataset = app.DataN.Value;
            end

            if ~iscell(dataset)
                dataset = {dataset};
            end

            % Check for non-existent datasets. If all of them all invalid, quit.
            tmp_dataset = {};
            for Di=1:numel(dataset)
                if isfield(app.data, dataset{Di})
                    tmp_dataset{end + 1} = dataset{Di};
                end
            end

            if isempty(tmp_dataset)
                gates = [];
                short_gates = [];
                return;
            else
                dataset = tmp_dataset;
            end

            % Somewhere something is making the elements of GateTags
            % cells...
            for Di = 1:numel(dataset)
                gates = app.data.(dataset{Di}).GateTags{1, :};
                for Gi =1:numel(gates)
                    if iscell(app.data.(dataset{Di}).GateTags{1, :}{Gi})
                        app.data.(dataset{Di}).GateTags{1, Gi} = app.data.(dataset{Di}).GateTags{1, :}{Gi};
                    end
                end
            end

            % I changed this so it would work with different Long gate
            % names, otherwise ismember returns all zeros
            gates = app.data.(dataset{1}).GateTags{1, :};
            short_gates = app.data.(dataset{1}).GateTags{2, :};

            for i = 2:numel(dataset)
                gates_i = app.data.(dataset{i}).GateTags{1, :};
                gates = gates(ismember(gates, gates_i));

                short_gates_i = app.data.(dataset{i}).GateTags{2, :};
                short_gates = short_gates(ismember(short_gates, short_gates_i));
            end
        end

        function tags = get_gate_tags(app, dataset)
            % GET_GATE_TAGS Gives easy access to tags of common gates 
            % across multiple datasets
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - dataset - cell, char, string - single or multiple names
            %       of datasets to get gate tags from
            %
            % Output:
            %   - tags - Tags of gates which have common paths for all 
            %       given datasets. Can be used to access AllCells, MFIs
            %       etc.
            
            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end

            if ~iscell(dataset)
                dataset = {dataset};
            end

            % Check for non-existent datasets. If all of them all invalid, quit.
            tmp_dataset = {};
            for Di=1:numel(dataset)
                if isfield(app.data, dataset{Di})
                    tmp_dataset{end + 1} = dataset{Di};
                end
            end

            if isempty(tmp_dataset)
                tags = [];
                return;
            else
                dataset = tmp_dataset;
            end

            % Then use an intersection of paths from all datasets to define tags
            tags = app.data.(dataset{1}).GateTags.Properties.VariableNames;
            paths = app.data.(dataset{1}).GateTags{1, :};

            for i = 2:numel(dataset)
                paths_i = app.data.(dataset{i}).GateTags{1, :};
                idxs = ismember(paths, paths_i);
                if isempty(idxs)
                    tags = {};
                    return;
                else
                    paths = paths(idxs);
                    tags = tags(idxs);
                end
            end
            tags = tags(~strcmp(strcat(Constants.gate_tag, '0'), tags));
        end

        function channels = get_channels(app, dataset)
            % GET_CHANNELS Gives easy access to common channels across
            % multiple datasets
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - dataset - cell, char, string - single or multiple names
            %       of datasets to get channels from
            %
            % Output:
            %   - channels - Channels which are common all given datasets.
            %       Can be used to access AllCells, MFIs etc.
            
            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end

            % Then dataset an intersection of channels from both datasets is taken
            if iscell(dataset)
                channels = app.data.(dataset{1}).AllCells.Properties.VariableNames;
                choices = strncmpi(channels, Constants.channel_tag, length(Constants.channel_tag));
                choices = choices | ismember(channels, {'X', 'Y', 'Z'});
                channels = channels(choices);

                % Intersect doesn't like strings
                for c = 1:numel(channels)
                    if isstring(channels{c})
                        channels{c} = char(channels{c});
                    end
                end

                for i = 2:numel(dataset)
                    channels_i = app.data.(dataset{i}).AllCells.Properties.VariableNames;
                    choices_i = strncmpi(channels_i, Constants.channel_tag, length(Constants.channel_tag));
                    choices_i = choices_i | ismember(channels_i, {'X', 'Y', 'Z'});
                    channels_i = channels_i(choices_i);

                    % Intersect doesn't like strings
                    for c = 1:numel(channels_i)
                        if isstring(channels_i{c})
                            channels_i{c} = char(channels_i{c});
                        end
                    end
                    channels = intersect(channels, channels_i, 'stable');
                end

            else
                channels = app.data.(dataset).AllCells.Properties.VariableNames;
                choices = strncmpi(channels, Constants.channel_tag, length(Constants.channel_tag));
                choices = choices | ismember(channels, {'X', 'Y', 'Z'});
                channels = channels(choices);
            end
        end

        function others = get_others(app, dataset, varargin)
            % GET_OTHERS Gives easy access to all common `others` across
            % multiple datasets
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - dataset - cell, char, string - single or multiple names
            %       of datasets to get others from
            %
            % Keyword Input (optional)
            %   - mfi - cell, char, string - default: {} - If mfi is given,
            %       then others are looke for in every mfi specified.
            %       In example if MFIRSN is given, then all returned others
            %       would have to be in MFIRSN in each sample it exists.
            %       If mfi is empty (default), then get_others looks only
            %       at AllCells.
            %       If some samples do not have mfi specified than those
            %       are ignored.
            %
            % Output:
            %   - others - `Others` which are common for all given datasets
            %       Can be used to access AllCells, MFIs etc.
            
            defaults.mfi = '';
            if mod(length(varargin), 2) ~= 0
                error('Uneven number of input arguments! populate_table needs propertyName/propertyValue pairs, after app argument')
            end

            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end

            if ~iscell(dataset)
                dataset = {dataset};
            end

            if ~iscell(defaults.mfi)
                defaults.mfi = {defaults.mfi};
            end

            for pair = reshape(varargin, 2, [])
                if isfield(defaults, lower(pair{1}))
                    defaults.(lower(pair{1})) = pair{2}; % Change default axis value to the given value
                else  % Key is not in the defaults
                    error('%s is not a recognized parameter name: ', pair{1})
                end
            end

            % Check for non-existent datasets. If all of them all invalid, quit.
            tmp_dataset = {};
            for Di=1:numel(dataset)
                if isfield(app.data, dataset{Di})
                    tmp_dataset{end + 1} = dataset{Di};
                end
            end

            if isempty(tmp_dataset)
                others = [];
                return;
            else
                dataset = tmp_dataset;
            end


            others_all = app.data.(dataset{1}).AllCells.Properties.VariableNames;
            choices_all = strncmpi(others_all, Constants.other_tag, length(Constants.other_tag));
            others = others_all(choices_all);
            if ismember({'MFIRSN'}, defaults.mfi) && isfield(app.data.(dataset{1}), 'MFIRSN')
                others_mfi = app.data.(dataset{1}).MFIRSN.Properties.VariableNames;
                choices_mfi = strncmpi(others_mfi, Constants.other_tag, length(Constants.other_tag));
                others = union(others, others_mfi(choices_mfi),'stable');
            end
            if ismember({'MFICCN'}, defaults.mfi) && isfield(app.data.(dataset{1}), 'MFICCN')
                others_mfi = app.data.(dataset{1}).MFICCN.Properties.VariableNames;
                choices_mfi = strncmpi(others_mfi, Constants.other_tag, length(Constants.other_tag));
                others = union(others, others_mfi(choices_mfi),'stable');
            end
            for i = 2:numel(dataset)
                others_all = app.data.(dataset{i}).AllCells.Properties.VariableNames;
                choices_all = strncmpi(others_all, Constants.other_tag, length(Constants.other_tag));
                others_i = others_all(choices_all);
                if ismember({'MFIRSN'}, defaults.mfi) && isfield(app.data.(dataset{i}), 'MFIRSN')
                    others_mfi = app.data.(dataset{i}).MFIRSN.Properties.VariableNames;
                    choices_mfi = strncmpi(others_mfi, Constants.other_tag, length(Constants.other_tag));
                    others_i = union(others_i, others_mfi(choices_mfi),'stable');
                end
                if ismember({'MFICCN'}, defaults.mfi) && isfield(app.data.(dataset{i}), 'MFICCN')
                    others_mfi = app.data.(dataset{i}).MFICCN.Properties.VariableNames;
                    choices_mfi = strncmpi(others_mfi, Constants.other_tag, length(Constants.other_tag));
                    others_i = union(others_i, others_mfi(choices_mfi),'stable');
                end

                others = intersect(others, others_i, 'stable');
            end
        end

        function neighs = get_neighs(app, dataset, mfi)
            % GET_NEIGHS Gives easy access to common neighborhoods across
            % multiple datasets
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - dataset - cell, char, string - single or multiple names
            %       of datasets to get neighborhoods from
            %   - mfi - char, string - Name of MFI in which to look for
            %       neighborhoods. Any sample provided in dataset, which
            %       does not contain this mfi as a field will be ignored.
            %
            % Output:
            %   - neighs - neighborhoods which are common for all given
            %       datasets. Can be used to access MFIs, process
            %       neighborhoods similarly to actual scans etc.
            
            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end
            if ~exist('mfi', 'var')  % Use most recent dataset
                mfi = 'MFIRSN';
            end
            if ~iscell(dataset)
                dataset = {dataset};
            end

            % Check for non-existent datasets. If all of them all invalid, quit.
            tmp_dataset = {};
            for Di=1:numel(dataset)
                if isfield(app.data, dataset{Di}) && isfield(app.data.(dataset{Di}), mfi)
                    tmp_dataset{end + 1} = dataset{Di};
                end
            end

            if isempty(tmp_dataset)
                neighs = [];
                return;
            else
                dataset = tmp_dataset;
            end

            neighs = app.data.(dataset{1}).(mfi).Properties.VariableNames;
            choices = startsWith(neighs, Constants.neigh_tag);
            neighs = neighs(choices);
            for i = 2:numel(dataset)
                neighs_i = app.data.(dataset{i}).(mfi).Properties.VariableNames;
                choices_i = startsWith(neighs_i, Constants.neigh_tag);
                neighs = intersect(neighs, neighs_i(choices_i), 'stable');
            end
        end

        function ignores = get_ignores(app, dataset, mfi)
            % GET_IGNORES Gives easy access to common columns which do not 
            % have any tag in front of them across multiple datasets
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - dataset - cell, char, string - single or multiple names
            %       of datasets to get `ignores` from
            %   - mfi - char, string - Name of MFI in which to look for
            %       `ignores`. Any sample provided in dataset, which does
            %       not contain this mfi as a field will be ignored.
            %
            % Output:
            %   - ignores - `ignores` which are common for all given
            %       datasets. Can be used to access AllCells, MFIs etc.
            
            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end
            if ~exist('mfi', 'var')  % Use most recent dataset
                mfi = 'MFIRSN';
            end
            if ~iscell(dataset)
                dataset = {dataset};
            end

            % Check for non-existent datasets. If all of them all invalid, quit.
            tmp_dataset = {};
            for Di=1:numel(dataset)
                if isfield(app.data, dataset{Di})
                    tmp_dataset{end + 1} = dataset{Di};
                end
            end

            if isempty(tmp_dataset)
                ignores = [];
                return;
            else
                dataset = tmp_dataset;
            end

            ignores = app.data.(dataset{1}).(mfi).Properties.VariableNames;
            choices = ismember(ignores, Constants.ignore_names);
            ignores = ignores(choices);
            for i = 2:numel(dataset)
                ignores_i = app.data.(dataset{i}).(mfi).Properties.VariableNames;
                choices_i = ismember(ignores_i, Constants.ignore_names);
                ignores = intersect(ignores, ignores_i(choices_i), 'stable');
            end
        end

        function name = valid(app, name, dataset)
            % VALID Enables to get a variables to the form which is
            % accepted by MatLab, while keeping some form of
            % reproductibility by the full function.
            %
            % Input:
            %   - app - Instance of CytoMAP
            %   - name - cell, string, char - Names which are going to be
            %       returned in a valid form
            %   - dataset - Names of sample which defines gate tags needed
            %       in for AllCells, MFIs etc.
            %
            % Output:
            %   - name - cell - Given names, but in a form understandable
            %       by CytoMAP in the given sample.
            %       Note that while it might still run on other samples,
            %       results might not be consistent due to tags denoting
            %       different things.
            
            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end
            if ~iscell(name)
                name = {name};
            end

            name = Helper.gate_full2tag(app, name, dataset);
            name = Helper.valid_channel(name);
            name = Helper.valid_other(name);
            name = Helper.valid_neigh(name);
            name = Helper.valid_var(name);
        end

        function name = valid_var(name)
            % VALID_VAR Enables to get a variables to the form which is
            % accepted by MatLab.
            %
            % Input:
            %   - name - cell, string, char - Names which are going to be
            %       returned in a valid form
            %
            % Output:
            %   - name - cell - Given names, but in a form understandable
            %       by MatLab.
            
            if isempty(name)
                return
            end

            % Makes a name into proper variable name
            name=strrep(name, ' ','_');
            name=strrep(name, ',','_');
            name=strrep(name, '''','');
            name=strrep(name, ':','_');
            name=strrep(name, '(','_');
            name=strrep(name, ')','_');
            name=strrep(name, '.','p');
            name=strrep(name, '*','X');
            name=strrep(name, '"','');
            name=strrep(name, '/','_');
            name=strrep(name, '\','_');
            name=strrep(name, 'Position_','');
            name=strrep(name, 'Position','');
            name=strrep(name, 'Parameter','Para');
            name=strrep(name, '_Para','');
            name=strrep(name, '+', 'POS');
            name=strrep(name, '-', 'NEG');
            name=matlab.lang.makeValidName(name);
            name=regexprep(name, "_+", "_");

            % Cannot start with _
            if startsWith(name, '_')
                name=name(2:end);
            end

            % Might be too long
            if iscell(name)  % Multiple variable names.
                for idx=1:numel(name)
                    if ~isvarname(name{idx})
                        name{idx} = name{idx}(1:namelengthmax);
                    end
                end
            else  % Just a regular, singular variable name
                if ~isvarname(name)
                    name = name(1:namelengthmax);
                end
            end
        end

        function name = valid_gate(name)
            % VALID_GATE Enables to convert exisiting names to form which
            % is acceptable by SAS3D's tree and GateTags table
            %
            % Input:
            %   - name - cell, string, char - Names which are going to be
            %       returned in a valid gate form
            %
            % Output:
            %   - name - cell - Given names, but in a form understandable
            %       by tree/GateTags. It will have a Constants.gate_tag
            %       attached to front of it.
            
            if isempty(name)
                return
            end

            wascell = iscell(name);

            if ~iscell(name)
                name = {name};
            end

            if iscell(name)
                for n_idx=1:numel(name)
                    name{n_idx} = char(name{n_idx});
                end
            end

            only_gates = ~startsWith(name, Constants.other_names);
            % only_gates = only_gates & ~startsWith(name, strcat(Constants.other_tag, Constants.other_names));
            only_gates = only_gates & ~ismember(name, Constants.ignore_names);
            only_gates = only_gates & ~startsWith(name, Constants.other_tag);
            only_gates = only_gates & ~startsWith(name, Constants.channel_tag);
            only_gates = only_gates & ~startsWith(name, Constants.neigh_tag);

            if ~startsWith(name(only_gates), Constants.gate_tag)
                name(only_gates) = strcat(Constants.gate_tag, name(only_gates));
            end

            name(only_gates) = Helper.valid_var(name(only_gates));

            name = regexprep(name, strcat('(', Constants.gate_tag, ')+'), Constants.gate_tag);

            if ~wascell
                name = name{1};
            end
        end

        function name = valid_channel(name)
            % VALID_CHANNEL Enables to convert exisiting names to form
            % which is acceptable by SAS3D
            %
            % Input:
            %   - name - cell, string, char - Names which are going to be
            %       returned in a valid channel form
            %
            % Output:
            %   - name - cell - Given names, but in a form understandable
            %       by SAS3D. It will have a Constants.channel_tag
            %       attached to front of it.
            
            if isempty(name)
                return
            end

            wascell = iscell(name);

            if ~iscell(name)
                name = {name};
            end

            if iscell(name)  % Sanity
                for n_idx=1:numel(name)
                    name{n_idx} = char(name{n_idx});
                end
            end

            only_channels = ~startsWith(name, Constants.other_names);
            only_channels = only_channels & ~ismember(name, Constants.ignore_names);
            only_channels = only_channels & ~startsWith(name, Constants.other_tag);
            only_channels = only_channels & ~startsWith(name, Constants.gate_tag);
            only_channels = only_channels & ~startsWith(name, Constants.neigh_tag);

            % Start Each channel with Ch
            for n=1:numel(name)
                if only_channels(n) && (startsWith(name(n), lower(Constants.channel_tag)) || startsWith(name(n), upper(Constants.channel_tag)))
                    name_n = name{n};
                    name_n = name_n(3:end);
                    name(n) = {name_n};
                end
            end

            % If it starts with some typical channel tag, which is not currently used, then remove that.
            for str_i=1:numel(Constants.channel_starts)
                if strcmp(Constants.channel_starts(str_i), Constants.channel_tag)
                    continue;
                end
                start_i = Constants.channel_starts{str_i};
                for n_i=1:numel(name)
                    if ~only_channels(n_i)
                        continue;
                    end

                    while startsWith(name{n_i}, string(start_i))
                        name{n_i} = name{n_i}(numel(start_i) + 1:end);
                    end
                end
            end
            % Add a tag if not yet added (only one though), and convert it to valid var in MatLab
            if ~startsWith(name(only_channels), Constants.channel_tag)
                name(only_channels) = strcat(Constants.channel_tag, name(only_channels));
            end

            name(only_channels) = Helper.valid_var(name(only_channels));
            name = regexprep(name, strcat('(', Constants.channel_tag, ')+'), Constants.channel_tag);
            
            % Unless it's X,Y,Z, etc.
            x = strcmp(name, strcat(Constants.channel_tag, 'X'));
            y = strcmp(name, strcat(Constants.channel_tag, 'Y'));
            z = strcmp(name, strcat(Constants.channel_tag, 'Z'));
            vol = strcmp(name, 'Neigh Volume');
            effvol = strcmp(name, 'Effective Neigh Volume');

            name(x) = {'X'};
            name(y) = {'Y'};
            name(z) = {'Z'};
            name(vol) = {'Neigh_Volume'};
            name(effvol) = {'Effective_Neigh_Volume'};
            

            if ~wascell
                name = name{1};
            end
        end

        function name = valid_other(name)
            % VALID_OTHER Enables to convert exisiting names to form
            % which is acceptable by SAS3D
            %
            % Input:
            %   - name - cell, string, char - Names which are going to be
            %       returned in a valid other form
            %
            % Output:
            %   - name - cell - Given names, but in a form understandable
            %       by SAS3D. It will have a Constants.other_tag attached
            %       to front of it.
            
            if isempty(name)
                return
            end

            wascell = iscell(name);

            if ~iscell(name)
                name = {name};
            end

            if iscell(name)  % Sanity
                for n_idx=1:numel(name)
                    name{n_idx} = char(name{n_idx});
                end
            end

            only_other = ~startsWith(name, Constants.gate_tag);
            only_other = only_other & ~ismember(name, Constants.ignore_names);
            only_other = only_other & ~startsWith(name, Constants.channel_tag);
            only_other = only_other & ~startsWith(name, Constants.neigh_tag);

            proc_name = name(only_other);
            for n=1:numel(proc_name)
                if any(startsWith(proc_name{n}, Constants.other_names))
                    % Start Each channel with Ch
                    if ~startsWith(proc_name{n}, Constants.other_tag)
                        proc_name{n} = strcat(Constants.other_tag, proc_name{n});
                    end
                    proc_name{n} = Helper.valid_var(proc_name{n});
                end
            end

            name(only_other) = proc_name;
            name = regexprep(name, strcat('(', Constants.other_tag, ')+'), Constants.other_tag);

            if ~wascell
                name = name{1};
            end
        end

        function name = valid_neigh(name)
            % VALID_NEIGH Enables to convert exisiting names to form
            % which is acceptable by SAS3D
            %
            % Input:
            %   - name - cell, string, char - Names which are going to be
            %       returned in a valid neigh form
            %
            % Output:
            %   - name - cell - Given names, but in a form understandable
            %       by SAS3D. It will have a Constants.neigh_tag attached
            %       to front of it.
            
            if isempty(name)
                return
            end

            wascell = iscell(name);

            if ~iscell(name)
                name = {name};
            end

            if iscell(name)  % Sanity
                for n_idx=1:numel(name)
                    name{n_idx} = char(name{n_idx});
                end
            end

            only_neighs = ~startsWith(name, Constants.other_names);
            % only_neighs = only_neighs & ~startsWith(name, strcat(Constants.other_tag, Constants.other_names));
            only_neighs = only_neighs & ~ismember(name, Constants.ignore_names);
            only_neighs = only_neighs & ~startsWith(name, Constants.other_tag);
            only_neighs = only_neighs & ~startsWith(name, Constants.channel_tag);
            only_neighs = only_neighs & ~startsWith(name, Constants.gate_tag);

            if ~startsWith(name(only_neighs), Constants.neigh_tag)
                name(only_neighs) = strcat(Constants.neigh_tag, name(only_neighs));
            end

            name(only_neighs) = Helper.valid_var(name(only_neighs));
            name = regexprep(name, strcat('(', Constants.neigh_tag, ')+'), Constants.neigh_tag);

            if ~wascell
                name = name{1};
            end
        end

        function name = full(app, name, dataset)
            % FULL Enables to get a variables to the form which is more
            % human readable, than the one accepted by MatLab, while
            % keeping some form of reproductibility by the valid function.
            %
            % Input:
            %   - app - Instance of SAS3D
            %   - name - cell, string, char - Names which are going to be
            %       returned in a full form
            %   - dataset - Names of sample which defines gate tags to be
            %       transferred back to user understandable format.
            %
            % Output:
            %   - name - cell - Given names, but in a human readable form.
            %       Note that while it might still run on other samples,
            %       results might not be consistent due to tags denoting
            %       different things.
            
            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end
            if ~iscell(name)
                name = {name};
            end
            % Gate tag to Cell name could be specific per sample
            name(startsWith(name, Constants.gate_tag)) = Helper.gate_tag2full(app,name(startsWith(name, Constants.gate_tag)),dataset);
            % remove the channel tag
            name(startsWith(name, Constants.channel_tag)) = Helper.full_channel(name(startsWith(name, Constants.channel_tag)));
            % remove the 'other' tag
            name(startsWith(name, Constants.other_tag)) = Helper.full_other(name(startsWith(name, Constants.other_tag)));
        end

        function name = full_var(name)
            % FULL_VAR Enables to get a variables to the human readable
            % form.
            %
            % Input:
            %   - name - cell, string, char - Names which are going to be
            %       returned in human readable form
            %
            % Output:
            %   - name - cell - Given names, but in human readable form.
            
            name = strrep(name, '_', ' ');
            name = strrep(name, 'NEG', '-');
            name = strrep(name, 'POS', '+');
        end

        function name = full_gate(name)
            % FULL_GATE Enables to get a gates to the human readable
            % form. Given an array of multiple different elements (other,
            % gates, channels etc.) it will convert only the ones starting
            % with the gate tag.
            %
            % Input:
            %   - name - cell, string, char - Gates which are going to be
            %       returned in human readable form
            %
            % Output:
            %   - name - cell - Given names, but in human readable form.
            
            if iscell(name)
                for n_idx=1:numel(name)
                    name{n_idx} = char(name{n_idx});
                end
            end

            idx = startsWith(name, Constants.gate_tag);
            name = strrep(name, Constants.gate_tag, '');
            name(idx) = Helper.full_var(name(idx));
        end

        function name = full_channel(name)
            % FULL_CHANNEL Enables to get a channel to the human readable
            % form. Given an array of multiple different elements (other,
            % gates, channels etc.) it will convert only the ones starting
            % with the channel tag.
            %
            % Input:
            %   - name - cell, string, char - Channel which are going to be
            %       returned in human readable form
            %
            % Output:
            %   - name - cell - Given names, but in human readable form.
            
            idx = startsWith(name, Constants.channel_tag);
            name = strrep(name, Constants.channel_tag, '');
            name(idx) = Helper.full_var(name(idx));
        end

        function name = full_other(name)
            % FULL_OTHER Enables to get a others to the human readable
            % form. Given an array of multiple different elements (other,
            % gates, channels etc.) it will convert only the ones starting
            % with the other tag.
            %
            % Input:
            %   - name - cell, string, char - Others which are going to be
            %       returned in human readable form
            %
            % Output:
            %   - name - cell - Given names, but in human readable form.
            
            idx = startsWith(name, Constants.other_tag);
            name = strrep(name, Constants.other_tag, '');
            name(idx) = Helper.full_var(name(idx));
        end

        function name = full_neigh(name)
            % FULL_NEIGH Enables to get a neighborhoods to the human
            % readable form. Given an array of multiple different elements
            % (other, gates, channels etc.) it will convert only the ones
            % starting with the neighborhood tag.
            %
            % Input:
            %   - name - cell, string, char - Neighborhood which are going
            %       to be returned in human readable form
            %
            % Output:
            %   - name - cell - Given names, but in human readable form.
            
            idx = startsWith(name, Constants.neigh_tag);
            name = strrep(name, Constants.neigh_tag, '');
            name(idx) = Helper.full_var(name(idx));
        end

        function name = gate_tag2full(app, name, dataset)
            % GATE_TAG2FULL Given a handle to the app,
            %                 dataset (or sample) name,
            %                 and full name of the gate,
            %                 returns a path corresponding to this gate name in this sample.
            %
            % Input:
            %   - app     - Handle to a main struct of SAS3D instance.
            %   - name    - A tag (or cell array) of gate, which names are to be
            %               returned.
            %   - dataset - A name of sample in which these gates exist.
            %               There can be only one dataset given,
            %               since otherwise behaviour would not be specified
            %               (tags can correspond to different gates across
            %                different samples).
            %
            % Output:
            %   - name    - Cell containing paths corresponding to the names
            %               of gates in the dataset.
            %               If no name corresponding to a tag exists a
            %               tag itself is kept.

            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end

            if iscell(dataset)
                if numel(dataset) == 1
                    dataset = dataset{1};
                else
                    error("Dataset cannot be a cell with size above 1, since tag is undefined for such case." + newline + ...
                        "Currently it is {" + string(join(dataset, ", ")) + "}.");
                end
            end

            if iscell(name)
                for i = 1:numel(name)
                    index = find(strcmp(name{i}, app.data.(dataset).GateTags.Properties.VariableNames));
                    if ~isempty(index)
                        name{i} = app.data.(dataset).GateTags{2, index}{1};
                    end
                end
            else
                index = find(strcmp({name}, app.data.(dataset).GateTags.Properties.VariableNames));
                if ~isempty(index)
                    name = app.data.(dataset).GateTags{2, index}{1};
                end
            end
        end

        function name = gate_full2tag(app, name, dataset)
            % GATE_FULL2TAG DEPRACATED. USE get_tag(app, name, dataset)
            %       (unless you explictly do not want to create new tags)
            %
            % gate_full2tag given a handle to the app,
            %                 dataset (or sample) name,
            %                 and full name of gate,
            %                 returns a tag corresponding to this gate in this sample.
            % Input:
            %   - app     - Handle to a main struct of SAS3D instance.
            %   - name    - A name (or cell array) of gate, which tags are to be
            %               returned
            %   - dataset - A name of sample in which these gates exist.
            %               There can be only one dataset given,
            %               since otherwise behaviour would not be specified
            %               (tags can correspond to different gates across
            %                different samples).
            %
            % Output:
            %   - name    - Cell containing tags corresponding to the names of dataset.
            %               If no tag corresponding to a name exists a new one is created.
            %
            % Notes:
            %   - If a gate name doesn't have a corresponding tag,
            %       then it's name is returned, unchanged.
            %       This contrasts with the behavior of get_tag,
            %       which would return just a new tag.

            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end

            if iscell(dataset)
                if numel(dataset) == 1
                    dataset = dataset{1};
                else
                    error("Dataset cannot be a cell with size above 1, since tag is undefined for such case." + newline + ...
                        "Currently it is {" + string(join(dataset, ", ")) + "}.");
                end
            end

            if iscell(name)
                for i = 1:numel(name)
                    index = find(strcmp(name{i}, app.data.(dataset).GateTags{2, :}));
                    if ~isempty(index)
                        name{i} = app.data.(dataset).GateTags.Properties.VariableNames{index};
                    end
                end
            else
                if strcmp(name, 'All Cells')
                    name = strcat(Constants.gate_tag, '1');
                    return;
                end
                index = find(strcmp({name}, app.data.(dataset).GateTags{2, :}));
                if ~isempty(index)
                    name = app.data.(dataset).GateTags.Properties.VariableNames{index};
                end
            end
        end

        function name = gate_full2path(app, name, dataset)
            % GATE_FULL2PATH Given a handle to the app,
            %                 dataset (or sample) name,
            %                 and full name of the gate,
            %                 returns a path corresponding to this gate name in this sample.
            %
            % Input:
            %   - app     - Handle to a main struct of SAS3D instance.
            %   - name    - A tag (or cell array) of gate, which paths are to be
            %               returned.
            %   - dataset - A name of sample in which these gates exist.
            %               There can be only one dataset given,
            %               since otherwise behaviour would not be specified
            %               (tags can correspond to different gates across
            %                different samples).
            %
            % Output:
            %   - name    - Cell containing paths corresponding to the names
            %               of gates in the dataset.
            %               If no paths corresponding to a name exists a
            %               name itself is kept.

            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end

            if iscell(name)
                for i = 1:numel(name)
                    index = find(strcmp(name{i}, app.data.(dataset).GateTags{2, :}));
                    if ~isempty(index)
                        name{i} = app.data.(dataset).GateTags{1, index};
                    end
                end
            else
                index = find(strcmp({name}, app.data.(dataset).GateTags{2, :}));
                if ~isempty(index)
                    name = app.data.(dataset).GateTags{1, index};
                end
            end
        end

        function name = gate_tag2parent_tag(app, name, dataset)
            % GATE_TAG2PARENT Given a handle to the app,
            %                 dataset (or sample) name,
            %                 and gate tag,
            %                 returns a full name corresponding to this tag in this sample.
            %
            % Input:
            %   - app     - Handle to a main struct of SAS3D instance.
            %   - name    - A tag (or cell array) of gate, which parents tag
            %               are to be returned
            %   - dataset - A name of sample in which these gates exist.
            %               There can be only one dataset given,
            %               since otherwise behaviour would not be specified
            %               (tags can correspond to different gates across
            %                different samples).
            %
            % Output:
            %   - name    - Cell containing tags of parents of the gates
            %               which have been given in the given dataset.
            
            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end

            if iscell(name)
                for i = 1:numel(name)
                    index = find(strcmp(name{i}, app.data.(dataset).GateTags.Properties.VariableNames));
                    if ~isempty(index)
                        path = app.data.(dataset).GateTags{1, index}{1};
                        path = split(path, '/');
                        path = path(1:end-1);
                        path = join(path, '/');

                        if iscell(path)
                            path = path{1};
                        end

                        if ~contains(path, '/')
                            name{i} = strcat(Constants.gate_tag, '0');
                        else
                            idx = strcmp(path, app.data.(dataset).GateTags{1, :});
                            name{i} = app.data.(dataset).GateTags.Properties.VariableNames{idx};
                        end

                    end
                end
            else
                index = find(strcmp({name}, app.data.(dataset).GateTags.Properties.VariableNames));
                if ~isempty(index)
                    path = app.data.(dataset).GateTags{1, index}{1};
                    path = split(path, '/');
                    path = path(1:(end-1));
                    path = join(path, '/');

                    if iscell(path)
                        path = path{1};
                    end
                    
                    if ~contains(path, '/')
                        name = strcat(Constants.gate_tag, '0');
                    else
                        idx = strcmp(path, app.data.(dataset).GateTags{1, :});
                        name = app.data.(dataset).GateTags.Properties.VariableNames{idx};
                    end

                end
            end
        end

        function tag = get_tag(app, name, dataset)
            % GET_TAG Given a handle to the app,
            %                 dataset (or sample) name,
            %                 and full name of gate,
            %                 returns a tag corresponding to this gate in this sample.
            %
            % Input:
            %   - app     - Handle to a main struct of SAS3D instance.
            %   - name    - A short/full name (or cell array) of gate, which tags are to be
            %               returned
            %   - dataset - A name of sample in which these gates exist.
            %               There can be only one dataset given,
            %               since otherwise behaviour would not be specified
            %               (tags can correspond to different gates across
            %                different samples).
            %
            % Output:
            %   - name    - Cell containing tags corresponding to the names of dataset.
            %               If no tag corresponding to a name exists a new one is created.
            %
            % Notes:
            %   - Name can be a cell of multiple gate names.
            %   - If a gate name doesn't have a corresponding tag,
            %       then it's name is returned, unchanged.
            %       This contrasts with the behavior of gate_full2tag,
            %       which would return just name it was given in the first place.

            if ~exist('dataset', 'var')  % Use most recent dataset
                dataset = app.DataN.Value;
            end
            if ~iscell(name)
                name = {name};
            end

            if iscell(dataset)
                if numel(dataset) == 1
                    dataset = dataset{1};
                else
                    error("Dataset cannot be a cell with size above 1, since tag is undefined for such case." + newline + ...
                        "Currently it is {" + string(join(dataset, ", ")) + "}.");
                end
            end

            % Create array of empty tag. Have to be chars, becuase otherwise ismember later on will not work.
            tag = cell(numel(name), 1);
            for t_i=1:numel(tag)
                tag(t_i) = {''};
            end

            for n=1:numel(name)
                n_name = name(n);
                if strcmp(n_name, 'All Cells')
                    tag{n} = strcat(Constants.gate_tag, '1');
                    continue;
                end
                idx = ~strcmp(n_name, app.data.(dataset).GateTags{2, :});

                if all(idx) % Make a new tag
                    tag_idx = numel(idx);
                    % Check if such tag already exists either in GateTags, or in already, newly created tags.
                    while ismember({strcat(Constants.gate_tag, num2str(tag_idx))}, app.data.(dataset).GateTags.Properties.VariableNames) ...
                            || ismember({strcat(Constants.gate_tag, num2str(tag_idx))}, tag)
                        tag_idx = tag_idx + 1;
                    end
                    tag{n} = strcat(Constants.gate_tag, num2str(tag_idx));
                else % Return already existing tag
                    tag{n} = app.data.(dataset).GateTags.Properties.VariableNames{~idx};
                end
            end
        end

        function tab = reorder_cols(tab)
            % REORDER_COLS Given table places X, Y, Z axis to be first 3 columns.
            % If Z column doesn't exist, then it is added to the table.
            
            if ~istable(tab) || any(~ismember({'X', 'Y'}, tab.Properties.VariableNames))
                return
            end
            if ~ismember('Z', tab.Properties.VariableNames)
                tab.Z = 0 .* tab.X;
            end

            tab = movevars(tab, 'X', 'Before', 1);
            tab = movevars(tab, 'Y', 'After', 'X');
            tab = movevars(tab, 'Z', 'After', 'Y');
        end

        function tab = populate_table(app, varargin)
            % POPULATE_TABLE
            % Input:
            %     - app - Handle to current instance of SAS3D
            %     - smpls - keyword/optional - which samples to use from app
            %     - mfi - keyword/optional -
            %              Whether to make columns corresponding to given MFI scan
            %              Either 'MFIRSN', 'MFICCN' or 'AllCells'
            %     - prev_table - keyword/optional -
            %              If it will be an update to table, it can be given
            %              to retain weights, and choices for Channels and MFI
            %              across update.
            %     - fill_weight - keyword/optional -
            %              Default value to fill 1st and 4th (if mfi is true) with.
            %              Has to be a number.
            %     - fill_checkbox - keyword/optional -
            %              Default value to fill 2nd and 5th (if mfi is true) with.
            %              Has to be a logical value.
            % Output:
            %     - tab - Table with either 5 (no MFI) or 8 (with MFI) columns.
            %             It contains weights, choices and names for Channels/MFI,
            %             and choices, names for Samples.
            % Notes:
            %     - All keywords are case insensitive
            
            %% Function intro - Read args
            tab = [];
            if isempty(app.data)
                errordlg('Data does not exist. Cannot produce a table');
                return;
            end

            defaults = struct;
            defaults.smpls = {app.DataN.Value};
            defaults.prev_table = {};
            defaults.mfi = {};
            defaults.fill_weight = 1;
            defaults.fill_checkbox = true;

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

            if ~iscell(defaults.smpls)
                defaults.smpls = {defaults.smpls};
            end

            use_smpl = isfield(app.data, defaults.smpls);
            defaults.smpls = defaults.smpls(use_smpl);

            if iscell(defaults.fill_weight)
                defaults.fill_weight = defaults.fill_weight{1};
            end
            if islogical(defaults.fill_weight)
                defaults.fill_weight = 1.0 * defaults.fill_weight;
            end
            if ~isnumeric(defaults.fill_weight)
                error('fill_weight argument has to be a number.');
            end

            if iscell(defaults.fill_checkbox)
                defaults.fill_checkbox = defaults.fill_checkbox{1};
            end
            if isnumeric(defaults.fill_checkbox) && defaults.fill_checkbox == 1
                defaults.fill_checkbox = true;
            elseif isnumeric(defaults.fill_checkbox) && defaults.fill_checkbox == 0
                defaults.fill_checkbox = false;
            end
            if ~islogical(defaults.fill_checkbox)
                error('fill_checkbox argument has to be a logical value.');
            end

            if ~iscell(defaults.prev_table)
                error('prev_table argument has to be a cell.');
            elseif ~isempty(defaults.prev_table)
                % Check cardinality
                if ~isempty(defaults.mfi)
                    if size(defaults.prev_table, 2) ~= 8
                        error('prev_table given mfi argument has to have 8 columns');
                    end
                elseif size(defaults.prev_table, 2) ~= 5
                    error('prev_table given NO (or empty) mfi argument has to have 5 columns');
                end
                % Check formatting
                if ~all(cellfun('isempty', defaults.prev_table(:, 1)) | cellfun(@isnumeric, defaults.prev_table(:, 1)))
                    error('1st column has to have all elements either numeric, or empty');
                end
                if ~all(cellfun('isempty', defaults.prev_table(:, 2)) | cellfun(@islogical, defaults.prev_table(:, 2)))
                    error('2nd column has to have all elements either logical, or empty');
                end
                if ~all(cellfun('isempty', defaults.prev_table(:, 3)) | cellfun(@isstring, defaults.prev_table(:, 3)) | cellfun(@ischar, defaults.prev_table(:, 3)))
                    error('3rd column has to have all elements either string, char, or empty');
                end
                if ~isempty(defaults.mfi)
                    if ~all(cellfun('isempty', defaults.prev_table(:, 4)) | cellfun(@isnumeric, defaults.prev_table(:, 4)))
                        error('Given mfi argument, 4th column has to have all elements either numeric, or empty');
                    end
                    if ~all(cellfun('isempty', defaults.prev_table(:, 5)) | cellfun(@islogical, defaults.prev_table(:, 5)))
                        error('Given mfi argument, 5th column has to have all elements either logical, or empty');
                    end
                    if ~all(cellfun('isempty', defaults.prev_table(:, 6)) | cellfun(@isstring, defaults.prev_table(:, 6)) | cellfun(@ischar, defaults.prev_table(:, 6)))
                        error('Given mfi argument, 6th column has to have all elements either string, char, or empty');
                    end
                else
                    if ~all(cellfun('isempty', defaults.prev_table(:, 4)) | cellfun(@islogical, defaults.prev_table(:, 4)))
                        error('Given NO (or empty) mfi argument, 4th column has to have all elements either logical, or empty');
                    end
                    if ~all(cellfun('isempty', defaults.prev_table(:, 5)) | cellfun(@isstring, defaults.prev_table(:, 5)) | cellfun(@ischar, defaults.prev_table(:, 5)))
                        error('Given NO (or empty) mfi argument, 5th column has to have all elements either string, char, or empty');
                    end
                end
            end

            %% Get gates for current table. Create a default values for the new table
            [~, cnames_pt] = Helper.get_gates(app, defaults.smpls);

            keep_cn = cell(numel(cnames_pt), 2);
            keep_cn(1:end, 1) = {defaults.fill_weight};
            if strcmp(defaults.fill_checkbox, "majority")
                if isempty(defaults.prev_table)
                    error("majority/minority options are invalid for the fill_checkbox argument" + ...
                          " if the prev_table argument is empty");
                end
                prev_idxs = ~cellfun('isempty', defaults.prev_table);
                ratio = sum(mat2num(defaults.prev_table(prev_idxs, 2))) * 1.0 / sum(prev_idxs);
                keep_cn(1:end, 2) = {ratio >= 0.5};
            elseif strcmp(defaults.fill_checkbox, "minority")
                if isempty(defaults.prev_table)
                    error("majority/minority options are invalid for the fill_checkbox argument" + ...
                          " if the prev_table argument is empty");
                end
                prev_idxs = ~cellfun('isempty', defaults.prev_table);
                ratio = sum(mat2num(defaults.prev_table(prev_idxs, 2))) * 1.0 / sum(prev_idxs);
                keep_cn(1:end, 2) = {ratio < 0.5};
            else
                if ~islogical(defaults.fill_checkbox)
                    error('fill_checkbox argument has to be a logical array.')
                end
                keep_cn(1:end, 2) = {logical(defaults.fill_checkbox)};
            end

            %% Process previous table
            % Retain info in the table, given that there was a previous one
            if ~isempty(defaults.prev_table)
                old_names = defaults.prev_table(~cellfun('isempty', defaults.prev_table(:, 3)), 3);
                keep_idxs = ismember(old_names, cnames_pt);

                % A_cn is basically a map between indices, weight and choice
                % Col 1: Index new
                % Col 2: Index old
                % Col 3: Old Weights
                % Col 4: Old Choices
                A_cn = cell(sum(keep_idxs, 'all'), 4);

                A_idx = 1;
                for i=1:numel(keep_idxs)
                    if keep_idxs(i) == 1
                        old_name = defaults.prev_table(i, 3);
                        new_idx = find(strcmp(cnames_pt, old_name));

                        A_cn(A_idx, 1) = {new_idx};
                        A_cn(A_idx, 2) = {i};
                        A_cn(A_idx, 3) = defaults.prev_table(i, 1);
                        A_cn(A_idx, 4) = defaults.prev_table(i, 2);

                        A_idx = A_idx + 1;
                    end
                end

                for i=1:size(A_cn, 1)
                    keep_cn(A_cn{i, 1}, 1) = A_cn(i, 3);
                    keep_cn(A_cn{i, 1}, 2) = A_cn(i, 4);
                end
            end

            % Exclude all of these things from the MFI names
            if ~isempty(defaults.mfi)
                use_smpl = logical(numel(defaults.smpls));
                for i=1:numel(defaults.smpls)
                    if iscell(defaults.smpls{i})
                        celled = defaults.smpls{i};
                        defaults.smpls{i} = celled{1};
                    end
                    if isfield(app.data, defaults.smpls{i}) && isfield(app.data.(defaults.smpls{i}), defaults.mfi)
                        use_smpl(i) = true;
                    else
                        warndlg(strcat("Chosen scan type was not found in sample ", defaults.smpls{i}, ". It will be omitted."));
                        use_smpl(i) = false;
                    end
                end
                defaults.smpls = defaults.smpls(use_smpl);
                if isempty(defaults.smpls)
                    errordlg("Neither of the given samples has a requested MFI. Aborting.");
                    return;
                end

                MFINMS = app.data.(defaults.smpls{1}).(defaults.mfi).Properties.VariableNames;
                for smpl_idx=1:numel(defaults.smpls)
                    MFINMS = intersect(MFINMS, app.data.(defaults.smpls{smpl_idx}).(defaults.mfi).Properties.VariableNames, 'stable');
                end

                MFINMS = MFINMS(~startsWith(MFINMS, Constants.gate_tag));
                MFINMS = Helper.full_channel(MFINMS);
                % Comment out those, since they do not work well with rest of program
                % MFINMS = Helper.full_other(MFINMS);
                % MFINMS = Helper.full_neigh(MFINMS);
                MFINMS(ismember(MFINMS, Constants.ignore_names)) = Helper.full_var(MFINMS(ismember(MFINMS, Constants.ignore_names)));

                keep_MFI = cell(numel(MFINMS), 2);
                keep_MFI(:, 1) = {defaults.fill_weight};
                if strcmp(defaults.fill_checkbox, "majority")
                    if isempty(defaults.prev_table)
                        error("majority/minority options are invalid for the fill_checkbox argument" + ...
                              " if the prev_table argument is empty");
                    elseif size(defaults.prev_table) < 8
                        error("majority/minority options require prev_table to contain MFI part, if mfi argument is true");
                    end
                    prev_idxs = ~cellfun('isempty', defaults.prev_table);
                    ratio = sum(mat2num(defaults.prev_table(prev_idxs, 5))) * 1.0 / sum(prev_idxs);
                    keep_MFI(1:end, 2) = {ratio >= 0.5};
                elseif strcmp(defaults.fill_checkbox, "minority")
                    if isempty(defaults.prev_table)
                        error("majority/minority options are invalid for the fill_checkbox argument" + ...
                              " if the prev_table argument is empty");
                    elseif size(defaults.prev_table) < 8
                        error("majority/minority options require prev_table to contain MFI part, if mfi argument is true");
                    end
                    prev_idxs = ~cellfun('isempty', defaults.prev_table);
                    ratio = sum(mat2num(defaults.prev_table(prev_idxs, 5))) * 1.0 / sum(prev_idxs);
                    keep_MFI(1:end, 2) = {ratio < 0.5};
                else
                    if ~islogical(defaults.fill_checkbox)
                        error('fill_checkbox argument has to be a logical array.')
                    end
                    keep_MFI(1:end, 2) = {logical(defaults.fill_checkbox)};
                end

                % Make sure that previous table is given
                if ~isempty(defaults.prev_table)
                    old_names = defaults.prev_table(~cellfun('isempty', defaults.prev_table(:, 6)), 6);
                    keep_idxs = ismember(old_names, MFINMS);

                    % A_cn is basically a map between indices, weight and choice
                    % Col 1: Index new
                    % Col 2: Index old
                    % Col 3: Old Weights
                    % Col 4: Old Choices
                    A_MFI = cell(sum(keep_idxs, 'all'), 4);

                    A_idx = 1;
                    for i=1:numel(keep_idxs)
                        if keep_idxs(i) == 1
                            old_name = defaults.prev_table(i, 6);
                            new_idx = find(strcmp(MFINMS, old_name));

                            A_MFI(A_idx, 1) = {new_idx};
                            A_MFI(A_idx, 2) = {i};
                            A_MFI(A_idx, 3) = defaults.prev_table(i, 4);
                            A_MFI(A_idx, 4) = defaults.prev_table(i, 5);

                            A_idx = A_idx + 1;
                        end
                    end

                    for i=1:size(A_MFI, 1)
                        keep_MFI(A_MFI{i, 1}, 1) = A_MFI(i, 3);
                        keep_MFI(A_MFI{i, 1}, 2) = A_MFI(i, 4);
                    end
                end

                tab = cell(max(max(numel(MFINMS), numel(cnames_pt)), numel(app.DataN.Items)),5);

                % Channel Names
                cnames_pt = Helper.full_channel(cnames_pt);
                tab(1:numel(cnames_pt), 1) = keep_cn(:, 1);
                tab(1:numel(cnames_pt), 2) = keep_cn(:, 2);
                tab(1:numel(cnames_pt), 3) = cnames_pt;

                % MFI's
                tab(1:numel(MFINMS), 4) = keep_MFI(:, 1);
                tab(1:numel(MFINMS), 5) = keep_MFI(:, 2);
                tab(1:numel(MFINMS), 6) = MFINMS;

                % Samples
                tab(1:numel(app.DataN.Items), 7) = num2cell(ismember(app.DataN.Items, defaults.smpls));
                tab(1:numel(app.DataN.Items), 8) = app.DataN.Items;
                
            else
                tab = cell(max(numel(cnames_pt), numel(app.DataN.Items)),5);

                % Channel Names
                cnames_pt = Helper.full_channel(cnames_pt);
                tab(1:numel(cnames_pt), 1) = keep_cn(:, 1);
                tab(1:numel(cnames_pt), 2) = keep_cn(:, 2);
                tab(1:numel(cnames_pt), 3) = cnames_pt;

                % Samples
                tab(1:numel(app.DataN.Items), 4) = num2cell(ismember(app.DataN.Items, defaults.smpls));
                tab(1:numel(app.DataN.Items), 5) = app.DataN.Items;
            end
        end

        function loaded = any_sample(app)
            % ANY_SAMPLE Given SAS3D instance, check if there is any sample/dataset
            % loaded in.
            
            loaded = ~isempty(app.data) && ~isempty(app.DataN.Value);
            if ~loaded
                errordlg("In order to use this option load any Sample or Dataset first.", "Data Not Loaded");
            end
        end

        function loaded = any_net(app)
            % ANY_NET Given SAS3D instance, check if there any model has been
            % trained yet.
            
            loaded = ~isempty(app.net) && ~isempty(fieldnames(app.net));
            if ~loaded
                errordlg("In order to use this option you have to train any model.", "Model Not Loaded");
            end
        end

        function net = get_net(app, net_tag)
            % GET_NET Gives functional access to SAS3Ds models which have
            % been trained by user in this session.
            %
            % Input:
            %   - app - Instance of SAS3D
            %   - net_tag - cell, char, string - Either single or multiple
            %       names of nets in valid form. If such name doesn't exist
            %       under app.net (location where trained models are put
            %       in), then it's gonna be ignored.
            %
            % Output:
            %   - net - cell, struct - Either struct corresponding to a
            %       name given, or a cell of structs corresponding to names
            %       given. If a name doesn't exist in app.net then an empty
            %       array is returned on that index of cell array, or
            %       instead of struct.
            
            if ~exist('net_tag', 'var') || isempty(net_tag)
                net = [];
                return;
            end
            wascell = iscell(net_tag);
            if ~iscell(net_tag)
                net_tag = {net_tag};
            end

            % Check if starts with other tag. If so, remove.
            for net_idx=1:numel(net_tag)
                if startsWith(net_tag{net_idx}, Constants.other_tag)
                    net_tag{net_idx} = net_tag{net_idx}(length(Constants.other_tag) + 1:end);
                end
            end

            % Check if starts with one of the other names (Basically either RegionCCN or RegionRSN)
            % If so remove.
            for net_idx=1:numel(net_tag)
                for n=1:numel(Constants.other_names)
                    if startsWith(net_tag{net_idx}, Constants.other_names{n})
                        net_tag{net_idx} = net_tag{net_idx}(length(Constants.other_names{n}) + 1:end);
                    end
                end
            end

            net = cell(numel(net_tag), 1);
            for net_idx = 1:numel(net)
                if isfield(app.net, net_tag{net_idx})
                    net{net_idx} = app.net.(net_tag{net_idx});
                else
                    net{net_idx} = [];
                end
            end

            if ~wascell
                net = net{1};
            end
        end
        
        function surf_names = get_common_surfaces(app, smpls, type)
            % GET_COMMON_SURFACES Returns common surfaces of given type
            % across all given samples
            
            if ~iscell(smpls)
                smpls = {smpls};
            end
            surf_names = fieldnames(app.data.(smpls{1}).Surfaces.(type));
            for idx=2:numel(smpls)
                surf_names = intersect(surf_names, fieldnames(app.data.(smpls{idx}).Surfaces.(type)), 'stable');
            end
        end

        function [mfi, sample] = find_MFI(app)
            % FIND_MFI Loops over all samples, starting at app.DataN.Value,
            % and finds return a valid MFI, sample combination (i.e. sample
            % and MFI such that app.data.(smpl).(mfi) exists).
            % If both RSN and CCN exists returns mfi='MFIRSN'.
            
            MFI_TYPE = '';
            if isfield(app.data.(app.DataN.Value), 'MFIRSN')
                MFI_TYPE = 'MFIRSN';
            elseif isfield(app.data.(app.DataN.Value), 'MFICCN')
                MFI_TYPE = 'MFICCN';
            else
                MFI_TYPE = 'AllCells';
            end
            mfi = '';
            sample = app.DataN.Value;
            smpl_idx = 1;
            while isempty(MFI_TYPE)
                if smpl_idx > numel(app.DataN.Items)
                    return;
                end
                sample = app.DataN.Items{smpl_idx};
                if isfield(app.data.(sample), 'MFIRSN')
                    MFI_TYPE = 'MFIRSN';
                elseif isfield(app.data.(sample), 'MFICCN')
                    MFI_TYPE = 'MFICCN';
                end
                smpl_idx = smpl_idx + 1;
            end

            mfi = MFI_TYPE;
        end

        function new_smpl = intersect_smpls(app, smpls, vPD)
            % INTERSECT_SMPLS Allows to combine samples into one, which
            % contains all common channels, gates etc. Such sample is ready
            % to be put in the app.data struct
            %
            % Input:
            %   - app - Instance of SAS3D
            %   - smpls - cell - Names of all samples to be combined into
            %       one. These have to exist under app.data, and will be
            %       left unchanged.
            %   - vPD - waitbar - optional - If given then this waitbar
            %       will be updated as samples are merged.
            %
            % Output:
            %   - new_smpl - struct - Struct which has a form of generic,
            %       newly loaded sample. It has AllCells, GateTags, tree
            %       and all of the other fields necessary for being
            %       processed by SAS3D. It's not yet added to the
            %       app.data however, which has to be done in place which
            %       calls this function.
            %
            % Note:
            %   - MFIs are not merged, due to a problem regarding
            %       interpolating neighborhoods with each other between
            %       samples. This is out of scope yet of the SAS3D, but
            %       may be added in later versions. For a current solution
            %       to the task look at Helper.merge_MFI.
            %   - There are plans of adding a union channels instead of
            %       intersect samples as a default behavior of combining
            %       samples, which would fill in missing columns with
            %       either NaNs or 0's. However it involves restructuring
            %       whole Plotting due to MatLab mishandling of such cases.
            
            if ~exist('vPD', 'var')
                vPD = [];
            end

            %% Loop over the previous sample AllCells and combine them into 1 AllCells
            % Get names of channels, others, gates etc.
            channels = Helper.get_channels(app, smpls);
            others = Helper.get_others(app, smpls);
            % Pull the cell names, making sure to exclude the file name
            paths_short = Helper.remove_smpl_from_paths(app.data.(smpls{1}).GateTags{1, :});  % Case of multiple samples

            for s_i=2:numel(smpls)
                paths_short_i = Helper.remove_smpl_from_paths(app.data.(smpls{s_i}).GateTags{1, :});
                                
                % Find the common cell types in the two samples
                [~,~, INDCons]  = intersect(paths_short, paths_short_i, 'stable');
                INDUniq = ~ismember(paths_short_i, paths_short_i(INDCons));
                if sum(INDUniq)~=0
                    paths_short = [paths_short, paths_short_i(INDUniq)];
                end                
            end
            
            % Pull the full path names using the first sample
            smpl_1_paths_short_to_full = Helper.removed_path_to_full(app, paths_short, smpls{1});
            
            if numel(channels) <= 3 && isempty(others) && isempty(paths_short)
                not_matching = questdlg( ...
                    strcat('Two samples are not matching at all.', newline, ...
                        'Only X, Y, Z coordinates will be kept in new sample.', newline, ...
                        'Do you want to proceed?'), ...
                        'Samples not matching', ...
                        'Yes', 'No' ...
                    );
                switch not_matching
                    case 'Yes'

                    otherwise
                        new_smpl = [];
                        return;
                end
            end

            % Main Table is table of channels and others
            [main_table, gate_table, gatetags] = Helper.merge_AllCells(app, smpls, vPD);

            allcells = horzcat(main_table, gate_table);

            [allcells, ~] = unique(allcells, 'rows');
            
            %% Make new Tree
            % We know it will be some subtree of tree from any sample.
            % Since it's intersection of all samples
            the_tree = app.data.(smpls{1}).tree.get_subtree(smpl_1_paths_short_to_full);
            % Get logicals of what to change later on
            old_full = the_tree.full_name;
            old_valid = the_tree.name;
            % Modify naming scheme so it's consistent
            the_tree.full_name = 'All';
            valid_root = Helper.valid_gate('All');
            if iscell(valid_root)
                valid_root = valid_root{1};
            end
            the_tree.name = valid_root;
            the_tree.tag = [Constants.gate_tag, '0'];
            % Adjust Path
            for gate_idx = 1:size(gatetags, 2)
                gate_path = gatetags{1, gate_idx};
                if iscell(gate_path)
                    gate_path = gate_path{1};
                end
                if ~contains(gate_path, '/') && strcmp(gate_path, old_valid)
                    % Replace with new one
                    gatetags{1, gate_idx} = {the_tree.name};
                elseif ~contains(gate_path, '/')
                    % Concatanate new one on top
                    gatetags{1, gate_idx} = {[char(the_tree.name), '/', gate_path]};
                elseif startsWith(gate_path, strcat(old_valid, '/'))
                    % Remove old root, and add new one at the top
                elseif ~startsWith(gate_path, strcat(the_tree.name, '/'))
                    % Add new one at the top
                    if ~startsWith(gate_path, strcat(the_tree.name, '/'))
                        gatetags{1, gate_idx} = {[the_tree.name '/' char(gate_path)]};
                    end
                end
            end
            % Adjust Immidiate Children Short Names
            for gate_idx = 1:size(gatetags, 2)
                gate_path = gatetags{1, gate_idx};
                gate_short = gatetags{2, gate_idx};
                if iscell(gate_path)
                    gate_path = gate_path{1};
                end
                if iscell(gate_short)
                    gate_short = gate_short{1};
                end
                if ~contains(gate_path, '/')
                    gatetags{2, gate_idx} = {the_tree.full_name};
                elseif numel(split(gate_path, '/')) == 2
                    gate_short = split(gate_short, '/');
                    gate_short = [the_tree.full_name '/' char(gate_short(2))];
                    gatetags{2, gate_idx} = {gate_short};
                end
            end
            % Adjust Tags
            if ~isempty(gatetags)
                the_tree = the_tree.retag_kid(gatetags.Properties.VariableNames, gatetags{1, :});
            end
            %% Create the new sample in the main app data structure
            new_smpl = Helper.make_empty_smpl();
            new_smpl.AllCells = allcells;
            new_smpl.GateTags = gatetags;
            new_smpl.tree = the_tree;
        end

        function [channel_other_table, gate_table, gate_tags] = merge_AllCells(app, smpls, vPD)
            % MERGE_ALLCELLS Allows combines AllCells table from multiple
            % samples into one, which contains all common channels, gates
            % etc. It returns gates, and everything else seperately so it
            % is more flexible, but there is a row correspondence which
            % means that one can simple horizontally concatanate them.
            % There is also gate tags table provided, to decode tags from
            % gate table.
            %
            % Input:
            %   - app - Instance of SAS3D
            %   - smpls - cell - Names of all samples to be combined into
            %       one. These have to exist under app.data, and will be
            %       left unchanged.
            %   - vPD - waitbar - optional - If given then this waitbar
            %       will be updated as samples are merged.
            %
            % Output:
            %   - channel_other_table - table - Table containing all the
            %       common channels/other fields between all the samples
            %       for the union of the cells. Same number of rows as
            %       gate_table.
            %   - gate_table - table - Table containing all the 
            %       gates for the samples given. Same number of rows as
            %       channel_other_table. Same number of columns as
            %       gate_tags.
            %   - gate_tags - table - Table with 2 rows. Decodes tags from
            %       gate_table. Column name containing the tag of the gate,
            %       first row path to the gate in a tree, and second
            %       containing a human readable name. Same number of
            %       columns as gate_table.
            
            tic  % Start function timer for troubleshooting
            if ~exist('vPD', 'var')
                vPD = [];
            end
            %% get the channel names that are common to all datasets
            channels = Helper.get_channels(app, smpls);
            others = Helper.get_others(app, smpls);
            [paths, short_names] = Helper.get_gates(app, smpls{1});  % Useful for 1 sample case
            paths_short = Helper.remove_smpl_from_paths(app.data.(smpls{1}).GateTags{1, :});  % Case of multiple samples
            for s_i=2:numel(smpls)
                [pathsi, short_namesi] = Helper.get_gates(app, smpls{s_i});
                paths_short_i = Helper.remove_smpl_from_paths(app.data.(smpls{s_i}).GateTags{1, :});

                % Find the common cell types in the two samples
                [~,~, INDCons]  = intersect(paths, pathsi, 'stable');
                INDUniq = ~ismember(pathsi, pathsi(INDCons));
                if sum(INDUniq)~=0
                    paths = [paths, pathsi(INDUniq)];
                end 
                
                % Find the common cell types in the two samples
                [~,~, INDCons]  = intersect(short_names, short_namesi, 'stable');
                INDUniq = ~ismember(short_namesi, short_namesi(INDCons));
                if sum(INDUniq)~=0
                    short_names = [short_names, short_namesi(INDUniq)];
                end 
                
                % Find the common cell types in the two samples
                [~,~, INDCons]  = intersect(paths_short, paths_short_i, 'stable');
                INDUniq = ~ismember(paths_short_i, paths_short_i(INDCons));
                if sum(INDUniq)~=0
                    paths_short = [paths_short, paths_short_i(INDUniq)];
                end    
                
            end
            
            %% Initialize Tables
            channel_other_table = app.data.(smpls{1}).AllCells(:, [channels(:)', others(:)']);
            % Tags are tags of gates in first sample (later intersect of tags till i-th sample)
            % Gate Table is a table of gates combined till i-th sample
            if ~isempty(vPD)
                waitbar(1 / numel(smpls), vPD, "Merging Samples 1/" + num2str(numel(smpls)));
            end
            
            % If there is only one sample
            if numel(smpls)==1
                % Pull the column names (gatetags) with the correct full path cell names
                tags = app.data.(smpls{1}).GateTags.Properties.VariableNames(...
                        ismember(app.data.(smpls{1}).GateTags{1, :}, paths));
                gate_table = app.data.(smpls{1}).AllCells(:, tags); 
                % Make Gate Tags
                gate_tags = table;
                for i=1:numel(paths)
                    gate_tags.(tags{i}) = [paths(i); short_names(i)];
                end
            else
                %% Build New Gate Tags
                paths_short_to_full = Helper.removed_path_to_full(app, paths_short, smpls{1});

                % size of the combined table
                nrows = size(channel_other_table, 1);
                gate_table = array2table(zeros(nrows, numel(paths_short)));
                                
                tags_new = cell(1, numel(paths_short));
                for gatei = 1:numel(paths_short)
                   tags_new{gatei} = [Constants.gate_tag num2str(gatei)];
                end
                gate_table.Properties.VariableNames = tags_new;
                
                % Pull the cells in this sample
                IND_Tags_i = ismember(app.data.(smpls{1}).GateTags{1, :}, paths_short_to_full(1, :));
                [IND_tags, IND_tagsb] = ismember(paths_short_to_full(1, :), app.data.(smpls{1}).GateTags{1, :});
                IND_tagsb = IND_tagsb(IND_tagsb ~= 0);

                tags_i = app.data.(smpls{1}).GateTags.Properties.VariableNames(IND_Tags_i);
                tags_i = tags_i(IND_tagsb);
                gate_table(:, IND_tags) = app.data.(smpls{1}).AllCells(:, tags_i);
                
                for s_i=2:numel(smpls)                
                    %% Build the Other channel merged table
                    % Extract current samples table, and save it's size, and merge with main table
                    channel_other_table_i = app.data.(smpls{s_i}).AllCells(:, [channels(:)', others(:)']);
                    channel_other_table = vertcat(channel_other_table, channel_other_table_i);
                
                    %% Build the Gate table
                    nrows = size(channel_other_table_i, 1);
                    gate_table_i = array2table(zeros(nrows, numel(paths_short)));
                    gate_table_i.Properties.VariableNames = tags_new;

                    IND_Tags_i = ismember(app.data.(smpls{s_i}).GateTags{1, :}, paths_short_to_full(1, :));
                    [IND_tags, IND_tagsb] = ismember(paths_short_to_full(1, :), app.data.(smpls{s_i}).GateTags{1, :});
                                        
                    IND_tagsb = IND_tagsb(IND_tagsb ~= 0);
                    tags_i = app.data.(smpls{s_i}).GateTags.Properties.VariableNames(IND_Tags_i);
                    tags_i = tags_i(IND_tagsb);
                    gate_table_i(:, IND_tags) = app.data.(smpls{s_i}).AllCells(:, tags_i);

                    % Add the sub table to the concatenated table
                    gate_table = vertcat(gate_table, gate_table_i);

                end
                %% Make Gate Tags
                gate_tags = table;
                for gatei=1:numel(paths_short)
                    gate_tags.(tags_new{gatei}) = [paths_short(gatei); short_names(gatei)];
                end
            end
        end

        function MFI_table = merge_MFI(app, smpls, MFI)
            % MERGE_MFI Allows to combine given MFI table from multiple
            % samples into one, which contains all common channels, gates
            % etc. It returns everything in one table, as it should be only
            % used as an additional function to the Helper.merge_AllCells.
            %
            % Input:
            %   - app - Instance of SAS3D
            %   - smpls - cell - Names of all samples to be combined into
            %       one. These have to exist under app.data, and will be
            %       left unchanged. Every sample has to contain that MFI.
            %   - MFI - char, string - Name of the MFI to be merged.
            %       Currently, either MFIRSN or MFICCN. It has to be
            %       contained in every sample.
            %
            % Output:
            %   - MFI_table - table - Table containing all the
            %       common channels, other, gates fields between all the
            %       samples for the union of the neighborhoods. Whether
            %       gates are common is determined on their path, and not
            %       gate tag. However gate tag corresponds to the tags from
            %       first sample given.
            %
            % Note:
            %   - A more correct approach would be to interpolate regions
            %   and join them together if they overlap. It would be more
            %   computationally intensive, and is not yet implemented due
            %   to lack of time.
            
            channels = Helper.get_channels(app, smpls);
            others = Helper.get_others(app, smpls, 'mfi', MFI);
            ignores = Helper.get_ignores(app, smpls, MFI);
            ignores = ignores(~ismember(ignores, channels));
            neighs = Helper.get_neighs(app, smpls, MFI);
            [paths, short_names] = Helper.get_gates(app, smpls);

            % Initialize Tables
            % Not all channels in AllCells will always be in MFI, i.e. distances etc.
            % Not all channels in MFIRSN will be in All Cells i.e. t-SNE
            Current_ChNames = app.data.(smpls{1}).(MFI).Properties.VariableNames;
            Target_Names = [channels(:)', others(:)', ignores(:)'];
            Target_Names = Target_Names(ismember(Target_Names, Current_ChNames));
            channel_other_table = app.data.(smpls{1}).(MFI)(:, Target_Names);

            % Neighborhood table
            if ~isempty(neighs)
                neigh_table = app.data.(smpls{1}).(MFI)(:, neighs);
            end

            % Tags are tags of gates in first sample (later intersect of tags till i-th sample)
            % Gate Table is a table of gates combined till i-th sample
            tags = app.data.(smpls{1}).GateTags.Properties.VariableNames(...
                ismember(app.data.(smpls{1}).GateTags{1, :}, paths));
            % New Tags might have been added after MFI was created
            idxs = ismember(tags, app.data.(smpls{1}).(MFI).Properties.VariableNames);
            paths = paths(idxs);
            tags = tags(idxs);

            gate_table = table2cell(app.data.(smpls{1}).(MFI)(:, tags));

            for s_i=2:numel(smpls)
                s = smpls{s_i};

                % Extract current samples table, and save it's size, and merge with main table
                Current_ChNames = app.data.(s).(MFI).Properties.VariableNames;
                Target_Names = [channels(:)', others(:)', ignores(:)'];
                Target_Names = Target_Names(ismember(Target_Names, Current_ChNames));
                channel_other_table_i = app.data.(s).(MFI)(:, Target_Names);
                channel_other_table = vertcat(channel_other_table, channel_other_table_i);

                % Merge with neighborhoods
                if ~isempty(neighs)
                    neigh_table = vertcat(neigh_table, app.data.(s).(MFI)(:, neighs));
                end

                [tags, tag_idxs] = ismember(app.data.(s).GateTags{1, :}, paths);
                tags = app.data.(s).GateTags.Properties.VariableNames(tags);  % Extract correct samples
                tag_idxs = tag_idxs(tag_idxs ~= 0);

                tags = tags(tag_idxs);  % Make ordering consistent across multiple samples

                % Convert to cell, as tags may be misaligned (i.e. tag1 in sample A may correspond to tag2 in sample B)
                gate_table_i = table2cell(app.data.(s).(MFI)(:, tags));
                gate_table = vertcat(gate_table, gate_table_i);
            end

            % Make Gate Tags
            if numel(smpls)==1
                % Deal with cases where MFI table is shorter than ALLCells
                INDPaths = ismember(app.data.(smpls{1}).GateTags{1,:}, paths);
                gate_tags = app.data.(smpls{1}).GateTags(:,INDPaths);
            else
                gate_tags = table;
                for i=1:numel(paths)
                    gate_tags.(strcat(Constants.gate_tag, num2str(i))) = [paths(i); short_names(i)];
                end
            end
                        
            % Combine Tables of channels, others and gates
            gate_table = cell2table(gate_table, 'VariableNames', gate_tags.Properties.VariableNames);

            if ~isempty(neighs)
                MFI_table = horzcat(channel_other_table, gate_table, neigh_table);
            else
                MFI_table = horzcat(channel_other_table, gate_table);
            end
            [MFI_table, ~] = unique(MFI_table, 'rows', 'stable');
            % Make sure NCells etc. in the MFI table are in the right order
            if ismember('NCells', MFI_table.Properties.VariableNames)
                MFI_table = movevars(MFI_table, 'NCells', 'After', 'Z');
            end
            if ismember('Effective_Neigh_Volume', MFI_table.Properties.VariableNames)
                MFI_table = movevars(MFI_table, 'Effective_Neigh_Volume', 'After', 'Z');
            end
            if ismember('Neigh_Volume', MFI_table.Properties.VariableNames)
                MFI_table = movevars(MFI_table, 'Neigh_Volume', 'After', 'Z');
            end
            if ismember('Effective_Neigh_Area', MFI_table.Properties.VariableNames)
                MFI_table = movevars(MFI_table, 'Effective_Neigh_Area', 'After', 'Z');
            end
            if ismember('Neigh_Area', MFI_table.Properties.VariableNames)
                MFI_table = movevars(MFI_table, 'Neigh_Area', 'After', 'Z');
            end
        end

        function smpl = make_empty_smpl()
            % MAKE_EMPTY_SMPL Creates an empty sample, which can be put under the app.data,
            % but contains no information whatsoever.
            
            smpl = struct;

            smpl.MetaData = struct;
            smpl.MetaData.Group = 1;
            smpl.MetaData.Annotation = 0;
            smpl.MetaData.Ver = Constants.CURR_VER;

            smpl.Surfaces = struct;
            smpl.Surfaces.UDS = struct;
            smpl.Surfaces.CCN = struct;
            smpl.Surfaces.RSN = struct;

            smpl.FNames = {};
            smpl.Path = "";

            smpl.AllCells = table;
            smpl.GateTags = table;
            smpl.tree = tree('All Cells', 'tag', strcat(Constants.gate_tag, '0'));
        end

        function paths = remove_smpl_from_paths(paths)
            % REMOVE_SMPL_FROM_PATHS Removes sample from path if such exists. It is the top-most
            % element in path, if it is common across all paths. It should
            % be used when comparing paths across different samples.
            %
            % Input:
            %   - paths - cell, char, string - Paths to be modified.
            %
            % Output:
            %   - paths - cell - Given paths, but with the top-most
            %       element of path (i.e. whatever is before first '/')
            %       removed. It is always a cell. If char/string where
            %       given then this cell contains 1 element.
            %
            % Note:
            %   - Due to Pre-August 2019 bug, this is necessary
            %       for compatibility reasons in Helper.merge_AllCells and
            %       Helper.merge_MFI. In particular paths would not be the
            %       same since tree root would have a different name, even
            %       though the rest of path was exactly the same.
            if isempty(paths)
                return;
            end
            if ~iscell(paths)
                paths = {paths};
            end
            first_elem = split(paths{1}, '/');
            first_elem = first_elem(1);

            % Ignore in house gates when removing filename
            firstignore = [Constants.gate_tag 'All'];
            logicignore = startsWith(paths, firstignore);
            
            % Ignore Radom points when removing filename
            firstignore = 'Rand_';
            logicignore = logicignore + startsWith(paths, firstignore);
            
            pathstmp = paths(~logicignore);

            if all(startsWith(pathstmp, first_elem))
                for p_idx=1:numel(pathstmp)
                    p = pathstmp{p_idx};
                    p = split(p, '/');
                    if numel(p) ~= 1
                        p = join(p(2:end), '/');
                    end
                    pathstmp(p_idx) = p;
                end
            end
            paths(~logicignore) = pathstmp;
%             paths = paths(2:end);
        end

        function full = removed_path_to_full(app, paths, smpls)
            % REMOVE_PATH_TO_FULL Adds sample to path. It will be placed in top-most element of
            % the path. Should be used when given paths were used to
            % compare paths across different samples, and now have to be
            % converted to form readable in smpls.
            %
            % Input:
            %   - paths - cell, char, string - Paths to be modified.
            %
            % Output:
            %   - full - cell - Cell of size:
            %       (number of samples, number of paths).
            %       Every element (x, y) is corresponding to a given path
            %       y, but with the top-most element added to correspond to
            %       the top of the tree of sample x. If there is no
            %       correspondance to path y, then it's kept in form it was
            %       passed in.
            %
            % Note:
            %   - Due to Pre-August 2019 bug, this is necessary
            %       for compatibility reasons in Helper.merge_AllCells and
            %       Helper.merge_MFI. In particular paths would not be the
            %       same since tree root would have a different name, even
            %       though the rest of path was exactly the same.
            
            if ~iscell(paths)
                paths = {paths};
            end
            if ~iscell(smpls)
                smpls = {smpls};
            end
            
            full = cell(numel(smpls), numel(paths));
            for s_i = 1:numel(smpls)
                s = smpls{s_i};
                full_paths = app.data.(s).GateTags{1, :};
                for p_i = 1:numel(paths)
                    p = paths{p_i};
                    ind = endsWith(full_paths, p);
                    if any(ind)
                        full{s_i, p_i} = full_paths{ind};
                    else
                        full{s_i, p_i} = p;
                    end
                end
            end
        end

        function result = assert_valid_net(net_in)
            % ASSERT_VALID_NET Asserts validity of net, with generic requirements for model.
            % i.e. checks whether modelcan be assigned to:
            %   app.net.(net_in.name) = net_in
            
            result = true;
            if ~isstruct(net_in)
                result = false;
                return;
            end
            if sum(~ismember(Constants.model_required_fields, fieldnames(net_in))) ~= 0
                result = false;
                return;
            end
            if strcmp(net_in.type, 'SOM')
                if ~isfield(net_in, 'Network') || ~isa(net_in.Network, 'network')
                    result = false;
                    return;
                end
            end
        end

        function cont = add_model(app, dat, m)
            % ADD_MODEL
            % !!! WORKING ONLY WITH Post-Feb 2019 VERSION OF MODELS. !!!
            %
            % Adds a model to the app. 
            % Makes interactive questions with the user,
            % if the model with the same name was found in the app.
            
            cont = true;
            new_m = m;
            if isfield(m, app.net)
                message = strcat('Model ', Helper.full_var(m), ' was found in current SAS3D instance. Do you want to override it with new model?');
                title = 'Duplicate Model found';
                answer = questdlg(message, title, ...
                    'Yes', 'No', 'Rename model', ...
                    'defbtn', 'Rename model');
                switch answer
                    case ''
                        cont = false;
                        return;
                    case 'Yes'
                        % In such case we just continue with normal loading of data.
                    case 'No'
                        return;
                    case 'Rename model'
                        prompt = {'Enter name of new model:'};
                        title = 'Rename model';
                        dims = [1 35];
                        definput = {strcat('Model', num2str(numel(fieldnames(app.net)) + 1))};
                        new_m = Helper.valid_var(inputdlg(prompt,title,dims,definput));
                        dat.(m).name = new_m;
                end
            end
            if Helper.assert_valid_net(dat.(m))
                app.net.(new_m) = dat.(m);
            end
        end

        function boolean = setequal(A, B)
            % SETEQUAL Compares two sets. Permutation invariant.

            if ~iscell(A)
                A = {A};
            end
            if ~iscell(B)
                B = {B};
            end
            % setdiff requires cell arrays of char vectors, NOT STRINGS!
            char_A = cell(size(A));
            char_B = cell(size(B));
            for i=1:numel(A)
                tmp = A{i};
                if islogical(tmp)
                    tmp = 1.0 * tmp;
                end
                if isnumeric(tmp)
                    tmp = num2str(tmp);
                end
                char_A{i} = char(tmp);
            end
            for i=1:numel(B)
                tmp = B{i};
                if islogical(tmp)
                    tmp = 1.0 * tmp;
                end
                if isnumeric(tmp)
                    tmp = num2str(tmp);
                end
                char_B{i} = char(tmp);
            end
            boolean = isempty(setdiff(char_A, char_B)) && isempty(setdiff(char_B, char_A));
        end

        function celled = logical2cell(A)
            % LOGICAL2CELL Converts logical array to cell array.
            
            celled = cell(size(A));
            for i=1:numel(celled)
                celled(i) = {logical(A(i))};
            end
        end

        function func_dataWDW(app, web)
            % FUNC_DATEWDW Front-End of Show Data Table function, and all Merge/Remove
            % Functions that go along with it.
            %
            % Input:
            %   - app - Instance of SAS3D
            %
            % Modifies:
            %   - app - Calling it itself doesn't change anything. However,
            %       user gains access to functions that tha can add,
            %       remove, and add populations or samples.
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};

            if ~Helper.any_sample(app)
                return;
            end
            smplstmp = fieldnames(app.data);
            
            % Build the clustering options menus
            UIfig = uifigure('Name', 'Phenotype Table', 'Scrollable', 'on');
            if web==1
                UIfig.Visible='OFF';
            end
            UIfig.Position = alpha*[10 10 1000 800];
            
            % Create the table of options, and populate a stub
            t = uitable(UIfig);
            
            tmp = Helper.populate_table(app, 'smpls', smplstmp{1});

            % Calculate the sums and percentages.
            tags = Helper.gate_full2tag(app, tmp(:, 3), smplstmp{1});
            tags = tags(~cellfun('isempty', tags));
            sums = sum(table2array(app.data.(smplstmp{1}).AllCells(:, tags)))';
            parent_sums = zeros(size(sums));
            parents = Helper.gate_tag2parent_tag(app, tags);
            for p=1:numel(parents)
                if strcmp(parents(p), strcat(Constants.gate_tag, '0'))
                    parent_sums(p) = size(app.data.(smplstmp{1}).AllCells, 1);
                else
                    try
                        parent_sums(p) = sum(table2array(app.data.(smplstmp{1}).AllCells(:, parents(p))));
                    catch
                        parent_sums(p) = NaN;
                    end
                end
            end
            percents = sums ./ parent_sums .* 100;

            t.Data = cell(size(tmp, 1)+1, 6);
            t.Position = alpha*[0 90 1000 710];
            t.ColumnName = {'', 'Phenotype', 'Total', 'Percent of Parent', '', 'Sample'};
            t.Data(1:numel(sums), 1) = tmp(1:numel(sums), 2);
            for idx =1:numel(t.Data(1:numel(sums), 1))
                if ~isempty(t.Data(idx, 1))
                    t.Data(idx, 1) = {false};
                end
            end
            % Deal with empty elements
            indempty = ~cellfun(@isempty,tmp(:, 3));
            t.Data(1:numel(sums), 2) = tmp(indempty, 3);
            t.Data(1:numel(sums), 3) = num2cell(sums);
            t.Data(1:numel(percents), 4) = num2cell(percents);
            t.Data(1:numel(smplstmp), 5) = tmp(1:numel(smplstmp), 4);
            t.Data(1:numel(smplstmp), 6) = tmp(1:numel(smplstmp), 5);
            
            % Add a select all button
            t.Data(end, :) = {false, 'Select All', [], [],false, 'Select All'};
            
            t.ColumnEditable = [true true false false true false];
            % t.ColumnWidth = {50, 100, 350, 50, 100, 125, 100, 125};
            t.CellEditCallback = @(dd, p) edited_table(app, dd, p);

            % Remove Population button
            btn = uibutton(UIfig, 'push', 'ButtonPushedFcn', @(btn,event) remove_PH_wrap(app, t));
            btn.Position = alpha*[10, 10, 117.5, 70];
            btn.Text = 'Remove Population';
            Helper.func_SetCLR(app, btn, 'button')

            % Merge populations button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) merge_PH_wrap(app, t));
            btn.Position = alpha*[137.5, 10, 117.5, 70];
            btn.BackgroundColor = app.GUIOPTS.bgclr;
            btn.Text = 'Merge Populations';
            Helper.func_SetCLR(app, btn, 'button')

            % Add Population button
            btn = uibutton(UIfig, 'push', 'ButtonPushedFcn', @(btn,event) add_PH_wrap(app, t));
            btn.Position = alpha*[265, 10, 117.5, 70];
            btn.Text = 'Add Population';
            Helper.func_SetCLR(app, btn, 'button')
            
            % Add Population button
            btn = uibutton(UIfig, 'push', 'ButtonPushedFcn', @(btn,event) add_SM_wrap(app, t));
            btn.Position = alpha*[265+130.5, 10, 117.5, 70];
            btn.Text = 'New Sample';
            Helper.func_SetCLR(app, btn, 'button')

            % Remove Sample button
            btn = uibutton(UIfig, 'push', 'ButtonPushedFcn', @(btn,event) remove_smpl_wrap(app, t));
            btn.Position = alpha*[t.Position(3) - 117.5 - 10, 10, 117.5, 70];
            btn.Text = 'Remove Sample';
            Helper.func_SetCLR(app, btn, 'button')

            % Merge Samples button
            btn = uibutton(UIfig, 'push', 'ButtonPushedFcn', @(btn,event) merge_smpl_wrap(app, t));
            btn.Position = alpha*[t.Position(3) - 117.5 - 137.5, 10, 117.5, 70];
            btn.Text = 'Merge Samples';
            Helper.func_SetCLR(app, btn, 'button')

            function add_PH_wrap(app, dd)
                ind = dd.Data(:, 5);
                ind = ind(~cellfun('isempty', ind));
                ind = cell2mat(ind);
                smpls = dd.Data(ind, 6);
                IO.func_AddPH(app, smpls);
                update_table(app, dd);
            end

            function remove_PH_wrap(app, dd)
                ind = dd.Data(:, 5);
                ind = ind(~cellfun('isempty', ind));
                ind = cell2mat(ind);
                smpls = dd.Data(ind, 6);
                
                ind = dd.Data(:, 1);
                ind = ind(~cellfun('isempty', ind));
                ind = cell2mat(ind);
                phns = dd.Data(ind, 2);
                
                IO.func_RemovePH(app, smpls, phns);
                update_table(app, dd);
            end

            function remove_smpl_wrap(app, dd)
                ind = dd.Data(:, 5);
                ind = ind(~cellfun('isempty', ind));
                ind = cell2mat(ind);
                smpls = dd.Data(ind, 6);
                IO.func_RemoveSmpl(app, smpls);
                update_table(app, dd);
            end

            function merge_PH_wrap(app, dd)
                ind = dd.Data(:, 5);
                ind = ind(~cellfun('isempty', ind));
                ind = cell2mat(ind);
                smpls = dd.Data(ind, 6);
                ind = dd.Data(:, 1);
                ind = ind(~cellfun('isempty', ind));
                ind = cell2mat(ind);
                phns = dd.Data(ind, 2);
                IO.func_MergePH(app, smpls, phns);
                update_table(app, dd);
            end

            function add_SM_wrap(app, dd)
                % Pull the indeces of the samples you want to merge
                ind = dd.Data(:, 5);
                ind = ind(~cellfun('isempty', ind));
                ind = cell2mat(ind);
                % Pull the sample names of the samples you want to merge
                smpls = dd.Data(ind, 6);
                
                
                % Pull the names of the phenotypes you want to include in the new samples
                ind = dd.Data(:, 1);
                ind = ind(~cellfun('isempty', ind));
                ind = cell2mat(ind);
                phnsKeep = dd.Data(ind, 2);
                
                % Ignore Select All
                IND = contains(phnsKeep, 'Select All');
                phnsKeep = phnsKeep(~IND);
                %Include cells from
                [ind,tf] = listdlg(...
                                   'PromptString',...
                                   {'All Gates will be in new sample',...
                                   'Cells not exclusive to selected', ...
                                   'will be excluded'},...
                                   'ListString', phnsKeep...
                                   );
                if tf==0
                    return
                end
                phnsKeep = phnsKeep(ind);
                                
                % Pull the names of the phenotypes you don't want to include in the new samples
                ind = dd.Data(:, 1);
                ind = ind(~cellfun('isempty', ind));
                ind = cell2mat(ind);
                phns = dd.Data(~ind, 2);
                % Don't remove All Cells Gate
                IND = contains(phns, 'All/All');
                phns = phns(~IND);
                % Ignore Select All
                IND = contains(phns, 'Select All');
                phns = phns(~IND);

                new_name = inputdlg("Enter name of new sample:");
                if isempty(new_name)
                    return
                end
                
                % Merge the two samples
                IO.func_NewSmpls(app, smpls, phns, Helper.valid_var(new_name), phnsKeep);
                update_table(app, dd);
            end
            
            function merge_smpl_wrap(app, dd)
                % Pull the indeces of the samples you want to merge
                ind = dd.Data(:, 5);
                ind = ind(~cellfun('isempty', ind));
                ind = cell2mat(ind);
                % Pull the sample names of the samples you want to merge
                smpls = dd.Data(ind, 6);
                new_name = inputdlg("Enter name of new sample:");
                % Merge the two samples
                IO.func_MergeSmpls(app, smpls, Helper.valid_var(new_name), 1);
                update_table(app, dd);
            end

            function update_table(app, dd, p)
                if nargin<3
                    p = struct;
                    p.Indices = [1,5];
                    p.NewData = 0;
                end
                
                % Process 'Select All' button
                if p.Indices(1) == find(strcmp(dd.Data(:, 2), 'Select All')) && p.Indices(2) == 1
                    ind = ~cellfun('isempty', dd.Data(:, 2));
                    dd.Data(ind, 1) = {logical(p.NewData)};  
                    return
                    
                elseif p.Indices(1) == find(strcmp(dd.Data(:, 6), 'Select All')) && p.Indices(2) == 5
                    ind = ~cellfun('isempty', dd.Data(:, 6));
                    dd.Data(ind, 5) = {logical(p.NewData)};
                    % Make sure that at least one thing is selected
                    if ~logical(p.NewData)
                        dd.Data{1, 5} = true;
                    end
                end
                
                if p.Indices(2) == 5    
                    % Remove Select All
                    tdat = dd.Data(1:(end-1), :);

                    ind = tdat(:, 5);
                    ind = ind(~cellfun('isempty', ind));
                    ind = cell2mat(ind);
                    smpls = tdat(ind, 6);
                    smpls = smpls(isfield(app.data, smpls));
                    if isempty(smpls)
                        if isempty(fieldnames(app.data))
                            % Last sample removed
                            return;
                        else
                            smpls = fieldnames(app.data);
                            smpls = smpls(1);
                        end
                    end

                    tmps = Helper.populate_table(app, 'smpls', smpls, 'fill_checkbox', false);

                    if ~isempty(tmps) % Sanity check.
                        % Figure out sums, percentages.
                        first_tags_tmp = Helper.gate_full2tag(app, tmps(:, 3), smpls{1});
                        first_parents_tmp = Helper.gate_tag2parent_tag(app, first_tags_tmp, smpls{1});

                        first_tags_tmp = first_tags_tmp(~cellfun('isempty', first_tags_tmp));
                        first_parents_tmp = first_parents_tmp(~cellfun('isempty', first_parents_tmp));

                        sums_tmp = sum(table2array(app.data.(smpls{1}).AllCells(:, first_tags_tmp)))';
                        parent_sums_tmp = zeros(size(sums_tmp));
                        for fp=1:numel(first_parents_tmp)
                            if strcmp(first_parents_tmp(fp), strcat(Constants.gate_tag, '0'))
                                parent_sums_tmp(fp) = size(app.data.(smpls{1}).AllCells, 1);
                            else
                                parent_sums_tmp(fp) = sum(table2array(app.data.(smpls{1}).AllCells(:, first_parents_tmp(fp))));
                            end
                        end

                        % If more than one sum is chosen, combine results.
                        if numel(smpls) > 1
                            for idx_smpl=2:numel(smpls)
                                s = smpls{idx_smpl};
                                curr_tags = Helper.gate_full2tag(app, tmps(:, 3), s);
                                curr_parent_tags = Helper.gate_tag2parent_tag(app, curr_tags, s);

                                curr_tags = curr_tags(~cellfun('isempty', curr_tags));
                                curr_parent_tags = curr_parent_tags(~cellfun('isempty', curr_parent_tags));

                                sums_tmp = sums_tmp + sum(table2array(app.data.(s).AllCells(:, curr_tags)))';

                                for fp=1:numel(curr_parent_tags)
                                    if strcmp(curr_parent_tags(fp), strcat(Constants.gate_tag, '0'))
                                        parent_sums_tmp(fp) = parent_sums_tmp(fp) + ...
                                            size(app.data.(s).AllCells, 1);
                                    else
                                        parent_sums_tmp(fp) = parent_sums_tmp(fp) + ...
                                            sum(table2array(app.data.(s).AllCells(:, curr_parent_tags(fp))));
                                    end
                                end
                            end
                        end

                        % TODO: Add columns which will specify distributions per sample.
                        percents_tmp = sums_tmp ./ parent_sums_tmp .* 100.0;
                        tdat = cell(size(tmps, 1), 6);
                        tdat(:, 1) = tmps(:, 2);
                        tdat(:, 2) = tmps(:, 3);
                        tdat(1:numel(sums_tmp), 3) = num2cell(sums_tmp);
                        tdat(1:numel(percents_tmp), 4) = num2cell(percents_tmp);
                        tdat(:, 5) = tmps(:, 4);
                        tdat(:, 6) = tmps(:, 5);
                        
                        dd.Data = [tdat; dd.Data(end, :)];
                    end
                end
            end

            function edited_table(app, dd, p)
                if p.Indices(2) == 5 || p.Indices(2) == 1
                    update_table(app, dd, p);
                elseif p.Indices(2) == 2

                    if strcmp(p.PreviousData, 'Select All')
                        dd.Data{p.Indices(1), p.Indices(2)} = 'Select All';
                        return
                    end
                    % Check if parent was not changed.
                    prev_parent = split(p.PreviousData, '/');
                    prev_parent = prev_parent(1);
                    curr_parent = split(p.EditData, '/');
                    if numel(curr_parent) <= 1 || ~strcmp(curr_parent(1), prev_parent)
                        dd.Data(p.Indices(1), p.Indices(2)) = {p.PreviousData};
                        warndlg('Sorry, but the name of parent (before /) cannot change.');
                        return;
                    end

                    % Check if already is a phenotype in one of the samples.
                    ind = dd.Data(:, 5);
                    ind = ind(~cellfun('isempty', ind));
                    ind = cell2mat(ind);
                    smpls = dd.Data(ind, 6);
                    insmpl = false(numel(smpls), 1);
                    for s=1:numel(smpls)
                        insmpl(s) = ismember(p.EditData, table2cell(app.data.(smpls{s}).GateTags(2, :)));
                    end
                    if any(insmpl)
                        dd.Data(p.Indices(1), p.Indices(2)) = {p.PreviousData};
                        warndlg(strcat('Sorry, but the new name must be unique. Such Phenotype already exists in sample(s):  ', ...
                                strjoin(smpls(insmpl), '  ')));
                        return;
                    end
                    
                    IO.func_Rename(app, smpls, p.PreviousData, p.EditData);
                    
                    update_table(app, dd, p);
                end
            end
        end

        function func_SetCLR(app, object, type)
            % FUNC_SETCLR Set the colors for given object
            % Sets color of the object to a currently chosen one in
            % SAS3D (app). type specifies what this object is.

            switch type
                case 'UIfigure'
                    object.Color = app.GUIOPTS.bgclr;
                case 'figure'
                    object.Color = app.GUIOPTS.bgclr;
                    object.InvertHardcopy = 'off';
                    ax = gca;
                    ax.XColor = app.GUIOPTS.txtclr;
                    ax.YColor = app.GUIOPTS.txtclr;
                    ax.ZColor = app.GUIOPTS.txtclr;
                    ax.FontName = app.GUIOPTS.FontName;
                    ax.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
                case 'label'
                    object.FontColor = app.GUIOPTS.txtclr;
                    object.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
                    object.FontName = app.GUIOPTS.FontName;
                    object.FontWeight = app.GUIOPTS.FontWeight;
                case 'labelSpecial'
                    object.FontName = app.GUIOPTS.FontName;
                    if ismember('FontColor', fieldnames(object))
                        object.FontColor = app.GUIOPTS.txtclr;
                    elseif ismember('ForegroundColor', fieldnames(object))
                        object.ForegroundColor = app.GUIOPTS.txtclr; 
                    end
                    object.BackgroundColor = app.GUIOPTS.bgclr;
                    % make an exception for links
                    if size(object.String, 1)==1
                        if any(contains(object.String, {'CytoMAP Wiki', 'Gerner Lab'}))
                            object.ForegroundColor = 'b';
                        elseif any(contains(object.String, {'SAS3D'}))
                            object.FontName = 'Arial Narrow';
                        end
                    end
                case 'UIField'
                    object.BackgroundColor = app.GUIOPTS.bgclr;
                    object.FontColor = app.GUIOPTS.txtclr;
                    object.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
                    object.FontName = app.GUIOPTS.FontName;
                case 'UICpopup'
                    object.BackgroundColor = app.GUIOPTS.bgclr;
                    object.ForegroundColor = app.GUIOPTS.txtclr;
                    object.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI}*0.7;
                    object.FontName = app.GUIOPTS.FontName;
                case 'UICpushbutton'
                    object.BackgroundColor = app.GUIOPTS.btnbgclr;
                    object.ForegroundColor = app.GUIOPTS.txtclr;
                    if strcmp(object.String, 'New Figure')
                        object.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
                    else
                        object.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI}*0.8;
                    end
                    object.FontName = app.GUIOPTS.FontName;
                case 'button' 
                    object.BackgroundColor = app.GUIOPTS.btnbgclr;
                    if ismember(fieldnames(object), 'FontColor')
                        object.FontColor = app.GUIOPTS.txtclr;
                    elseif ismember(fieldnames(object), 'ForegroundColor')
                        object.ForegroundColor = app.GUIOPTS.txtclr; 
                    end
                    object.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
                    object.FontName = app.GUIOPTS.FontName;
%                     object.Position(4) = app.GUIOPTS.btnheight;
                case 'table'
%                     object.BackgroundColor = app.GUIOPTS.bgclr;
%                     object.ForegroundColor = app.GUIOPTS.txtclr;
                    object.FontSize = app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI};
                    object.FontName = app.GUIOPTS.FontName;
%                     object.RowStriping = 'on';
            end
        end

        function func_EditGUIInterface(app, web)
            % FUNC_EDITGUIINTERFACE Front-End of Edit Main GUI Configuration.
            % It lets the user change the fonts colors etc.
            %
            % Input:
            %   - app - Instance of SAS3D
            %
            % Modifies:
            %   - app - Way the SAS3D looks. For example switches colors
            %   of background, text etc.
            
            if nargin<2
                web=0;
            end
            alpha=app.GUIOPTS.Alpha{app.GUIOPTS.AlphaI};
            
            %% Make a interface
            UIfig = uifigure('Name', 'User Interface Options', 'Resize', 'on');
            if web==1
                UIfig.Visible='OFF';
            end
            UIfig.Position = alpha*[10 10 300 400];
            Helper.func_SetCLR(app, UIfig, 'UIfigure')

            % Edit Primary color
            lbl = uilabel(UIfig); lbl.Text = sprintf('Primary Background Color');
            lbl.Position = alpha*[10 UIfig.Position(4)-15 400 15];
            Helper.func_SetCLR(app, lbl, 'label')
            edt = uieditfield(UIfig,'text');
            edt.Value = '';
            edt.Position = alpha*[10+175 UIfig.Position(4)-20 50 20];
            edt.Editable = 'off';
            Helper.func_SetCLR(app, edt, 'UIField')
            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Helper.EditGUIInterface_backend(app, 'bgclr', app.GUIOPTS.bgclr,UIfig));
            btn.Position = alpha*[10+225, UIfig.Position(4)-20, 50, 20];
            btn.Text = 'Edit';
            Helper.func_SetCLR(app, btn, 'button')

            % Edit Secondary color
            lbl = uilabel(UIfig); lbl.Text = sprintf('Secondary Background Color');
            lbl.Position = alpha*[10 UIfig.Position(4)-35 400 15];
            Helper.func_SetCLR(app, lbl, 'label')
            edt = uieditfield(UIfig,'text');
            edt.Value = '';
            edt.Position = alpha*[10+175 UIfig.Position(4)-40 50 20];
            edt.Editable = 'off';
            Helper.func_SetCLR(app, edt, 'UIField')
            edt.BackgroundColor = app.GUIOPTS.btnbgclr;
            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Helper.EditGUIInterface_backend(app, 'btnbgclr', app.GUIOPTS.btnbgclr,UIfig));
            btn.Position = alpha*[10+225, UIfig.Position(4)-40, 50, 20];
            btn.Text = 'Edit';
            Helper.func_SetCLR(app, btn, 'button')

            % Edit Secondary color
            lbl = uilabel(UIfig); lbl.Text = sprintf('Secondary Background Color');
            lbl.Position = alpha*[10 UIfig.Position(4)-35 400 15];
            Helper.func_SetCLR(app, lbl, 'label')
            edt = uieditfield(UIfig,'text');
            edt.Value = '';
            edt.Position = alpha*[10+175 UIfig.Position(4)-40 50 20];
            edt.Editable = 'off';
            Helper.func_SetCLR(app, edt, 'UIField')
            edt.BackgroundColor = app.GUIOPTS.btnbgclr;
            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Helper.EditGUIInterface_backend(app, 'btnbgclr', app.GUIOPTS.btnbgclr,UIfig));
            btn.Position = alpha*[10+225, UIfig.Position(4)-40, 50, 20];
            btn.Text = 'Edit';
            Helper.func_SetCLR(app, btn, 'button')

            % Edit Font color
            lbl = uilabel(UIfig); lbl.Text = sprintf('Font Color');
            lbl.Position = alpha*[10 UIfig.Position(4)-55 400 15];
            Helper.func_SetCLR(app, lbl, 'label')
            edt = uieditfield(UIfig,'text');
            edt.Value = '';
            edt.Position = alpha*[10+175 UIfig.Position(4)-60 50 20];
            edt.Editable = 'off';
            Helper.func_SetCLR(app, edt, 'UIField')
            edt.BackgroundColor = app.GUIOPTS.txtclr;
            
            % Create a push button
            btn = uibutton(UIfig,'push', 'ButtonPushedFcn', @(btn,event) Helper.EditGUIInterface_backend(app, 'txtclr', app.GUIOPTS.txtclr,UIfig));
            btn.Position = alpha*[10+225, UIfig.Position(4)-60, 50, 20];
            btn.Text = 'Edit';
            Helper.func_SetCLR(app, btn, 'button')

            % Edit Font type
            lbl = uilabel(UIfig); lbl.Text = sprintf('Font Type');
            lbl.Position = alpha*[10 UIfig.Position(4)-75 400 15];
            Helper.func_SetCLR(app, lbl, 'label')
            dropdown = uidropdown(UIfig);
            dropdown.Items = app.GUIOPTS.Fonts;
            dropdown.Value = dropdown.Items{strcmp(dropdown.Items, app.GUIOPTS.FontName)};
            dropdown.Position = alpha*[10+175, UIfig.Position(4)-80, 100, 20];
            dropdown.ValueChangedFcn = @(~, ~) Helper.EditGUIInterface_backend(app, 'FontName', dropdown.Value, UIfig);
            Helper.func_SetCLR(app, dropdown, 'UIField')

            % Edit Font Size
            lbl = uilabel(UIfig); lbl.Text = sprintf('Font Size');
            lbl.Position = alpha*[10 UIfig.Position(4)-95 400 15];
            Helper.func_SetCLR(app, lbl, 'label')
            dropdown2 = uidropdown(UIfig);
            dropdown2.Items = app.GUIOPTS.AlphaL;
            dropdown2.Value = app.GUIOPTS.AlphaL{app.GUIOPTS.AlphaI};
            dropdown2.Position = alpha*[10+175 UIfig.Position(4)-100 100 20];
            dropdown2.ValueChangedFcn = @(~, ~) Helper.EditGUIInterface_backend(app, 'txtsz', dropdown2.Value, UIfig);
            Helper.func_SetCLR(app, dropdown2, 'UIField')

            % Edit button height
            lbl = uilabel(UIfig); lbl.Text = sprintf('Button Height');
            lbl.Position = alpha*[10 UIfig.Position(4)-115 400 15];
            Helper.func_SetCLR(app, lbl, 'label')
            edt = uieditfield(UIfig,'numeric');
            edt.Value = app.GUIOPTS.btnheight;
            edt.Position = alpha*[10+175 UIfig.Position(4)-120 100 20];
            edt.Editable = 'on';
            edt.ValueChangedFcn = @(~, ~) Helper.EditGUIInterface_backend(app, 'btnheight', edt.Value, UIfig);
            Helper.func_SetCLR(app, edt, 'UIField')

            % Default Colorschemes
            lbl = uilabel(UIfig); lbl.Text = sprintf('Defaults');
            lbl.Position = alpha*[10 UIfig.Position(4)-135 400 15];
            Helper.func_SetCLR(app, lbl, 'label')
            dropdown = uidropdown(UIfig);
            dropdown.Items = {'Light', 'Dark', 'Party', 'Random'};
            dropdown.Value = app.GUIOPTS.Scheme;
            dropdown.Position = alpha*[10+175, UIfig.Position(4)-140, 100, 20];
            dropdown.ValueChangedFcn = @(~, ~) Helper.EditGUIInterface_backend(app, dropdown.Value, dropdown.Value, UIfig);
            Helper.func_SetCLR(app, dropdown, 'UIField')
            
            % Change Main UI
            Helper.func_SetCLR(app, app.GUI.Main, 'UIfigure')
            fieldnms = fieldnames(app.GUI.Labels);
            for obji = 1:numel(fieldnms)
                Helper.func_SetCLR(app, app.GUI.Labels.(fieldnms{obji}), 'labelSpecial')
            end

            fieldnms = fieldnames(app.GUI.Buttons);
            for obji = 1:numel(fieldnms)
                Helper.func_SetCLR(app, app.GUI.Buttons.(fieldnms{obji}), 'UICpushbutton')
            end

            fieldnms = fieldnames(app.GUI.Dropdowns);
            for obji = 1:numel(fieldnms)
                Helper.func_SetCLR(app, app.GUI.Dropdowns.(fieldnms{obji}), 'UIField')
            end

        end

        function EditGUIInterface_backend(app, field, option, UIfig)
            % EDITGUIINTERFACE_BACKEND Back-End of Edit Main GUI Configuration.
            % Changes the fonts colors etc.
            %
            % Input:
            %   - app - Instance of SAS3D
            %   - field - What to change. Can either be a specific thing to
            %       change such as background, text, font, or a theme such
            %       as light/dark etc.
            %   - option - If field specified what to change, and not
            %       theme, specifies in what way to change this thing.
            %       i.e. color name, font size etc.
            %   - UIfig - Handle to interface front-end, which will be
            %       closed after this function executes.
            %
            % Modifies:
            %   - app - Way the SAS3D looks. For example switches colors
            %   of background, font size, overall theme etc.
            switch field
                case 'bgclr'
                    c = uisetcolor(option);
                    app.GUIOPTS.bgclr = c;
                case 'btnbgclr'
                    c = uisetcolor(option);
                    app.GUIOPTS.btnbgclr = c;
                case 'txtclr'
                    c = uisetcolor(option);
                    app.GUIOPTS.txtclr = c;
                case 'FontName'
                    app.GUIOPTS.FontName = option;
                case 'txtsz'
                    index = find(strcmp(app.GUIOPTS.AlphaL,option));
                    app.GUIOPTS.AlphaI = index;
                case 'btnheight'
                    app.GUIOPTS.btnheight = option;
                case 'Light'
                    % Color Scheme
                    app.GUIOPTS.Scheme = 'Light';
                    % background window color
                    app.GUIOPTS.bgclr = [1,1,1];
                    % button background color
                    app.GUIOPTS.btnbgclr = [1,1,1];
                    % text color
                    app.GUIOPTS.txtclr = [0,0,0];
                    % text font
                    app.GUIOPTS.FontName = 'Helvetica'; % Use listfonts command to see what fonts are available
                    % text size
                    app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI} = 11;
                    % button height
                    app.GUIOPTS.btnheight = 20;
                    % Font weight
                    app.GUIOPTS.FontWeight = 'normal';
                case 'Dark'
                    % Color Scheme
                    app.GUIOPTS.Scheme = 'Dark';
                    % background window color
                    app.GUIOPTS.bgclr = [0,0,0];
                    % button background color
                    app.GUIOPTS.btnbgclr = [0,0,0];
                    % text color
                    app.GUIOPTS.txtclr = [1,1,1];
                    % text font
                    app.GUIOPTS.FontName = 'Helvetica';
                    % text size
                    app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI} = 11;
                    % button height
                    app.GUIOPTS.btnheight = 20;
                    % Font weight
                    app.GUIOPTS.FontWeight = 'normal';
                case 'Party'
                    % Color Scheme
                    app.GUIOPTS.Scheme = 'Party';
                    % background window color
                    app.GUIOPTS.bgclr = [0,1,0];
                    % button background color
                    app.GUIOPTS.btnbgclr = [1,0,0];
                    % text color
                    app.GUIOPTS.txtclr = [0,0,1];
                    % text font
                    app.GUIOPTS.FontName = 'Harrington';
                    % text size
                    app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI} = 11;
                    % button height
                    app.GUIOPTS.btnheight = 20;
                    % Font weight
                    app.GUIOPTS.FontWeight = 'normal';
                case 'Random'
                    % Color Scheme
                    app.GUIOPTS.Scheme = 'Random';
                    % background window color
                    app.GUIOPTS.bgclr = rand(1,3);
                    % button background color
                    app.GUIOPTS.btnbgclr = rand(1,3);
                    % text color
                    app.GUIOPTS.txtclr = rand(1,3);
                    % text font
                    app.GUIOPTS.FontName = app.GUIOPTS.Fonts{round(numel(app.GUIOPTS.Fonts).*rand)};
                    % text size
                    app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI} = 11;
                    % button height
                    app.GUIOPTS.btnheight = 20;
                    % Font weight
                    app.GUIOPTS.FontWeight = 'normal';
            end
            if ~isempty(UIfig)
                Helper.func_EditGUIInterface(app)
                close(UIfig)
            end
        end

        function [Dat_Pre, Dat_PreALL, INDDatPre, INDzrs, INDons] = func_loaddata(app, smplnms, phenotypes, MFIList, pheno_weights, MFI_weights, DataType, DataPrep, NormOpt, RmvZeros, DataPrepMFI)
            % FUNC_LOADDATA Loads up and/or prepares data in a user
            % defined way.
            %
            % Input:
            %   - app - Instance of SAS3D
            %   - smpls - cell - Names of all samples to be loaded.
            %       These have to exist under app.data, and will be
            %       left unchanged. Every sample has to contain that MFI.
            %   - phenotypes - cell - Names of all phenotypes to be combined into
            %       one. These have to exist under app.data, and will be
            %       left unchanged. Every sample has to contain that phenotype.
            %   - MFIList - cell - Names of all channels to be combined into
            %       one. These have to exist under app.data, and will be
            %       left unchanged. Every sample has to contain that channel.
            %
            %   - pheno_weights - cell - weights for the selected cell
            %       types
            %   - MFI_weights - cell - weights for the selected channels
            %       types
            %       left unchanged. Every sample has to contain that channel.
            %   - DataType - string - Type of neighborhoods being loaded
            %   - NormOpt - string - Option to either normalize data per 
            %       sample or per whole dataset 
            %   - Dat_Pre - double - Data which is supposed to be
            %       processed, and prepared.
            %   - DataPrep - string - What techinque to use when modifiying
            %       the data.
            %   - RmvZeros - logical - option to either remove empty
            %       neighborhoods or leave them in the dataset
            %
            % Output:
            %   - Dat_Pre - double - Data which is supposed to be
            %       processed, and prepared.
            %   - Dat_PreALL - table - Indexed portion of whole table which Dat_Pre came from.
            %       It will be used to define variables which will be used
            %       in data preparation process.
            
            %% Load in the first dataset
            if ~iscell(smplnms)
                smplnms = {smplnms};
            end
            % If I didn't pass the load data function a MFI norm option
            if nargin < 11
                DataPrepMFI = 'Standardize: subtract Mean, divide by standard deviation';
            else
                switch DataPrepMFI
                    % Pull the raw data and don't normalize it
                    case 'Sum MFI per neighborhood'
                        DataPrepMFI = 'Cellularity: Number of Cells / Neighborhood';
                    case 'Raw MFI per cell'
                        DataPrepMFI = 'Cellularity: Number of Cells / Neighborhood';
                    % normalize to the total number of cells in neighborhood
                    case 'Average MFI per cell per neighborhood'
                        DataPrepMFI = 'Composition: Number of Cells / Total cells in Neighborhood';
                    % Normalize to the mean MFI in the sample
                    case 'MFI normalized to mean MFI per neighborhood'
                        DataPrepMFI = 'Global Fold Change: Number of Cells / Mean Cells in Tissue Neighborhoods';
                    case 'MFI normalized to mean MFI of all cells'
                        DataPrepMFI = 'Global Fold Change: Number of Cells / Mean Cells in Tissue Neighborhoods';
                    % Normalize to the maximum MFI in the sample
                    case 'MFI normalized to max MFI per neighborhood'
                        DataPrepMFI = 'Global Composition: Number of Cells / Max Cells in Tissue Neighborhoods';
                    case 'MFI normalized to max MFI of all cells'
                        DataPrepMFI = 'Global Composition: Number of Cells / Max Cells in Tissue Neighborhoods';
                    % normalize to neighborhood volume/area
                    case 'Density: sum(MFI) / Volume (Area if 2D) Per Neighborhood'
                        DataPrepMFI = 'Density: Number of Cells / Volume (Area if 2D) Per Neighborhood';
                    % Binary
                    case  'Binary: If MFI of cell > 1 = 1'
                        DataPrepMFI = 'Binary: If cell is in neighborhood = 1';
                    case 'Binary: If MFI in neighborhood > 1 = 1'
                        DataPrepMFI = 'Binary: If cell is in neighborhood = 1';
                    case 'Corrected Density: sum(MFI) / Volume (Area if 2D) of Neighborhood Inside Tissue'
                        DataPrepMFI = 'Corrected Density: Number of Cells / Volume (Area if 2D) of Neighborhood Inside Tissue';
                    case 'Composition+Density: sum(MFI) / Total cells in Neighborhood; Density of Neighborhood'
                        DataPrepMFI = 'Composition+Density: sum(MFI) / Total cells in Neighborhood; Density of Neighborhood';                          
                end
            end
            
            switch DataType
                case 'Raster Scanned Neighborhoods'
                    DType = 'MFIRSN';
                case 'MFIRSN'
                    DType = 'MFIRSN';
                case 'Cell Centered Neighborhoods'
                    DType = 'MFICCN';
                case 'MFICCN'
                    DType = 'MFICCN';
                case 'Individual Cells'
                    DType = 'AllCells';
            end
            
            % This should't be necesary but without it I get errors..
            app.DataN.Value = smplnms{1};

            % Prepare an array to keep track of which data goes with which
            % sample, make it [start, finish, sample name i;]
            INDDatPre = cell(numel(smplnms), 1);
            INDzrs = cell(numel(smplnms), 1);
            INDons = cell(numel(smplnms), 1);

            % Pull the phenotype names from the weights list
            Phenotypes = phenotypes;
            SortNames = Helper.gate_full2tag(app, Phenotypes, smplnms{1});
            
            Dat_Pre = [];
            Dat_MFI = [];
            
            % Load the first samples data
            Dat_PreALL = app.data.(smplnms{1}).(DType);
            if ~all(ismember(SortNames, Dat_PreALL.Properties.VariableNames))
                errordlg(["Not all chosen Phenotypes are available in " DataType]);
                return;
            end
            Dat_Pre = Dat_PreALL(:, SortNames);
            Dat_Pre = table2array(Dat_Pre);
            
            % Pull any MFI data
            Dat_MFI = table2array(Dat_PreALL(:, MFIList));

            if strcmp(DType, 'AllCells') && ~isempty(SortNames) 
                % Pull out the cell vectors, leaving only MFI
                logic = sum(Dat_Pre, 2);
                logic = logic~=0;
                % Pull out indeces to keep track of where things are
                INDzrs{1} = all(logic.*[Dat_Pre, Dat_MFI]==0,2)==1;
                INDons{1} = all(logic.*[Dat_Pre, Dat_MFI]==0,2)==0;
            else
                % Pull out indeces to keep track of where things are
                INDzrs{1} = all([Dat_Pre, Dat_MFI]==0,2)==1;
                INDons{1} = all([Dat_Pre, Dat_MFI]==0,2)==0;
            end

            if RmvZeros ==1
                % Delete any zero cell elements
                Dat_Pre(INDzrs{1},:)=[];
            end

            % Pull the length of the data
            INDDatPre{1} = {[1, size(Dat_Pre, 1)], smplnms{1}};

            if RmvZeros ==1
                Dat_MFI(INDzrs{1},:)=[];
            end
            %% Load the rest of the samples
            
            % If the user selected to normalize the data per sample
            switch NormOpt
                case 'Sample'
                    if RmvZeros ==1
                        % Prepare the data as selected by the user
                        Dat_Pre = Helper.Func_DataPrep(Dat_PreALL(INDons{1},:), Dat_Pre, DataPrep);
                        Dat_MFI = Helper.Func_DataPrep(Dat_PreALL(INDons{1},:), Dat_MFI, DataPrepMFI);
                    else
                        % Prepare the data as selected by the user
                        Dat_Pre = Helper.Func_DataPrep(Dat_PreALL, Dat_Pre, DataPrep);
                        Dat_MFI = Helper.Func_DataPrep(Dat_PreALL, Dat_MFI, DataPrepMFI);
                    end
                case 'Dataset'
                    if RmvZeros ==1
                        if any(strcmp(Dat_PreALL.Properties.VariableNames, 'Neigh_Volume'))
                            Dat_PreALL_Norm = Dat_PreALL(INDons{1},{'NCells', 'Neigh_Volume'});
                        elseif any(strcmp(Dat_PreALL.Properties.VariableNames, 'Neigh_Area'))
                            Dat_PreALL_Norm = Dat_PreALL(INDons{1},{'NCells', 'Neigh_Area'});
                        elseif strcmp(DType, 'AllCells')
                            Dat_PreALL_Norm = Dat_PreALL(:,1);
                            Dat_PreALL_Norm.NCells = ones(size(Dat_PreALL, 1),1);
                            Dat_PreALL_Norm.Neigh_Volume = ones(size(Dat_PreALL, 1),1);
                        end
                    else
                        if any(strcmp(Dat_PreALL.Properties.VariableNames, 'Neigh_Volume'))
                            Dat_PreALL_Norm = Dat_PreALL(:,{'NCells', 'Neigh_Volume'});
                        elseif any(strcmp(Dat_PreALL.Properties.VariableNames, 'Neigh_Area'))
                            Dat_PreALL_Norm = Dat_PreALL(:,{'NCells', 'Neigh_Area'});
                        elseif strcmp(DType, 'AllCells')
                            Dat_PreALL_Norm = Dat_PreALL(:,1);
                            Dat_PreALL_Norm.NCells = ones(size(Dat_PreALL, 1),1);
                            Dat_PreALL_Norm.Neigh_Volume = ones(size(Dat_PreALL, 1),1);
                        end
                    end
            end
            
            % Load in the other selected data and append it to Dat_Pre
            if numel(smplnms)>1
                for i=2:numel(smplnms)
                    % ToDo fix this so that if different samples have
                    % different Gate_XX definitions it still works
                    SortNames_i = Helper.gate_full2tag(app, Phenotypes, smplnms{i});
                    
                    %Pull the data in, in the same order as in sample 1
                    datTEMPALL = app.data.(smplnms{i}).(DType);
                    datTEMP = datTEMPALL(:, SortNames_i);
                    
                    datTEMP = table2array(datTEMP);
                    % Pull the MFI data
                    Temp_MFI = table2array(datTEMPALL(:, MFIList));

                    if strcmp(DType, 'AllCells') && ~isempty(SortNames)
                        % Pull out the cell vectors, leaving only MFI
                        logic = sum(datTEMP, 2);
                        logic = logic~=0;
                        % Pull out indeces to keep track of where zero elements are
                        INDzrs{i} = all(logic.*[datTEMP, Temp_MFI]==0,2)==1;
                        INDons{i} = all(logic.*[datTEMP, Temp_MFI]==0,2)==0;
                    else
                        % Pull out indeces to keep track of where zero elements are
                        INDzrs{i} = all([datTEMP, Temp_MFI]==0,2)==1;
                        INDons{i} = all([datTEMP, Temp_MFI]==0,2)==0;
                    end

                    if RmvZeros ==1
                        % Delete any zero cell elements
                        datTEMP(INDzrs{i},:)=[];
                    end
                    
                    if RmvZeros == 1
                        Temp_MFI(INDzrs{i},:)=[];
                    end
                    
                    switch NormOpt
                        case 'Sample'
                            if RmvZeros == 1
                                % Prepare the data as selected by the user
                                datTEMP = Helper.Func_DataPrep(datTEMPALL(INDons{i},:), datTEMP, DataPrep);
                                Temp_MFI = Helper.Func_DataPrep(datTEMPALL(INDons{i},:), Temp_MFI, DataPrepMFI);
                            else
                                % Prepare the data as selected by the user
                                datTEMP = Helper.Func_DataPrep(datTEMPALL, datTEMP, DataPrep);
                                Temp_MFI = Helper.Func_DataPrep(datTEMPALL, Temp_MFI, DataPrepMFI);
                            end
                        case 'Dataset'
                            if RmvZeros == 1
                                if any(strcmp(datTEMPALL.Properties.VariableNames, 'Neigh_Volume'))
                                    Dat_PreALL_Norm = [Dat_PreALL_Norm; datTEMPALL(INDons{i},{'NCells', 'Neigh_Volume'})];
                                elseif any(strcmp(datTEMPALL.Properties.VariableNames, 'Neigh_Area'))
                                    Dat_PreALL_Norm = [Dat_PreALL_Norm; datTEMPALL(INDons{i},{'NCells', 'Neigh_Area'})];
                                end
                            else
                                if any(strcmp(datTEMPALL.Properties.VariableNames, 'Neigh_Volume'))
                                    Dat_PreALL_Norm = [Dat_PreALL_Norm; datTEMPALL(:,{'NCells', 'Neigh_Volume'})];
                                elseif any(strcmp(datTEMPALL.Properties.VariableNames, 'Neigh_Area'))
                                    Dat_PreALL_Norm = [Dat_PreALL_Norm; datTEMPALL(:,{'NCells', 'Neigh_Area'})];
                                end
                            end
                    end
                    
                    % Keep track of where you put the new data
                    INDDatPre{i} = {[size(Dat_Pre, 1)+1, size(Dat_Pre, 1)+size(datTEMP, 1)], smplnms{i}};
                    %Append the new data on the end of Dat_Pre
                    Dat_Pre = [Dat_Pre; datTEMP];
                    Dat_MFI = [Dat_MFI; Temp_MFI];
                end
            end
            
            % If the user selected to normalize the whole dataset do that
            switch NormOpt
                case 'Dataset'
                    Dat_Pre = Helper.Func_DataPrep(Dat_PreALL_Norm, Dat_Pre, DataPrep);
                    Dat_MFI = Helper.Func_DataPrep(Dat_PreALL_Norm, Dat_MFI, DataPrepMFI);
            end
            
            % Multiply by all of the data by the assigned cell weights
            if ~strcmp(DataPrep, 'Composition+Density: Number of Cells / Total cells in Neighborhood; Density of Neighborhood')
                if ~isempty(pheno_weights)
                    Dat_Pre = Dat_Pre .* pheno_weights;
                end
            elseif ~isempty(pheno_weights)
                Dat_Pre(:, 1:(end-1)) = Dat_Pre(:, 1:(end-1)) .* pheno_weights;
            end
            if ~isempty(MFI_weights)
                Dat_MFI = Dat_MFI .* MFI_weights;
            end
            Dat_Pre = [Dat_Pre, Dat_MFI];
        end

        function Dat_Pre = Func_DataPrep(Dat_PreALL, Dat_Pre, DataPrep)
            % FUNC_DATAPREP Normalizes and/or prepares data in a user
            % defined way.
            %
            % Input:
            %   - Dat_PreALL - table - Indexed portion of whole table which Dat_Pre came from.
            %       It will be used to define variables which will be used
            %       in data preparation process.
            %   - Dat_Pre - double - Data which is supposed to be
            %       processed, and prepared.
            %   - DataPrep - string - What techinque to use when modifiying
            %       the data.
            %
            % Output:
            %   - Dat_Pre - double - Given Dat_Pre, but prepared in a way
            %       specified by user.
            %% Normalize the data    
            if strcmp(DataPrep, 'Cellularity: Number of Cells / Neighborhood')
                % Use the Raw cell numbers from CellularityM
            elseif strcmp(DataPrep, 'Density: Number of Cells / Volume (Area if 2D) Per Neighborhood')
                % Divide by the volume/area of the neighborhoods (Will probably look identical to Cell numbers map)
                % {'X', 'Y', 'Z', V, ['Effective_' V], 'NCells'}
                if any(strcmp(Dat_PreALL.Properties.VariableNames, 'Neigh_Volume'))
                    Dat_Pre = Dat_Pre./Dat_PreALL.Neigh_Volume;
                elseif any(strcmp(Dat_PreALL.Properties.VariableNames, 'Neigh_Area'))
                    Dat_Pre = Dat_Pre./Dat_PreALL.Neigh_Area;
                end
            elseif strcmp(DataPrep, 'Global Composition: Number of Cells / Max Cells in Tissue Neighborhoods')
                    % Divide by the neighborhood with the max number of cells
                    if size(Dat_Pre, 1)~=1
                        Dat_Pre = Dat_Pre./max(Dat_Pre);
                    end
            elseif strcmp(DataPrep, 'Global Fold Change: Number of Cells / Mean Cells in Tissue Neighborhoods')
                    % Divide by the neighborhood with the mean number of cells
                    if size(Dat_Pre, 1)~=1
                        Dat_Pre = Dat_Pre./mean(Dat_Pre);   
                    end
            elseif strcmp(DataPrep, 'Binary: If cell is in neighborhood')
                    % Choose neighborhoods which have more than 1 cell.
                    Dat_Pre = Dat_Pre > 1;
            elseif strcmp(DataPrep, 'Standardize: subtract Mean, divide by standard deviation')
                    % Put data into form of (mu=0, sigma=1)
                    if size(Dat_Pre, 1)~=1
                        Dat_Pre = Dat_Pre - mean(Dat_Pre);
                        stdeva = std(Dat_Pre);
                        stdeva(stdeva==0) = 1; % account for single cells
                        Dat_Pre = Dat_Pre./stdeva;
                    else
                        Dat_Pre = Dat_Pre.*0;
                    end
            elseif strcmp(DataPrep, 'Corrected Density: Number of Cells / Volume (Area if 2D) of Neighborhood Inside Tissue')
                % Divide by the volume/area of the neighborhood that is inside the tissue (This should correct for edge effects)
                if any(strcmp(Dat_PreALL.Properties.VariableNames, 'Neigh_Volume'))
                    Vol = Dat_PreALL.Effective_Neigh_Volume;
                    Vol(Vol==0) = 1;
                    Dat_Pre = Dat_Pre./Vol;
                elseif any(strcmp(Dat_PreALL.Properties.VariableNames, 'Neigh_Area'))
                    Vol = Dat_PreALL.Effective_Neigh_Area;
                    Vol(Vol==0) = 1;
                    Dat_Pre = Dat_Pre./Vol;
                end
            elseif strcmp(DataPrep, 'Composition: Number of Cells / Total cells in Neighborhood')
                % Divide by the total number of cells in the neighborhoods (This should also correct for edge effects)
                NCell = Dat_PreALL.NCells;
                NCell(NCell==0) = 1;
                Dat_Pre = Dat_Pre./NCell;
            elseif strcmp(DataPrep, 'Composition+Density: Number of Cells / Total cells in Neighborhood; Density of Neighborhood')
                % Divide by the total number of cells in the neighborhoods (This should also correct for edge effects)
                NCell = Dat_PreALL.NCells;
                NCell(NCell==0) = 1;
                Dat_Pre = Dat_Pre./NCell;
                % Scale the nuber of cells/neighborhood
                NCell = NCell./max(NCell);
                % Add that as a relative density column
                Dat_Pre = [Dat_Pre, 0.5.*NCell];
            end            
        end

        function mem = get_free_memory()
            if ispc
                mem = 500;
            elseif ismac
                % We do not have Mac.
                mem = 500;
            elseif isunix
                [~, cmdout] = system('cat /proc/meminfo');
                cmdout = split(cmdout, newline);
                cmdout = cmdout{startsWith(cmdout, 'MemAvailable')};
                cmdout = split(cmdout, ':');
                cmdout = cmdout{2};
                cmdout = strtrim(cmdout);
                magnitude = cmdout(end - 1);
                cmdout = split(cmdout, ' ');
                cmdout = cmdout{1};

                mem = str2double(cmdout);
                switch magnitude
                    case 'k'
                        mem = mem * 10e3;
                    case 'm'
                        mem = mem * 10e6;
                    case 'g'
                        mem = mem * 10e9;
                    case 'p'
                        mem = mem * 10e12;
                    otherwise
                end
            else
                mem = 500;
            end
            
        end

        function func_RunExtensions(app, btn, ExtenMenu)
            %% Populate the drop down:
            % Pull the names of the functions in the Custom Functions folder
            MyFolderInfo = dir(fullfile([app.CurrentPath filesep 'Extensions']));
            % remove weird stuff
            MyFolderInfo = MyFolderInfo(~strcmp('.', {MyFolderInfo.name}));
            MyFolderInfo = MyFolderInfo(~strcmp('..', {MyFolderInfo.name}));
% % %             MyFolderInfo = dir(fullfile([app.CurrentPath filesep 'Extensions'], '*.m'));
            MyFolderInfo = MyFolderInfo(~cellfun('isempty', {MyFolderInfo.date}));
            extennames = cell(size(MyFolderInfo,1),1);
            % if there are no extensions
            if isempty(MyFolderInfo)
                MyFolderInfo(1).name = 'No Extensions Found';
            end
            
            for i=1:size(MyFolderInfo,1)
                extennames{i}=MyFolderInfo(i).name;
                % Check to see if there is anything in the drop down menu
                if isempty(ExtenMenu.Children)
                    % Create an entry in the dropdown
                    ext = uimenu(ExtenMenu);
                    ext.MenuSelectedFcn = @(ext,event) Helper.func_RunExtensions(app, ext, ExtenMenu);
                    ext.Text = MyFolderInfo(i).name;  
                end
            end
            currentNames = {ExtenMenu.Children.Text};
            % See if there are any new extensions to add
            IND = find(~ismember(extennames,currentNames));
            if ~isempty(IND)
                for i = 1:numel(IND)
                    % Create an entry in the dropdown
                    ext = uimenu(ExtenMenu);
                    ext.MenuSelectedFcn = @(ext,event) Helper.func_RunExtensions(app, ext, ExtenMenu);
                    ext.Text = extennames{IND(i)};            
                end
            end
            % Run the extension
            if ~isempty(btn)
                if strcmp(btn.Text, 'No Extensions Found')
                    app.CurrentPath = uigetdir(app.CurrentPath, 'Select Extensions Folder');
                    if ~isempty(app.CurrentPath)
                        Helper.func_RunExtensions(app, [], ExtenMenu);
                    end
                    return
                end
                temp=split(btn.Text,'.');
                func=string(temp(1))+'(app)';
                addpath(fullfile([app.CurrentPath filesep 'Extensions']));  
                eval(func);
            end
        end
        
        function func_exit(app)
            % FUNC_EXIT Exit Program function
% %             % Close all of your figures TODO Make this an option dialog when
% closing
% %             FigN = get(groot, 'Children');
% %             if ~isempty(FigN)
% %                 FigN = [FigN.Number];
% %                 for i=1:numel(FigN)
% %                     close(figure(FigN(i)))
% %                 end
% %             end
            % Close the SAS3D and PhList objects
            delete(app.PhList)
            delete(app.GUI.Main)
        end
  
        function decVect = bi2de(MAT)
            % convert binary matrix to unique decimals
            Ncols = size(MAT, 2);
            for col_i = 1:Ncols
                col_i_dec = (2^(col_i-1));
                MAT(:, col_i) = MAT(:, col_i).*col_i_dec;
            end
            decVect = sum(MAT,2);
        end
        
        function MAT = de2bi(decVect)
            MAXdec = max(decVect);    
            ncol = 0;
            while MAXdec~=0
                ncol = ncol+1;
                R = rem(MAXdec,2);
                MAXdec = (MAXdec-R)/2;
                if ncol>500
                    break
                end
            end
            MAT = zeros(numel(decVect), ncol);
            for col_i = 1:ncol
                R = rem(decVect,2);
                decVect = (decVect-R)/2;
                MAT(:,col_i) = R;
            end
        end
        
        %% Helpers related to clusterd, loading, and normalizing data
        
        function [smplnms, type] = check_model(app, smplnms, DataType, MFIList)
            %% This function makes sure the selected model or regions are in all samples
            % It takes the sample names and leaves out any samples that
            % don't include the selected model
            is_in = ones(numel(smplnms), 1);
            for smpl_idx=1:numel(smplnms)
                if strcmp(DataType, 'Manually Defined Regions (RSN)')
                    % for now manual regions only works ith raster scanned
                    % neighborhoods
                    type = 'Raster Scanned Neighborhoods';
% % %                     IND_reg = startsWith(MFIList, Constants.neigh_tag);
% % %                     ManRegs = MFIList(IND_reg);
% % %     % %                 MFIList = MFIList(~IND_reg);
% % %                     if ~ismember(...
% % %                             ManRegs, ...
% % %                             app.data.(smplnms{smpl_idx}).MFIRSN.Properties.VariableNames...
% % %                     )
% % %                         is_in(smpl_idx) = 0;
% % %                     end
                elseif strcmp(DataType, 'Manually Defined Regions (CCN)')
                    % for now manual regions only works ith raster scanned
                    % neighborhoods
                    type = 'Cell Centered Neighborhoods';
% % %                     IND_reg = startsWith(MFIList, Constants.neigh_tag);
% % %                     ManRegs = MFIList(IND_reg);
% % %     % %                 MFIList = MFIList(~IND_reg);
% % %                     if ~ismember(...
% % %                             ManRegs, ...
% % %                             app.data.(smplnms{smpl_idx}).MFICCN.Properties.VariableNames...
% % %                     )
% % %                         is_in(smpl_idx) = 0;
% % %                     end 
                else % if the user selected a model
                    if strcmp(app.net.(DataType).userdata.DataType, 'Raster Scanned Neighborhoods')
                        type = 'Raster Scanned Neighborhoods';
                        if ~ismember(...
                                [Constants.other_tag, DataType], ...
                                app.data.(smplnms{smpl_idx}).MFIRSN.Properties.VariableNames...
                        )
                            is_in(smpl_idx) = 0;
                        end
                    elseif strcmp(app.net.(DataType).userdata.DataType, 'Cell Centered Neighborhoods')
                        type = 'Cell Centered Neighborhoods';
                        if ~ismember(...
                                [Constants.other_tag, DataType], ...
                                app.data.(smplnms{smpl_idx}).MFICCN.Properties.VariableNames...
                        )
                            is_in(smpl_idx) = 0;
                        end
                    elseif strcmp(app.net.(DataType).userdata.DataType, 'Individual Cells')
                        type = 'Individual Cells';
                    end

                end
            end

            is_in = logical(is_in);
            if ~all(is_in)
                if any(is_in)
                    warndlg( ...
                        "Samples: " + join(smplnms(~is_in), ", ") + ...
                        " do not contain the chosen model." + newline, ...
                        "However samples: " + join(smplnms(is_in), ",") + ...
                        " contain this model, and process will continue with only these samples."...
                    );
                    smplnms = smplnms(is_in);
                else
                    errordlg( ...
                        "All of the samples chosen do not contain the chosen model." + newline + ...
                        "The program will abort." + newline + ...
                        "In order to choose these samples, reuse the model on these datasets."...
                    );
                    return;
                end
            end
        end % end function check model
        
        function func_annotate_regions(app, UIfig, mdl)
            %% Build a user interface
            % annotate clusters function builds a user interface for the user to
            % define the names of clustered cells
            %
            % Input:
            %   - app - Instance of SAS3D
            if ~Helper.any_sample(app) || ~Helper.any_net(app)
                return;
            end
            
            if nargin ~= 3
                % Build main UI window
                UIfig = uifigure('Name', 'Annotate Regions', 'Scrollable', 'on');
                Helper.func_SetCLR(app, UIfig, 'UIfigure');
                UIfig.Position = [10 10 350 700];
            end
            
            % Select clustering Model
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Clustering Model:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = [5 700-20 400 15];
            ModelNames = uidropdown(UIfig);
            ModelNames.Items = fieldnames(app.net);
            if nargin == 3
                ModelNames.Value = mdl;
            else
                ModelNames.Value = ModelNames.Items{1};
            end
            ModelNames.Position = [5 700-30-20 150 30];
            ModelNames.BackgroundColor = app.GUIOPTS.bgclr;
        %%%%%%%%%%%%%%%%%%%

            if ismember('CellNames', fieldnames(app.net.(ModelNames.Value)))
                CellNamesTMP = app.net.(ModelNames.Value).CellNames;
                if ismember('N', fieldnames(CellNamesTMP))
                    CellNamesTMP.N = app.net.(ModelNames.Value).NR;
                end
                CellNamesTMP.names = struct;
                CellNamesTMP.clusters = struct;
                CellNamesTMP.btnrmv = struct;
                CellNamesTMP.btnclr = struct;
            else
                % Make a local copy of these tied to this figure
                CellNamesTMP = struct;
                CellNamesTMP.N = app.net.(ModelNames.Value).NR;
                CellNamesTMP.names = struct;
                CellNamesTMP.clusters = struct;
                CellNamesTMP.btnrmv = struct;
                CellNamesTMP.btnclr = struct;
            end

        %%%%%%%%%%%%%%%%%%%

            % Assign Region Names
            lbl = uilabel(UIfig); lbl.Text = ...
                sprintf('Assign Names to Regions:               Regions: 1,3, etc.');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = [5 700-40-3*20 400 15];

            for btni = 1:CellNamesTMP.N
                % assign first name                
                nme = ['name_' num2str(btni)];
                CellNamesTMP.names.(nme) = uidropdown(UIfig,'Editable','on');
                CellNamesTMP.names.(nme).Position = [5 700-40-4*20-20*(btni-1) 175 20];
                if ismember('RegNames', fieldnames(app.net.(ModelNames.Value)))
                    CellNamesTMP.names.(nme).Items = app.net.(ModelNames.Value).RegNames(btni);
                    CellNamesTMP.names.(nme).Value = app.net.(ModelNames.Value).RegNames{btni};
                else
                    CellNamesTMP.names.(nme).Items = {nme};
                    CellNamesTMP.names.(nme).Value = nme;
                end

                % assign clusters to this name
                clust = ['clust_' num2str(btni)];
                CellNamesTMP.clusters.(clust) = uieditfield(UIfig,'text');
                CellNamesTMP.clusters.(clust).Position = [5+175 700-40-4*20-20*(btni-1) 80 20];
                CellNamesTMP.clusters.(clust).Value = num2str(btni);

                % add a set region Color button
                clr = ['btnclr_' num2str(btni)];
                CellNamesTMP.btnclr.(clr) = uibutton(UIfig,'push');
                CellNamesTMP.btnclr.(clr).Position = [5+175+80, 700-40-4*20-20*(btni-1), 40, 20];
                CellNamesTMP.btnclr.(clr).Text = 'Color';
                Helper.func_SetCLR(app, CellNamesTMP.btnclr.(clr), 'button')
                CellNamesTMP.btnclr.(clr).BackgroundColor = app.net.(ModelNames.Value).cmap(btni+1, :);
                CellNamesTMP.btnclr.(clr).ButtonPushedFcn = ...
                    @(~,~) func_btncolor(app, CellNamesTMP.btnclr, btni);
            end

            % Reset Colors
            rstclr = uibutton(UIfig,'push');
            rstclr.Position = [5+150 700-30-20 150 30];
            rstclr.Text = 'Select Colorscheme';
            Helper.func_SetCLR(app, rstclr, 'button')
            
            % Save Annotations
            svbtn = uibutton(UIfig,'push');
            svbtn.Position = [5 700-40-2*20 150 30];
            svbtn.Text = 'Save Annotations';
            Helper.func_SetCLR(app, svbtn, 'button')

            % Export Neighborhoods
            exbtn = uibutton(UIfig,'push');
            exbtn.Position = [5+150 700-40-2*20 150 30];
            exbtn.Text = 'Save Neighborhoods';
            Helper.func_SetCLR(app, exbtn, 'button')

            % add a new name
            btn = uibutton(UIfig,'push');
            btn.Position = [5, 700-40-5*20-20*(btni-1), 20, 20];
            btn.Text = '+';
            Helper.func_SetCLR(app, btn, 'button')
            
            % remove name button
            CellNamesTMP.btnrmv.btnrmv_1 = uibutton(UIfig,'push');
            CellNamesTMP.btnrmv.btnrmv_1.Position = [5+175+80+20+20, 700-40-4*20-20*(btni-1), 20, 20];
            CellNamesTMP.btnrmv.btnrmv_1.Text = '-';
            Helper.func_SetCLR(app, CellNamesTMP.btnrmv.btnrmv_1, 'button')

            % Set the button pushed functions
            rstclr.ButtonPushedFcn = @(~,~) func_reset_colors(app, CellNamesTMP);
            svbtn.ButtonPushedFcn = @(~,~) func_save_annotations(app,ModelNames, CellNamesTMP);
            exbtn.ButtonPushedFcn = @(~,~) func_export_annotations(app,ModelNames, CellNamesTMP);

            btn.ButtonPushedFcn = @(~,~) func_new_name(app, ModelNames, UIfig, btn, svbtn, rstclr, CellNamesTMP);
            
            CellNamesTMP.btnrmv.btnrmv_1.ButtonPushedFcn = ...
                @(~,~) func_btnrmv(app, ModelNames, UIfig, btn, svbtn, rstclr, num2str(btni), CellNamesTMP, true);

            ModelNames.ValueChangedFcn = @(btn, event) ChangeDataType(app, ModelNames, UIfig);
            
            function func_reset_colors(app, CellNamesTMP)
                clrlist = {'jet', ...
                        'bone', ...
                        'autumn', ...
                        'prism', ...
                        'hsv', ...
                        'hot', ...
                        };
                [indx,tf] = listdlg('ListString',clrlist, 'SelectionMode', 'single');
                if tf==0
                    return
                end
                exprsn = [clrlist{indx} '(CellNamesTMP.N);'];
                clrmp = eval(exprsn);
                if indx~=1
                    clrmp = flipud(clrmp);
                end
                for clri = 1:CellNamesTMP.N
                    % add a set region Color button
                    clrnmi = ['btnclr_' num2str(clri)];
                    CellNamesTMP.btnclr.(clrnmi).BackgroundColor = clrmp(clri, :);
                    CellNamesTMP.btnclr.(clrnmi).ButtonPushedFcn = ...
                        @(~,~) func_btncolor(app, CellNamesTMP.btnclr, clri);
                end
            end

            function func_new_name(app, ModelNames, UIfig, btn, svbtn, rstclr, CellNamesTMP)

                CellNamesTMP.N = CellNamesTMP.N+1;

                nms = strrep(fieldnames(CellNamesTMP.names), 'name_', '');
                S = sprintf('%s ', nms{:});
                numi = num2str(max(sscanf(S, '%f')+1));

                name_i = ['name_' numi];
                CellNamesTMP.names.(name_i) = uidropdown(UIfig,'Editable','on');
                CellNamesTMP.names.(name_i).Position = [5 btn.Position(2) 175 20];
                CellNamesTMP.names.(name_i).Items ...
                    = fieldnames(CellNamesTMP.names);
                CellNamesTMP.names.(name_i).Value = name_i;

                % assign clusters to this name
                clust_i = ['clust_' numi];
                CellNamesTMP.clusters.(clust_i) = uieditfield(UIfig,'text');
                CellNamesTMP.clusters.(clust_i).Position = [5+175 btn.Position(2) 80 20];
                CellNamesTMP.clusters.(clust_i).Value = '1';

                btnclr_i = ['btnclr_' numi];
                CellNamesTMP.btnclr.(btnclr_i) = uibutton(UIfig,'push');
                CellNamesTMP.btnclr.(btnclr_i).Position = [5+175+80, btn.Position(2), 40, 20];
                CellNamesTMP.btnclr.(btnclr_i).Text = 'Color';
                CellNamesTMP.btnclr.(btnclr_i).BackgroundColor = [0.8, 0.8, 0.8];

                CellNamesTMP.btnrmv.btnrmv_1.Position(2) = ...
                    CellNamesTMP.btnrmv.btnrmv_1.Position(2)-20;

                CellNamesTMP.btnrmv.btnrmv_1.ButtonPushedFcn = ...
                    @(~,~) func_btnrmv(app, ModelNames, UIfig, btn, svbtn, rstclr, numi, CellNamesTMP, true);

                CellNamesTMP.btnclr.(btnclr_i).ButtonPushedFcn = ...
                    @(~,~) func_btncolor(app,CellNamesTMP.btnclr, str2double(numi));


                btn.Position(2) = btn.Position(2)-20;

                btn.ButtonPushedFcn = @(~,~) func_new_name(app,ModelNames, UIfig, btn, svbtn, rstclr, CellNamesTMP);
                
                svbtn.ButtonPushedFcn = @(~,~) func_save_annotations(app,ModelNames, CellNamesTMP);
                rstclr.ButtonPushedFcn = @(~,~) func_reset_colors(app, CellNamesTMP);

                ModelNames.ValueChangedFcn = @(btn, event) ...
                    ChangeDataType(app, ModelNames, UIfig);
            end

            function ChangeDataType(app, ModelNames, UIfig)
% % %                 if ismember('CellNames', fieldnames(app.net.(ModelNames.Value)))
                    mnm = ModelNames.Value;
                    h=findall(UIfig);
                    delete(h(2:end))
                    Helper.func_annotate_regions(app, UIfig, mnm)
% % %                 else
                    % if there are no current cell names do nothing
% % %                 end

            end

            function func_btnrmv(app, ModelNames, UIfig, btn, svbtn, rstclr, numi, CellNamesTMP, chngdata)
                %%
                %get current name
                name_i = ['name_' numi];
                clust_i = ['clust_' numi];
                btnrmv_i = ['btnrmv_' numi];
                btnclr_i = ['btnclr_' numi];

                nms = fieldnames(CellNamesTMP.names);
                cls = fieldnames(CellNamesTMP.clusters);
                nmind = find(ismember(nms,name_i));

                CellNamesTMP.N = CellNamesTMP.N-1;
                % turn off buttons to be removed
                CellNamesTMP.names.(name_i).Visible = 'off';
                CellNamesTMP.clusters.(clust_i).Visible = 'off';
                CellNamesTMP.btnclr.(btnclr_i).Visible = 'off';

        % % %         btnrmv_i.Visible = 'off';

                CellNamesTMP.btnrmv.btnrmv_1.Position(2) = ...
                    CellNamesTMP.btnrmv.btnrmv_1.Position(2)+20;
% % % 
% % %                 CellNamesTMP.btnrmv.btnrmv_1.ButtonPushedFcn = ...
% % %                     @(~,~) func_btnrmv(app,ModelNames, UIfig, btn, svbtn, num2str(str2double(numi)-1), CellNamesTMP);


                %% Remove the selected name and clusters

        % % %         CellNamesTMP.names.(name_i) = [];
        % % %         CellNamesTMP.clusters.(clust_i) = [];
                CellNamesTMP.names = ...
                    rmfield(CellNamesTMP.names, name_i);
                CellNamesTMP.clusters = ...
                    rmfield(CellNamesTMP.clusters, clust_i);
                CellNamesTMP.btnclr = ...
                    rmfield(CellNamesTMP.btnclr, btnclr_i);

                % move the + button up one
                btn.Position(2) = btn.Position(2)+20;

                %% move other buttons up
                
                if nmind ~= numel(nms)
                    for nm_i = 1:(numel(nms)-nmind)
                        if contains((['name_' num2str(nm_i + nmind)]), fieldnames(CellNamesTMP.names))
                            if ~isempty(CellNamesTMP.names.(['name_' num2str(nm_i + nmind)]))                        
                                CellNamesTMP.names.(['name_' num2str(nm_i + nmind)]).Position(2) = ...
                                 CellNamesTMP.names.(['name_' num2str(nm_i + nmind)]).Position(2)+20; 
                                CellNamesTMP.clusters.(['clust_' num2str(nm_i + nmind)]).Position(2) = ...
                                 CellNamesTMP.clusters.(['clust_' num2str(nm_i + nmind)]).Position(2)+20;
% %                                 CellNamesTMP.btnrmv.(['btnrmv_' num2str(nm_i + nmind)]).Position(2) = ...
% %                                  CellNamesTMP.btnrmv.(['btnrmv_' num2str(nm_i + nmind)]).Position(2)+20;
                                CellNamesTMP.btnrmv.btnrmv_1.Position(2) = ...
                                 CellNamesTMP.btnrmv.btnrmv_1.Position(2)+20;
                                CellNamesTMP.btnclr.(['btnclr_' num2str(nm_i + nmind)]).Position(2) = ...
                                 CellNamesTMP.btnclr.(['btnclr_' num2str(nm_i + nmind)]).Position(2)+20;
                            end
                        end
                    end
                end
                   
                CellNamesTMP.btnrmv.btnrmv_1.ButtonPushedFcn = ...
                    @(~,~) func_btnrmv(app,ModelNames, UIfig, btn, svbtn, rstclr, num2str(str2double(numi)-1), CellNamesTMP, true);
                
                btn.ButtonPushedFcn = @(~,~) func_new_name(app, ModelNames, UIfig, btn, svbtn, rstclr, CellNamesTMP);
                rstclr.ButtonPushedFcn = @(~,~) func_reset_colors(app, CellNamesTMP);
                svbtn.ButtonPushedFcn = @(~,~) func_save_annotations(app,ModelNames, CellNamesTMP);
                    ModelNames.ValueChangedFcn = @(btn, event) ...
                        ChangeDataType(app, ModelNames, UIfig);
            % %         Nrmv = str2num(strrep(name_i, 'name_', ''))
            % %         fieldnames(names) 
            end

            function func_btncolor(app, btnclr , ind)
                %%
                btnclr_i = ['btnclr_' num2str(ind)];
                selected = btnclr.(btnclr_i).BackgroundColor;
                c = uisetcolor(selected);
                btnclr.(btnclr_i).BackgroundColor = c;

            end

            function func_save_annotations(app, ModelNames, CellNamesTMP)
                
                RegNames = cell(CellNamesTMP.N, 1);
                Regs = cell(CellNamesTMP.N, 1);
                rgb = ones(CellNamesTMP.N+1, 3);
                newORDER  = zeros(1, app.net.(ModelNames.Value).NR+1);

                for Reg_i = 1:CellNamesTMP.N
                    name_i = ['name_' num2str(Reg_i)];
                    clust_i = ['clust_' num2str(Reg_i)];
                    btnclr_i = ['btnclr_' num2str(Reg_i)];


                    RegNames{Reg_i} = CellNamesTMP.names.(name_i).Value;
                    Regs{Reg_i} = str2double(strsplit(CellNamesTMP.clusters.(clust_i).Value, ','));            
                    for reg_i = Regs{Reg_i}
                        if (reg_i+1) <=numel(newORDER)
                            newORDER(reg_i+1) = Reg_i;
                        else
                            warning(['Region ' num2str(reg_i)  ' not found'])
                        end
                    end
                    rgb(Reg_i+1, :) = CellNamesTMP.btnclr.(btnclr_i).BackgroundColor;

                end
                save_annotation(app, ModelNames.Value, RegNames, newORDER, rgb, CellNamesTMP)

            end

            function func_export_annotations(app,ModelNames, CellNamesTMP)

                for cell_i = 1:CellNamesTMP.N
                    name_i = ['name_' num2str(cell_i)];
                    clust_i = ['clust_' num2str(cell_i)];

                    CellNames = CellNamesTMP.names.(name_i).Value;
                    Clusters = str2double(strsplit(CellNamesTMP.clusters.(clust_i).Value, ','));
                    gate_annotation(app, ModelNames, CellNames, Clusters)
                end

            end

            function save_annotation(app, ModelName, RegNames, newORDER, rgb, CellNamesTMP)
                % ask the user to name the model
                prompt = {'Enter name of annotated model:'};
                title = 'Name model';
                dims = [1 35];
                definput = {[ModelName '_Annotated' ]};
                NewName = inputdlg(prompt,title,dims,definput);
                if isempty(NewName)
                    return;
                end
                % Make sure the user actually wants to overwrite an existing model
                if any(ismember(fieldnames(app.net), NewName))
                    answer = questdlg(['Are you certain you want to overwrite' NewName], ...
                        'Dessert Menu', ...
                        'Yes', 'No', 'No');
                    if strcmp(answer, 'No') || isempty(answer)
                        % User decided not to override data, abort 
                        return;
                    end
                    % If the user chose to override an existing model do
                    % that
                    app.net.(NewName{1}) = app.net.(ModelName);
                else
                    % If the model does not exist add it now
                    app.net.(NewName{1}) = app.net.(ModelName);
                end

                app.net.(NewName{1}).CellNames = struct;
                app.net.(NewName{1}).CellNames.N = CellNamesTMP.N;
                app.net.(NewName{1}).CellNames.ClusterOrder = newORDER;
                app.net.(NewName{1}).CellNames.Parent = ModelName;

                % Edit the region colormap
                app.net.(NewName{1}).cmap = rgb;

                % Add the new region names
                app.net.(NewName{1}).RegNames = RegNames;

                % Pull all smaple names
                smpls = fieldnames(app.data);

                % Define the original numbering scheme
                oldORDER = 0:app.net.(ModelName).NR;
                %% Check to see if the region numbers have been chaged
                for sample_i = 1:numel(smpls)
                    % Now that we are sure the user really wants to do this,
                    % Check to see if the region numbers have been chaged
                    
                    % rearrange the numbers in all of the data
                    % initialize some stuff
                    RowOrder = zeros(1,2);
                    ORDER = newORDER;
                    if ismember([Constants.other_tag ModelName], ...
                            fieldnames(app.data.(smpls{sample_i}).AllCells))
                        ROWAC = app.data.(smpls{sample_i}).AllCells.([Constants.other_tag ModelName]);
                        IND_AC = zeros(numel(ROWAC), numel(ORDER));
                    end
                    if ismember('MFIRSN', fieldnames(app.data.(smpls{sample_i})))
                        if ismember([Constants.other_tag ModelName], ...
                                fieldnames(app.data.(smpls{sample_i}).MFIRSN))
                        ROWRSN = app.data.(smpls{sample_i}).MFIRSN.([Constants.other_tag ModelName]);
                        IND_RSN = zeros(numel(ROWRSN), numel(ORDER));
                        end
                    end
                    if ismember('MFICCN', fieldnames(app.data.(smpls{sample_i})))
                        if ismember([Constants.other_tag ModelName], ...
                                fieldnames(app.data.(smpls{sample_i}).MFICCN))
                        ROWCCN = app.data.(smpls{sample_i}).MFICCN.([Constants.other_tag ModelName]);
                        IND_CCN = zeros(numel(ROWCCN), numel(ORDER));
                        end
                    end

                    %% Loop through the region numbers and re-define them
                    if any(oldORDER~=newORDER)
                        for i=1:numel(ORDER)
                            % If the model is in this AllCells in sample
                            if ismember([Constants.other_tag ModelName], ...
                                    fieldnames(app.data.(smpls{sample_i}).AllCells))
                                if ORDER(i) ~= i-1
                                    RowOrder = [RowOrder; i-1, ORDER(i)];
                                    % Pull the indeces of the elements
                                    IND_AC(:, i) = ROWAC==(i-1);
                                end
                            end
                            % If the model is in MFIRSN
                            if ismember('MFIRSN', fieldnames(app.data.(smpls{sample_i})))
                                if ismember([Constants.other_tag ModelName], ...
                                        fieldnames(app.data.(smpls{sample_i}).MFIRSN))
                                    if ORDER(i) ~= i-1
                                        IND_RSN(:, i) = ROWRSN==(i-1);
                                    end
                                end
                            end
                            % If the model is in MFICCN
                            if ismember('MFICCN', fieldnames(app.data.(smpls{sample_i})))
                                if ismember([Constants.other_tag ModelName], ...
                                        fieldnames(app.data.(smpls{sample_i}).MFICCN))
                                    if ORDER(i) ~= i-1
                                        IND_CCN(:, i) = ROWCCN==(i-1);
                                    end
                                end
                            end
                        end % end loop through region numbers
                    end % end if there were redefined numbers
                    %% Save the new vectors
                    for i=1:numel(ORDER)
                        % If the model is in this AllCells in sample
                        if ismember([Constants.other_tag ModelName], ...
                                fieldnames(app.data.(smpls{sample_i}).AllCells))
                            if ORDER(i) ~= i-1
                                % Change the elements the indeces of the elements
                                ROWAC(IND_AC(:, i)==1) = ORDER(i);
                            end
                            app.data.(smpls{sample_i}).AllCells.([Constants.other_tag NewName{1}]) = ROWAC;
                        end
                        % If the model is in MFIRSN
                        if ismember('MFIRSN', fieldnames(app.data.(smpls{sample_i})))
                            if ismember([Constants.other_tag ModelName], ...
                                    fieldnames(app.data.(smpls{sample_i}).MFIRSN))
                                if ORDER(i) ~= i-1
                                    ROWRSN(IND_RSN(:, i)==1) = ORDER(i);
                                end
                                app.data.(smpls{sample_i}).MFIRSN.([Constants.other_tag NewName{1}]) = ROWRSN;
                            end
                        end
                        % If the model is in MFICCN
                        if ismember('MFICCN', fieldnames(app.data.(smpls{sample_i})))
                            if ismember([Constants.other_tag ModelName], ...
                                    fieldnames(app.data.(smpls{sample_i}).MFICCN))
                                if ORDER(i) ~= i-1
                                    ROWCCN(IND_CCN(:, i)==1) = ORDER(i);
                                end
                                app.data.(smpls{sample_i}).MFICCN.([Constants.other_tag NewName{1}]) = ROWCCN;
                            end
                        end
                    end % end second loop through region numbers
                    %reset the number of regions
                    app.net.(NewName{1}).NR = app.net.(NewName{1}).CellNames.N;
                end % end loop through samples
            end

            function gate_annotation(app, ModelNames, name, Clusters)
                % save_annotation It is a function that saves the annotations as new cell populations.
                %
                % Input:
                %   - app - Instance of SAS3D
                %   - num - Identifier of figure from which to pull gate and
                %       data which will be used in saving process.
                %
                % Modifies:
                %   - app - Specifically fields in app.data.(smpl), where smpl
                %       are all samples on which this gate was defined.\

                valid_gate_name = Helper.valid_gate(name);
                if startsWith(valid_gate_name, 'Other')
                    valid_gate_name = 'Gate_Other';
                end
                if iscell(valid_gate_name)
                    valid_gate_name = valid_gate_name{1};
                end

                wb = waitbar(0, 'Please wait. Adding Annotations.');

                % import sample names that contain the chosen model
                smpls = fieldnames(app.data);
                for smpl_i = 1:numel(smpls)
                    if ~contains([Constants.other_tag ModelNames.Value], fieldnames(app.data.(smpls{smpl_i}).AllCells))
                        smpls{smpl_i} = [];
                    end
                end

                switch app.net.(ModelNames.Value).userdata.DataType
                    case 'Raster Scanned Neighborhoods'
                        phenos = {'Density/MFI RSN'};
                    case 'Cell Centered Neighborhoods'
                        phenos = {'Density/MFI CCN'};
                end

                type = 'Histogram';

                % Precompute stuff for non-global neighborhoods (ones, starting with a tag).
                in_rsn = true(numel(phenos), 1);
                if all(startsWith(phenos, Constants.neigh_tag))
                    % Determine if it's MFI RSN or CCN. Consistently for all Samples.
                    for phn_idx = 1:numel(phenos)
                        phn = phenos{phn_idx};
                        for smpl_idx = 1:numel(smpls)
                            smpl = smpls{smpl_idx};
                            if ~isfield(app.data.(smpl), 'MFIRSN') || ...
                                    ~ismember({phn}, app.data.(smpl).MFIRSN.Properties.VariableNames)
                                in_rsn(phn_idx) = false;
                            end
                        end
                    end
                end

                waitbar(.05, wb, 'Saving to the samples');
                for smpl_idx=1:numel(smpls)

                    smpl = smpls{smpl_idx};
                    if numel(phenos) == 1 && strcmp(phenos, 'Density/MFI RSN')
                        smpl_dat = app.data.(smpl).MFIRSN;
                    elseif numel(phenos) == 1 && strcmp(phenos, 'Density/MFI CCN')
                        smpl_dat = app.data.(smpl).MFICCN;
                    elseif numel(phenos) == 1 && strcmp(phenos, 'All Cells')
                        smpl_dat = app.data.(smpl).AllCells;
                    elseif all(startsWith(phenos, Constants.neigh_tag))
                        smpl_dat = {};
                        for phn_idx = 1:numel(phenos)
                            phn = phenos{phn_idx};
                            % Load all the neighborhoods into dat.
                            if in_rsn(phn_idx)
                                smpl_dat{phn_idx} = app.data.(smpl).MFIRSN(app.data.(smpl).MFIRSN.(phn), :);
                            else
                                smpl_dat{phn_idx} = app.data.(smpl).MFICCN(app.data.(smpl).MFICCN.(phn), :);
                            end
                        end
                    else
                        smpl_dat = {};
                        % Pull the data that is to be plotted
                        for phn_idx = 1:numel(phenos)
                            phn = phenos{phn_idx};
                            if ismember({phn}, {'Density/MFI RSN', 'Density/MFI CCN', 'All Cells'}) || ...
                                startsWith({phn}, Constants.neigh_tag)
                                % Just a sanity check it shouldn't be possible to come in here.
                                warndlg("Something went wrong. Neighborhoods and phenotypes shouldn't be mixed. You might need to restart CytoMap.");
                                continue;
                            end
                            tag = Helper.gate_full2tag(app, phn, smpl);
                            if iscell(tag)
                                tag = tag{1};
                            end
                            smpl_dat{phn_idx} = app.data.(smpl).AllCells(app.data.(smpl).AllCells.(tag)==1, :);
                        end
                    end

                    switch type
                        case 'Histogram'
                            points = Clusters;
                            logical_smpl = ismember(smpl_dat.([Constants.other_tag ModelNames.Value]),Clusters);
                    end

                    if ~iscell(smpl_dat)
                        if numel(phenos) == 1 && strcmp(phenos, 'Density/MFI RSN')
                            % Add the gated data to the sample
                            valid_name = Helper.valid_neigh(name);
                            if iscell(valid_name)
                                valid_name = valid_name{1};
                            end
                            if startsWith(valid_name, 'Other')
                                valid_name = 'Nei_Other';
                            end
                            app.data.(smpl).MFIRSN.(valid_name) = logical_smpl;

                            continue;
                        elseif numel(phenos) == 1 && strcmp(phenos, 'Density/MFI CCN')
                            % Add the gated data to the sample
                            valid_name = Helper.valid_neigh(name);
                            if iscell(valid_name)
                                valid_name = valid_name{1};
                            end
                            app.data.(smpl).MFICCN.(valid_name) = logical_smpl;
                            continue;
                        elseif numel(phenos) == 1 && strcmp(phenos, 'All Cells')
                            short_name = strcat('All/', name);
                            if iscell(short_name)
                                short_name = short_name{1};
                            end
                            tag = Helper.get_tag(app, short_name, smpl);
                            if iscell(tag)
                                tag = tag{1};
                            end
                            path = strcat('Gate_All/', Helper.valid_gate(name));
                            if iscell(path)
                                path = path{1};
                            end

                            new_tree = tree(name, 'gate_points', points, ...
                                'gate_axes', {[Constants.other_tag ModelNames.Value]}, ...
                                'gate_type', lower(type), ...
                                'tag', tag);

                            app.data.(smpl).AllCells.(tag) = logical_smpl;
                            app.data.(smpl).GateTags.(tag) = {path; short_name};
                            app.data.(smpl).tree = app.data.(smpl).tree.add_kid(new_tree);
                        end
                    else
                        %% Phenotype Case or Neighborhood
                        valid_name = Helper.valid_neigh(name);
                        if all(startsWith(phenos, Constants.neigh_tag))
                            %% Neighborhoods
                            if iscell(valid_name)
                                valid_name = valid_name{1};
                            end

                            if isfield(app.data.(smpl), 'MFIRSN')
                                logical_rsn = zeros(size(app.data.(smpl).MFIRSN, 1), 1);
                            end

                            if isfield(app.data.(smpl), 'MFICCN')
                                logical_ccn = zeros(size(app.data.(smpl).MFICCN, 1), 1);
                            end

                            for phn_idx=1:numel(phenos)
                                phn = phenos{phn_idx};

                                if in_rsn(phn_idx)
                                    logical_rsn = logical_rsn | app.data.(smpl).MFIRSN.(phn);
                                else
                                    logical_ccn = logical_ccn | app.data.(smpl).MFICCN.(phn);
                                end

                                if any(logical_rsn)
                                    INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                                    INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).MFIRSN.Properties.VariableNames);
                                    INDIgnore = dat.Properties.VariableNames(INDIgnore);
                                    logical_rsn(logical_rsn == 1) = ismember(smpl_dat{phn_idx}, dat(Logical, INDIgnore));
                                    app.data.(smpl).MFIRSN.(valid_name) = logical_rsn;
                                elseif any(logical_ccn)
                                    INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                                    INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).MFICCN.Properties.VariableNames);
                                    INDIgnore = dat.Properties.VariableNames(INDIgnore);
                                    logical_ccn(logical_ccn == 1) = ismember(smpl_dat{phn_idx}, dat(Logical, INDIgnore));
                                    app.data.(smpl).MFICCN.(valid_name) = logical_ccn;
                                end

                            end

                            continue;
                        else
                            %% Phenotype
                            for phn_idx=1:numel(phenos)
                                %deal with some issues with mismatch between channel
                                %indeces between samples
                                INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                                INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, smpl_dat{phn_idx}.Properties.VariableNames);
                                taginclude = dat.Properties.VariableNames(INDIgnore);
                                for taginclude_i = 1:(numel(taginclude)-3)
                                    if ~any(ismember(smpl_dat{phn_idx}(:,taginclude), dat(Logical,taginclude)))
                                        taginclude = taginclude(1:end-1);
                                    elseif any(ismember(smpl_dat{phn_idx}(:,taginclude), dat(Logical,taginclude)))
                                        break 
                                    end
                                end
                                logical_smpl = ismember(smpl_dat{phn_idx}(:,taginclude), dat(Logical, taginclude));

                                phn = phenos{phn_idx};
                                short_name = split(phn, '/');
                                if strcmp(short_name(2), name)
                                    continue;
                                end
                                short_name = strcat(short_name(2), '/', name);
                                if iscell(short_name)
                                    short_name = short_name{1};
                                end
                                tag = Helper.get_tag(app, short_name, smpl);
                                if iscell(tag)
                                    tag = tag{1};
                                end

                                new_path = Helper.gate_full2path(app, phn, smpl);
                                new_path = strcat(new_path, '/', valid_name);
                                if iscell(new_path)
                                    new_path = new_path{1};
                                end

                                new_tree = tree(name, 'gate_points', points, ...
                                    'gate_axes', {app.NewGate.(valid_gate_name).GateAxis1; ...
                                                  app.NewGate.(valid_gate_name).GateAxis2}, ...
                                    'gate_type', lower(app.NewGate.(valid_gate_name).type), ...
                                    'tag', tag ...
                                );

                                parent_tag = Helper.gate_full2tag(app, phn, smpl);
                                if iscell(parent_tag)
                                    parent_tag = parent_tag{1};
                                end
                                new_logical = app.data.(smpl).AllCells.(parent_tag);
                                new_logical(new_logical == 1) = logical_smpl;
                                app.data.(smpl).AllCells.(tag) = new_logical;
                                app.data.(smpl).GateTags.(tag) = {new_path; short_name};
                                app.data.(smpl).tree = app.data.(smpl).tree.add_kid(new_tree, new_path);
                            end
                        end
                    end
                end
                close(wb);
            end

        end

        function func_annotate_clustered_cells(app, UIfig, mdl, ClassChNM)
            %% Build a user interface
            % annotate clusters function builds a user interface for the user to
            % define the names of clustered cells
            %
            % Input:
            %   - app - Instance of CytoMAP


            if nargin < 3
                % Build main UI window
                UIfig = uifigure('Name', 'Annotate Cells', 'Scrollable', 'on');
                Helper.func_SetCLR(app, UIfig, 'UIfigure');
                UIfig.Position = [10 10 320 700];
            end
            
            % Select clustering Model
            lbl = uilabel(UIfig); lbl.Text = sprintf('Select Clustering Model:');
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = [5 700-20 400 15];
            ModelNames = uidropdown(UIfig);
            fnms = fieldnames(app.net);
            ModelNames.Items = {fnms{:}, 'Use Classifier'};
            if nargin > 1
                ModelNames.Value = mdl;
            else
                ModelNames.Value = ModelNames.Items{1};
            end
            ModelNames.Position = [5 700-30-20 175 30];
            ModelNames.BackgroundColor = app.GUIOPTS.bgclr;
            nmstmp = [];
            if strcmp('Use Classifier', ModelNames.Value)
                smpnms = fieldnames(app.data);
                CellNamesTMP = struct;
                CellNamesTMP.N = 1;
                if nargin < 3
                    VNamestmp = Helper.get_channels(app);
                    indx = listdlg('PromptString', 'Select Classification Channel:','SelectionMode','single', 'ListString', VNamestmp);
                    ClassChNM = VNamestmp{indx};
                end
                
                for smp_i = 1:numel(smpnms)
                    if contains(ClassChNM, fieldnames(app.data.(smpnms{smp_i}).AllCells))
                        nmstmp = unique([nmstmp; app.data.(smpnms{smp_i}).AllCells.(ClassChNM)]);
                        numtmp = numel(nmstmp);
                        CellNamesTMP.N = max(CellNamesTMP.N, numtmp);
                    end
                end
                if CellNamesTMP.N > 50
                    prompt = {'More than 50 unique classes detected. Confirm number of cell types:'};
                    dlgtitle = 'Input';
                    dims = [1 35];
                    definput = {num2str(CellNamesTMP.N)};
                    answer = inputdlg(prompt,dlgtitle,dims,definput);
                    if isempty(answer)
                        CellNamesTMP.N = 1;
                    else
                        CellNamesTMP.N = str2double(answer{1});
                    end
                end
                CellNamesTMP.names = struct;
                CellNamesTMP.clusters = struct;
                CellNamesTMP.btnrmv = struct;
                CellNamesTMP.btnclr = struct;
            else
                if ismember('CellNames', fieldnames(app.net.(ModelNames.Value)))
                    CellNamesTMP = app.net.(ModelNames.Value).CellNames;
                    if ismember('N', fieldnames(CellNamesTMP))
                        CellNamesTMP.N = app.net.(ModelNames.Value).NR;
                    end
                    CellNamesTMP.names = struct;
                    CellNamesTMP.clusters = struct;
                    CellNamesTMP.btnrmv = struct;
                    CellNamesTMP.btnclr = struct;
                else
                    % Make a local copy of these tied to this figure
                    CellNamesTMP = struct;
                    CellNamesTMP.N = app.net.(ModelNames.Value).NR;
                    CellNamesTMP.names = struct;
                    CellNamesTMP.clusters = struct;
                    CellNamesTMP.btnrmv = struct;
                    CellNamesTMP.btnclr = struct;
                end
            end

            % Assign Cell Names
            if strcmp('Use Classifier', ModelNames.Value)
                lbl = uilabel(UIfig); lbl.Text = ...
                    sprintf(['Assign Cell Names:                   ' ClassChNM  ' Values:']);
            else
                lbl = uilabel(UIfig); lbl.Text = ...
                    sprintf('Assign Cell Names to Clusters:        Clusters: 1,3, etc.');
            end
            Helper.func_SetCLR(app, lbl, 'label')
            lbl.Position = [5 700-30-2*20 400 15];

            for btni = 1:CellNamesTMP.N
                % assign first name                
                nme = ['name_' num2str(btni)];
                CellNamesTMP.names.(nme) = uidropdown(UIfig,'Editable','on');
                CellNamesTMP.names.(nme).Position = [5 700-30-3*20-20*(btni-1) 175 20];
                if strcmp('Use Classifier', ModelNames.Value)
                    if isempty(nmstmp)
                        CellNamesTMP.names.(nme).Items = {nme};
                        CellNamesTMP.names.(nme).Value = nme;
                    else
                        if iscell(nmstmp)
                            CellNamesTMP.names.(nme).Items = nmstmp;
                            CellNamesTMP.names.(nme).Value = nmstmp{btni};
                        else
                            CellNamesTMP.names.(nme).Items = {nme};
                            CellNamesTMP.names.(nme).Value = nme;
                        end

                    end  
                else
                    if ismember('RegNames', fieldnames(app.net.(ModelNames.Value)))
                        CellNamesTMP.names.(nme).Items = app.net.(ModelNames.Value).RegNames(btni);
                        CellNamesTMP.names.(nme).Value = app.net.(ModelNames.Value).RegNames{btni};
                    else
                        CellNamesTMP.names.(nme).Items = {nme};
                        CellNamesTMP.names.(nme).Value = nme;
                    end
                end

                % assign clusters to this name
                clust = ['clust_' num2str(btni)];
                CellNamesTMP.clusters.(clust) = uieditfield(UIfig,'text');
                CellNamesTMP.clusters.(clust).Position = [5+175 700-30-3*20-20*(btni-1) 80 20];
                
                if isempty(nmstmp)
                    CellNamesTMP.clusters.(clust).Value = num2str(btni);
                else
                    if iscell(nmstmp)
                        CellNamesTMP.clusters.(clust).Value = nmstmp{btni};
                    else
                        CellNamesTMP.clusters.(clust).Value = num2str(btni);
                    end
                    
                end

            end

            % add a remove name button
            CellNamesTMP.btnrmv.btnrmv_1 = uibutton(UIfig,'push');
            CellNamesTMP.btnrmv.btnrmv_1.Position = [5+175+80, 700-30-3*20-20*(btni-1), 20, 20];
            CellNamesTMP.btnrmv.btnrmv_1.Text = '-';
            Helper.func_SetCLR(app, CellNamesTMP.btnrmv.btnrmv_1, 'button')

            % Save Annotations
            svbtn = uibutton(UIfig,'push');
            svbtn.Position = [10+175 700-30-20 100 30];
            svbtn.Text = 'Save Annotations';
            Helper.func_SetCLR(app, svbtn, 'button')

            % add a new name
            btn = uibutton(UIfig,'push');
            btn.Position = [5, 700-30-4*20-20*(btni-1), 20, 20];
            btn.Text = '+';
            Helper.func_SetCLR(app, btn, 'button')
            if strcmp('Use Classifier', ModelNames.Value)
                svbtn.ButtonPushedFcn = @(~,~) func_save_annotations(app, ModelNames, ClassChNM);
                btn.ButtonPushedFcn = @(~,~) func_new_name(app,ModelNames, UIfig, btn, svbtn, ClassChNM);
            else
                svbtn.ButtonPushedFcn = @(~,~) func_save_annotations(app, ModelNames);
                btn.ButtonPushedFcn = @(~,~) func_new_name(app,ModelNames, UIfig, btn, svbtn);
            end
            
            CellNamesTMP.btnrmv.btnrmv_1.ButtonPushedFcn = ...
                @(~,~) func_btnrmv(app,ModelNames, UIfig, btn, svbtn, nme, clust, CellNamesTMP.btnrmv.btnrmv_1);
            ModelNames.ValueChangedFcn = @(btn, event) ChangeDataType(app, ModelNames, UIfig);



            function func_new_name(app, ModelNames, UIfig, btn, svbtn, ClassChNM)
                CellNamesTMP.N = CellNamesTMP.N+1;

                nms = strrep(fieldnames(CellNamesTMP.names), 'name_', '');
                S = sprintf('%s ', nms{:});
                numi = num2str(max(sscanf(S, '%f')+1));

                name_i = ['name_' numi];
                CellNamesTMP.names.(name_i) = uidropdown(UIfig,'Editable','on');
                CellNamesTMP.names.(name_i).Position = [5 btn.Position(2) 175 20];
                CellNamesTMP.names.(name_i).Items ...
                    = fieldnames(CellNamesTMP.names);
                CellNamesTMP.names.(name_i).Value = name_i;

                % assign clusters to this name
                clust_i = ['clust_' numi];
                CellNamesTMP.clusters.(clust_i) = uieditfield(UIfig,'text');
                CellNamesTMP.clusters.(clust_i).Position = [5+175 btn.Position(2) 80 20];
                CellNamesTMP.clusters.(clust_i).Value = '1';

                CellNamesTMP.btnrmv.btnrmv_1.Position(2) = ...
                    CellNamesTMP.btnrmv.btnrmv_1.Position(2)-20;

                CellNamesTMP.btnrmv.btnrmv_1.ButtonPushedFcn = ...
                    @(~,~) func_btnrmv(app, ModelNames, UIfig, btn, svbtn, name_i, clust_i, CellNamesTMP.btnrmv.btnrmv_1);

                btn.Position(2) = btn.Position(2)-20;

% % %                 btn.ButtonPushedFcn = @(~,~) func_new_name(app,ModelNames, UIfig, btn, svbtn);
% % %                 svbtn.ButtonPushedFcn = @(~,~) func_save_annotations(app, ModelNames);
                if strcmp('Use Classifier', ModelNames.Value)
                    svbtn.ButtonPushedFcn = @(~,~) func_save_annotations(app, ModelNames, ClassChNM);
                    btn.ButtonPushedFcn = @(~,~) func_new_name(app,ModelNames, UIfig, btn, svbtn, ClassChNM);
                else
                    svbtn.ButtonPushedFcn = @(~,~) func_save_annotations(app, ModelNames);
                    btn.ButtonPushedFcn = @(~,~) func_new_name(app,ModelNames, UIfig, btn, svbtn);
                end

            end

            function ChangeDataType(app, ModelNames, UIfig)
                if strcmp('Use Classifier', ModelNames.Value)
                    VNamestmp = Helper.get_channels(app);
                    indx = listdlg('PromptString', 'Select Classification Channel:','SelectionMode','single', 'ListString', VNamestmp);
                    ClassCh = VNamestmp{indx};
                    
                    mnm = ModelNames.Value;
                    h=findall(UIfig);
                    delete(h(2:end))
                    Helper.func_annotate_clustered_cells(app, UIfig, mnm, ClassCh)
                        
                else
% % %                     if ismember('CellNames', fieldnames(app.net.(ModelNames.Value)))
                        mnm = ModelNames.Value;
                        h=findall(UIfig);
                        delete(h(2:end))
                        Helper.func_annotate_clustered_cells(app, UIfig, mnm)
% % %                     else
                        % if there are no current cell names do nothing
% % %                     end
                end

            end

            function func_btnrmv(app, ModelNames, UIfig, btn, svbtn, name_i, clust_i, btnrmv_i)
                %%
                %get current name
                nms = fieldnames(CellNamesTMP.names);
                cls = fieldnames(CellNamesTMP.clusters);
                nmind = find(ismember(nms,name_i));

                CellNamesTMP.N = CellNamesTMP.N-1;
                % turn off buttons to be removed
                CellNamesTMP.names.(name_i).Visible = 'off';
                CellNamesTMP.clusters.(clust_i).Visible = 'off';

                btnrmv_i.Position(2) = ...
                    btnrmv_i.Position(2)+20;

                btnrmv_i.ButtonPushedFcn = ...
                    @(~,~) func_btnrmv(app,ModelNames, UIfig, btn, svbtn, nms{end-1}, cls{end-1}, btnrmv_i);


                %% Remove the selected name and clusters

                CellNamesTMP.names = ...
                    rmfield(CellNamesTMP.names, name_i);
                CellNamesTMP.clusters = ...
                    rmfield(CellNamesTMP.clusters, clust_i);

                % move the + button up one
                btn.Position(2) = btn.Position(2)+20;
                %% move other buttons up
                if nmind ~= numel(nms)
                    for nm_i = 1:(numel(nms)-nmind)
                        if contains((['name_' num2str(nm_i + nmind)]), fieldnames(CellNamesTMP.names))
                            if strcmp('Use Classifier', ModelNames.Value)
                                CellNamesTMP.names.(['name_' num2str(nm_i + nmind)]).Position(2) = ...
                                 CellNamesTMP.names.(['name_' num2str(nm_i + nmind)]).Position(2)+20; 
                                CellNamesTMP.clusters.(['clust_' num2str(nm_i + nmind)]).Position(2) = ...
                                 CellNamesTMP.clusters.(['clust_' num2str(nm_i + nmind)]).Position(2)+20;
                                CellNamesTMP.btnrmv.(['btnrmv_' num2str(nm_i + nmind)]).Position(2) = ...
                                 CellNamesTMP.btnrmv.(['btnrmv_' num2str(nm_i + nmind)]).Position(2)+20;
                            else
                                if ~isempty(app.net.(ModelNames.Value).CellNames.names.(['name_' num2str(nm_i + nmind)]))                        
                                    CellNamesTMP.names.(['name_' num2str(nm_i + nmind)]).Position(2) = ...
                                     CellNamesTMP.names.(['name_' num2str(nm_i + nmind)]).Position(2)+20; 
                                    CellNamesTMP.clusters.(['clust_' num2str(nm_i + nmind)]).Position(2) = ...
                                     CellNamesTMP.clusters.(['clust_' num2str(nm_i + nmind)]).Position(2)+20;
                                    CellNamesTMP.btnrmv.(['btnrmv_' num2str(nm_i + nmind)]).Position(2) = ...
                                     CellNamesTMP.btnrmv.(['btnrmv_' num2str(nm_i + nmind)]).Position(2)+20;
                                end
                                
                            end
                        end
                    end
                end
            end

            function func_save_annotations(app, ModelNames, ClassChNM)
                if ~strcmp('Use Classifier', ModelNames.Value)
                    app.net.(ModelNames.Value).CellNames = CellNamesTMP;
                end
                
                for cell_i = 1:CellNamesTMP.N
                    name_i = ['name_' num2str(cell_i)];
                    clust_i = ['clust_' num2str(cell_i)];

                    CellNames = CellNamesTMP.names.(name_i).Value;
                    Clusters = str2double(strsplit(CellNamesTMP.clusters.(clust_i).Value, ','));
                    rng = false;
                    if isnan(Clusters)
                        Clusters = str2double(strsplit(CellNamesTMP.clusters.(clust_i).Value, '-'));
                        rng = true;
                        if isnan(Clusters)
                            % this means it is clusters string (or the wrong format)
                            Clusters = CellNamesTMP.clusters.(clust_i).Value;
                        end
                    end                    
                    if isempty(CellNames)
                        CellNames = ' ';
                    end
                    if strcmp('Use Classifier', ModelNames.Value)
                        save_annotation(app, ModelNames, CellNames, Clusters, ClassChNM, rng)
                    else
                        save_annotation(app, ModelNames, CellNames, Clusters,[], rng)
                    end
                end

            end

            function save_annotation(app, ModelNames, name, Clusters, ClassChNM, rng)
                % save_annotation It is a function that saves the annotations as new cell populations.
                %
                % Input:
                %   - app - Instance of CytoMAP
                %   - num - Identifier of figure from which to pull gate and
                %       data which will be used in saving process.
                %
                % Modifies:
                %   - app - Specifically fields in app.data.(smpl), where smpl
                %       are all samples on which this gate was defined.\
                
                valid_gate_name = Helper.valid_gate(name);
                if iscell(valid_gate_name)
                    valid_gate_name = valid_gate_name{1};
                end

                wb = waitbar(0, 'Please wait. Adding Cell Annotations.');

                % import sample names that contain the chosen model
                smpls = fieldnames(app.data);
                indrmv = zeros(numel(smpls), 1);
                for smpl_i = 1:numel(smpls)
                    if strcmp('Use Classifier', ModelNames.Value)
                        if ~contains(ClassChNM, fieldnames(app.data.(smpls{smpl_i}).AllCells))
                            indrmv(smpl_i) = 1;
                        end
                    else
                        if ~contains([Constants.other_tag ModelNames.Value], fieldnames(app.data.(smpls{smpl_i}).AllCells))
                            indrmv(smpl_i) = 1;
                        end
                    end
                end
                smpls(indrmv==1) = [];

                type = 'Histogram';
                phenos = {'All Cells'};

                % Precompute stuff for non-global neighborhoods (ones, starting with a tag).
                in_rsn = true(numel(phenos), 1);
                if all(startsWith(phenos, Constants.neigh_tag))
                    % Determine if it's MFI RSN or CCN. Consistently for all Samples.
                    for phn_idx = 1:numel(phenos)
                        phn = phenos{phn_idx};
                        for smpl_idx = 1:numel(smpls)
                            smpl = smpls{smpl_idx};
                            if ~isfield(app.data.(smpl), 'MFIRSN') || ...
                                    ~ismember({phn}, app.data.(smpl).MFIRSN.Properties.VariableNames)
                                in_rsn(phn_idx) = false;
                            end
                        end
                    end
                end

                waitbar(.05, wb, 'Saving to the samples');
                for smpl_idx=1:numel(smpls)

                    smpl = smpls{smpl_idx};
                    if numel(phenos) == 1 && strcmp(phenos, 'Density/MFI RSN')
                        smpl_dat = app.data.(smpl).MFIRSN;
                    elseif numel(phenos) == 1 && strcmp(phenos, 'Density/MFI CCN')
                        smpl_dat = app.data.(smpl).MFICCN;
                    elseif numel(phenos) == 1 && strcmp(phenos, 'All Cells')
                        smpl_dat = app.data.(smpl).AllCells;
                    elseif all(startsWith(phenos, Constants.neigh_tag))
                        smpl_dat = {};
                        for phn_idx = 1:numel(phenos)
                            phn = phenos{phn_idx};
                            % Load all the neighborhoods into dat.
                            if in_rsn(phn_idx)
                                smpl_dat{phn_idx} = app.data.(smpl).MFIRSN(app.data.(smpl).MFIRSN.(phn), :);
                            else
                                smpl_dat{phn_idx} = app.data.(smpl).MFICCN(app.data.(smpl).MFICCN.(phn), :);
                            end
                        end
                    else
                        smpl_dat = {};
                        % Pull the data that is to be plotted
                        for phn_idx = 1:numel(phenos)
                            phn = phenos{phn_idx};
                            if ismember({phn}, {'Density/MFI RSN', 'Density/MFI CCN', 'All Cells'}) || ...
                                startsWith({phn}, Constants.neigh_tag)
                                % Just a sanity check it shouldn't be possible to come in here.
                                warndlg("Something went wrong. Neighborhoods and phenotypes shouldn't be mixed. You might need to restart CytoMap.");
                                continue;
                            end
                            tag = Helper.gate_full2tag(app, phn, smpl);
                            if iscell(tag)
                                tag = tag{1};
                            end
                            smpl_dat{phn_idx} = app.data.(smpl).AllCells(app.data.(smpl).AllCells.(tag)==1, :);
                        end
                    end

                    switch type
                        case 'Histogram'
                            points = Clusters;
                            if strcmp('Use Classifier', ModelNames.Value)
                                if rng
                                    if iscell(smpl_dat.(ClassChNM))
                                        logical_smpl = strcmp(smpl_dat.(ClassChNM), Clusters);
                                    else
                                        logical_smpl = min(Clusters) <= smpl_dat.(ClassChNM) ...
                                            & smpl_dat.(ClassChNM) <= max(Clusters);
                                    end
                                else
                                    if iscell(smpl_dat.(ClassChNM))
                                        logical_smpl = zeros(size(smpl_dat.(ClassChNM)))==0;
                                    else
                                        logical_smpl = smpl_dat.(ClassChNM) == Clusters;
                                        logical_smpl = sum(logical_smpl, 2)~=0;
                                    end
                                end
                            else
                                if rng
                                   Clusters = [min(Clusters):1:max(Clusters)];
                                end
                                logical_smpl = ismember(smpl_dat.([Constants.other_tag ModelNames.Value]),Clusters);
                            end
                    end

                    if ~iscell(smpl_dat)
                        if numel(phenos) == 1 && strcmp(phenos, 'Density/MFI RSN')
                            % Add the gated data to the sample
                            valid_name = Helper.valid_neigh(name);
                            if iscell(valid_name)
                                valid_name = valid_name{1};
                            end
                            app.data.(smpl).MFIRSN.(valid_name) = logical_smpl;

                            continue;
                        elseif numel(phenos) == 1 && strcmp(phenos, 'Density/MFI CCN')
                            % Add the gated data to the sample
                            valid_name = Helper.valid_neigh(name);
                            if iscell(valid_name)
                                valid_name = valid_name{1};
                            end
                            app.data.(smpl).MFICCN.(valid_name) = logical_smpl;
                            continue;
                        elseif numel(phenos) == 1 && strcmp(phenos, 'All Cells')
                            short_name = strcat('All/', name);
                            if iscell(short_name)
                                short_name = short_name{1};
                            end
                            tag = Helper.get_tag(app, short_name, smpl);
                            if iscell(tag)
                                tag = tag{1};
                            end
                            path = strcat('Gate_All/', Helper.valid_gate(name));
                            if iscell(path)
                                path = path{1};
                            end
                            if strcmp('Use Classifier', ModelNames.Value)
                                new_tree = tree(name, 'gate_points', points, ...
                                    'gate_axes', {ClassChNM}, ...
                                    'gate_type', lower(type), ...
                                    'tag', tag);
                            else
                                new_tree = tree(name, 'gate_points', points, ...
                                    'gate_axes', {[Constants.other_tag ModelNames.Value]}, ...
                                    'gate_type', lower(type), ...
                                    'tag', tag);
                            end

                            app.data.(smpl).AllCells.(tag) = logical_smpl;
                            app.data.(smpl).GateTags.(tag) = {path; short_name};
                            app.data.(smpl).tree = app.data.(smpl).tree.add_kid(new_tree);
                        end
                    else
                        %% Phenotype Case or Neighborhood
                        valid_name = Helper.valid_neigh(name);
                        if all(startsWith(phenos, Constants.neigh_tag))
                            %% Neighborhoods
                            if iscell(valid_name)
                                valid_name = valid_name{1};
                            end

                            if isfield(app.data.(smpl), 'MFIRSN')
                                logical_rsn = zeros(size(app.data.(smpl).MFIRSN, 1), 1);
                            end

                            if isfield(app.data.(smpl), 'MFICCN')
                                logical_ccn = zeros(size(app.data.(smpl).MFICCN, 1), 1);
                            end

                            for phn_idx=1:numel(phenos)
                                phn = phenos{phn_idx};

                                if in_rsn(phn_idx)
                                    logical_rsn = logical_rsn | app.data.(smpl).MFIRSN.(phn);
                                else
                                    logical_ccn = logical_ccn | app.data.(smpl).MFICCN.(phn);
                                end

                                if any(logical_rsn)
                                    INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                                    INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).MFIRSN.Properties.VariableNames);
                                    INDIgnore = dat.Properties.VariableNames(INDIgnore);
                                    logical_rsn(logical_rsn == 1) = ismember(smpl_dat{phn_idx}, dat(Logical, INDIgnore));
                                    app.data.(smpl).MFIRSN.(valid_name) = logical_rsn;
                                elseif any(logical_ccn)
                                    INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                                    INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, app.data.(smpl).MFICCN.Properties.VariableNames);
                                    INDIgnore = dat.Properties.VariableNames(INDIgnore);
                                    logical_ccn(logical_ccn == 1) = ismember(smpl_dat{phn_idx}, dat(Logical, INDIgnore));
                                    app.data.(smpl).MFICCN.(valid_name) = logical_ccn;
                                end

                            end

                            continue;
                        else
                            %% Phenotype
                            for phn_idx=1:numel(phenos)
                                %deal with some issues with mismatch between channel
                                %indeces between samples
                                INDIgnore = ~startsWith(dat.Properties.VariableNames, Constants.gate_tag);
                                INDIgnore = INDIgnore & ismember(dat.Properties.VariableNames, smpl_dat{phn_idx}.Properties.VariableNames);
                                taginclude = dat.Properties.VariableNames(INDIgnore);
                                for taginclude_i = 1:(numel(taginclude)-3)
                                    if ~any(ismember(smpl_dat{phn_idx}(:,taginclude), dat(Logical,taginclude)))
                                        taginclude = taginclude(1:end-1);
                                    elseif any(ismember(smpl_dat{phn_idx}(:,taginclude), dat(Logical,taginclude)))
                                        break 
                                    end
                                end
                                logical_smpl = ismember(smpl_dat{phn_idx}(:,taginclude), dat(Logical, taginclude));

                                phn = phenos{phn_idx};
                                short_name = split(phn, '/');
                                if strcmp(short_name(2), name)
                                    continue;
                                end
                                short_name = strcat(short_name(2), '/', name);
                                if iscell(short_name)
                                    short_name = short_name{1};
                                end
                                tag = Helper.get_tag(app, short_name, smpl);
                                if iscell(tag)
                                    tag = tag{1};
                                end

                                new_path = Helper.gate_full2path(app, phn, smpl);
                                new_path = strcat(new_path, '/', valid_name);
                                if iscell(new_path)
                                    new_path = new_path{1};
                                end

                                new_tree = tree(name, 'gate_points', points, ...
                                    'gate_axes', {app.NewGate.(valid_gate_name).GateAxis1; ...
                                                  app.NewGate.(valid_gate_name).GateAxis2}, ...
                                    'gate_type', lower(app.NewGate.(valid_gate_name).type), ...
                                    'tag', tag ...
                                );

                                parent_tag = Helper.gate_full2tag(app, phn, smpl);
                                if iscell(parent_tag)
                                    parent_tag = parent_tag{1};
                                end
                                new_logical = app.data.(smpl).AllCells.(parent_tag);
                                new_logical(new_logical == 1) = logical_smpl;
                                app.data.(smpl).AllCells.(tag) = new_logical;
                                app.data.(smpl).GateTags.(tag) = {new_path; short_name};
                                app.data.(smpl).tree = app.data.(smpl).tree.add_kid(new_tree, new_path);
                            end
                        end
                    end
                end
                close(wb);
            end

        end
                
        function [result, model] = func_cluster(n_cluster, alg, data, model)
            switch alg
                case 'NN Self Organizing Map'
                    if isempty(model)
                        model = selforgmap([n_cluster 1]);
                        %Don't show the training window
                        model.trainParam.showWindow = false;
                        % Train the Network
                        [model, ~] = train(model, data');
                    end
                    result = vec2ind(model(data'))';
                case 'k-means'
                    model = [];
                    result = kmeans(data, n_cluster);
                case 'Gaussian Distribution Model'
                    if isempty(model)
                        % Generate a gaussian distribution model
                        % I think this is wrong here TODO: fit it
                        model = fitgmdist(data, n_cluster, 'RegularizationValue', 0.1);
                    end
                    % Implement hard clustering with the model
                    result = cluster(model, data);
                case 'DBSCAN'
                    if isempty(model)
                        model = cell(2,2);
                        % Options for DBSCAN
                        model(:,1) = {'epsilon', ...
                            'Minimum Points within epsilon'};

                        % Defaults
                        model(:,2) = {'2', ...
                            '10'};

                        dlg_title = 'DBSCAN Options';
                        num_lines = 1;
                        vAnswer = inputdlg(model(:,1),dlg_title,num_lines,model(:,2));
                        if isempty(vAnswer)
                            return
                        end
                        model(:,2) = vAnswer;
                    end
                    result = dbscan(data,...
                        str2double(model{1, 2}), ...
                        str2double(model{2, 2}));
            end            
        end
        
        function func_closeVPD(app, vPD)

            if ~isvalid(vPD)
                return;
            end
            switch app.GUI.Dropdowns.FileMenu.Pref_Sound.Text
                case 'Sound ''off'''

                    
                case 'Sound ''on'''
                    if sempty(app.GUIOPTS.dat_audio)
                        app.GUIOPTS.dat_audio = [sin(0.2*(1:1000)), ...
                                                zeros(1,100), ...
                                                sin(0.5*(1:1000)), ...
                                                zeros(1,100), ...
                                                sin(0.8*(1:1000)) ...
                                                    ];
                    end
                    sound(app.GUIOPTS.dat_audio) 
                    
            end  
            close(vPD)
            
% % %             %% Save audio file
% % %             fnm = 'C:\Users\calebst\Desktop\Local Temp\audiofile.wav'
% % %             Fs = 8192;
% % %             audiowrite(fnm,dat_audio,Fs);
% % %             %% load audio file ... should I?
% % %             [y,Fs] = audioread(fnm)
% % %             sound(y, Fs)
        end
        
        function func_prefs(app)

            switch app.GUI.Dropdowns.FileMenu.Pref_Sound.Text
                case 'Sound ''off'''
                    % Sound is currently off, turn it on
                    app.GUI.Dropdowns.FileMenu.Pref_Sound.Text = 'Sound ''on''';
                    
                case 'Sound ''on'''
                    % Sound is currently on, turn it off
                    app.GUI.Dropdowns.FileMenu.Pref_Sound.Text = 'Sound ''off''';
                    
            end
            
        end
    end % end of methods
end % end of class definition

