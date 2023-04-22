classdef tree
% TREE Class representing a singular gate. Allows for representation of multiple nested gates.

    properties
        kids;  % All children gates of this gate.
        name;  % Name of this gate, which can be used as variable.
        full_name;  % Full name of this gate, as passed in to the constructor.
        gate_points;  % Points describing positioning of the gate.
        gate_type;  % Type of the gate (ex. 'rectangle', 'polygon' etc.)
        gate_axes;  % Dimensions along which this gate is specified.
        tag;  % Tag used in AllCells to represent this gate.
        n_gated;  % Number of cells selected on last gate_cells call
    end

    methods
        function obj = tree(full_name, varargin)
% % %             import Helper.*;

            %Define defualt arguments for varargin (all lowercase!)
            defaults = struct('kids', table, ...
                'gate_points', [], ...
                'gate_axes', [], ...
                'gate_type', 'all', ...
                'tag', '');

            % Check that varargin is of even length
            if mod(length(varargin), 2) ~= 0
                error('new_figure_layout needs propertyName/propertyValue pairs, after app argument')
            end

            % Process key, value.
            for pair = reshape(varargin, 2, [])
                if isfield(defaults, lower(pair{1})) % case invariant
                    defaults.(pair{1}) = pair{2}; % Change value at the key to given value
                else  % Key is not in the defaults
                    error('%s is not a recognized parameter name', pair{1})
                end
            end
            obj.kids = defaults.kids;
            if iscell(full_name)
                full_name = full_name{1};
            end
            obj.full_name = full_name;
            obj.name = Helper.valid_gate(full_name);
            if iscell(obj.name)
                obj.name = obj.name{1};
            end
            obj.gate_points = defaults.gate_points;
            obj.gate_axes = defaults.gate_axes;
            obj.gate_type = defaults.gate_type;
            obj.tag = defaults.tag;
            obj.n_gated = 0.0;
        end

        function obj = remove_kid(obj, paths)
            % REMOVE_KID Removes multiple kids from this tree.
            %
            % Input:
            %   - obj - pointer to this tree
            %   - paths - string, char, cell - Paths of nodes in this tree
            %       which will be removed.
            %
            % Output:
            %   - obj - Removes specific columns from kids, or columns from
            %       those kids of those kids etc.
            
            if ~exist('paths', 'var')
                error('Paths have to be given.');
            elseif ~iscell(paths)
                paths = {paths};
            end

            for idx=1:numel(paths)
                parent = split(paths{idx}, '/');
                kid = parent(end);
                if iscell(kid)
                    kid = kid{1};
                end
                parent = join(parent(1:end-1), '/');

                parent_tree = path2tree(obj, parent);
                parent_tree.kids = removevars(parent_tree.kids, {kid});
                if strcmp(parent, obj.name)
                    obj = parent_tree;
                else
                    obj.add_kid(parent_tree, parent); % Re-add it, since MATLAB copies by value.
                end
            end
        end

        function obj = rename_kid(obj, new_kids, paths)
            % RENAME_KID Renames multiple kids in this tree.
            %
            % Input:
            %   - obj - pointer to this tree
            %   - new_kids - string, char, cell - Same size as paths. New
            %       names of those nodes.
            %   - paths - string, char, cell - default: '' - Paths of nodes
            %       in this tree which will be renamed.
            %
            % Output:
            %   - obj - Renames specific columns from kids fields and those
            %       kids (recursive, can be applied to grandkids etc.)
            
            if ~exist('paths', 'var')
                paths = cell(size(new_kids));
                paths(:) = {''};
            elseif ~iscell(paths)
                paths = {paths};
            end
            if ~iscell(new_kids)
                new_kids = {new_kids};
            end
            for idx=1:numel(paths)
                % Figure out specific names, get parent and old kid.
                valid_kid = Helper.valid_gate(new_kids{idx});
                parent = split(paths{idx}, '/');
                kid = parent(end);
                parent = join(parent(1:end-1), '/');
                parent_tree = path2tree(obj, parent);
                if iscell(kid)
                    kid = kid{1};
                end

                kid_tree = parent_tree.kids.(kid);
                % Remove old kid
                parent_tree.kids = removevars(parent_tree.kids, {kid});

                % Change and read new kid
                kid_tree.name = valid_kid;
                kid_tree.full_name = new_kids{idx};
                parent_tree.kids.(valid_kid) = kid_tree;

                obj = obj.add_kid(parent_tree, parent); % Re-add parent to root, since MATLAB copies by value.
            end
        end
        
        function obj = retag_kid(obj, new_tag, paths)
            % RETAG_KID Given multiple kids in this tree (obj) changes
            % their tag to new ones.
            %
            % Input:
            %   - obj - pointer to this tree
            %   - new_tag - string, char, cell - Same size as paths. New
            %       tags of those nodes.
            %   - paths - string, char, cell - default: '' - Paths of nodes
            %       in this tree which will be retagged.
            %
            % Output:
            %   - obj - Retags specific columns from kids fields and those
            %       kids (recursive, can be applied to grandkids etc.)
            
            if ~exist('paths', 'var')
                paths = cell(size(new_tag));
                paths(:) = {''};
            elseif ~iscell(paths)
                paths = {paths};
            end
            for idx=1:numel(paths)
                % Figure out specific names, get parent and old kid.
                parent = split(paths{idx}, '/');
                kid = parent(end);
                parent = join(parent(1:end-1), '/');
                parent_tree = path2tree(obj, parent);
                if iscell(kid)
                    kid = kid{1};
                end

                % Change Tag
                if isfield(parent_tree.kids, kid)
                    parent_tree.kids.(kid).tag = new_tag{idx};
                else
                    full_kid = Helper.full_gate(kid);
                    parent_tree.kids.(kid) = tree(full_kid, 'tag', new_tag{idx}, 'gate_type', 'logic');
                end

                obj = obj.add_kid(parent_tree, parent); % Re-add parent to root, since MATLAB copies by value.
            end
        end

        function obj = add_kid(obj, kid, paths)
            % ADD_KID Adds new subtrees to a given tree (obj), under given
            % paths.
            %
            % Input:
            %   - obj - pointer to this tree
            %   - kid - tree, cell(trees) - Same size as paths. Subtrees to
            %       be added to this tree.
            %   - paths - string, char, cell - default: '' - Paths of
            %       places to which put new subtrees.
            %
            % Returns:
            %   - obj - Given obj, but with added new subtrees at specified
            %       paths.
            
            
            % Add a kid to the tree. Can be multiple of them in an array/cell.
            % Paths specifies path to that kid (inclusive with kid).
            % If not given defaults to being immidately under obj.
            % There should be 1-to-1 mapping from path to kid.
            if ~exist('paths', 'var')
                paths = cell(size(kid));
                paths(:) = {''};
            elseif ~iscell(paths)
                paths = {paths};
            end
            if ~isa(kid, 'tree')
                error('Kid needs to be a tree');
            else
                for idx=1:numel(kid)
                    i = kid(idx);
                    i.name = Helper.valid_var(i.name);
                    root_child = add_single_kid(i, paths{idx});
                    if strcmp(root_child.name, obj.name)
                        obj = root_child;
                    else
                        obj.kids.(root_child.name) = root_child;
                    end
                end
            end

            % Basically traverses the tree from top to bottom adding a tree on the bottom.
            % Is a workaround since there is no such thing as pointers in MatLab
            % In other words: if we have a tree a -> b -> c, and we want to add 'd' under 'c'.
            % Then this function will recursively iterate such that:
            % 1. c -> d
            % 2. b -> c -> d
            % Then a Caller to this function should override the variable at the root of the tree (as it does)
            % NOTE, SPECIAL CASE:
            % if a child is added right under root, then it returns root, with that added child.
            function tmp = add_single_kid(kid_tmp, path)
                if contains(path, '/')
                    new_path = split(path, '/');
                    new_path = new_path(1:end - 1);
                    new_path = join(new_path, '/');
                    tmp = obj.path2tree(new_path);
                    tmp.kids.(kid_tmp.name) = kid_tmp;

                    if contains(new_path, '/')
                        tmp = add_single_kid(tmp, new_path);
                    end
                else
                    tmp = kid_tmp;
                    return;
                end
            end
        end

        function kid_names = get_kid_names(obj, recurse)
            % GET_KID_NAMES It will return valid names of all the kids of a
            % given tree.
            %
            % Input:
            %   - obj - pointer to this tree
            %   - recurse - bool - defualt: false - whether to return names
            %       only from one generation below (false), or from all
            %       generations below (true).
            %
            % Output:
            %   - kid_names - string array, cell - All name of children
            %       from first generation (string array, if recurse) or all
            %       generations below (cell, if not recurse).
            
            if ~exist('recurse', 'var')
                recurse = false;
            end
            % Get names of all the kids.
            if recurse
                kid_names = {};
                for i=obj.kids
                    kid_names(end + 1) = {i{1, 1}.name};
                    kids_i = i{1, 1}.get_kid_names(recurse);
                    kid_names = [kid_names(:)', strcat(i{1, 1}.name, '/', kids_i(:))'];
                end
            else
                % We can pre-allocate memory
                kid_names = zeros(length(obj.kids));
                j=1;
                for i=obj.kids
                    kid_names(j) = i.name;
                    j=j+1;
                end
            end
        end

        function sub_tree = get_subtree(obj, paths)
            % GET_SUBTREE
            %
            % Input:
            %   - obj - pointer to this tree
            %   - paths - string, char, cell - default: '' - Paths of
            %       subtrees which are to be returned.
            %
            % Output:
            %   - sub_tree - cell of trees - All the subtrees corresponding
            %       to given paths.
            
            to_remove = obj.get_kid_names(true);
            to_remove = strcat(obj.name, '/', to_remove);
            to_remove = to_remove(~ismember(to_remove, paths));

            idx = 1;
            % Get rid of any gates which are children of gates we will have to remove either way
            while idx <= numel(to_remove)
                logic = startsWith(to_remove, to_remove(idx));
                logic(idx) = false;  % Do not remove string itself.

                % Update to remove.
                to_remove = to_remove(~logic);

                % Increment idx by 1.
                % Lower down if any of the elements of the to_remove occuring
                % before the idx were removed.
                idx = idx + 1 - sum(logic(1:idx));
            end

            sub_tree = obj.remove_kid(to_remove);
        end

        function logic = iskid(obj, path)
            % ISKID Returns a bool array of whether there exists a child
            % with given path in obj tree
            
            if ~iscell(path)
                path = {path};
            end
            logic = ismember(path, obj.get_kid_names(true));
        end

        function leaf = is_leaf(obj)
            % IS_LEAF
            % Whether this gate has any kids or not.
            
            leaf = isempty(obj.kids);
        end

        function new_cell_table = gate_cells(obj, cell_table)
            % GATE_CELLS Given a cell table, modifies it and adds logical
            % columns corresponding to the gating of each gate in that tree.
            %
            % Input:
            %   - obj - Pointer to this tree
            %   - cell_table - table - Table of cells to gate on this tree
            %       (i.e. app.data.(smpl).AllCells)
            %
            % Output:
            %   - new_cell_table - table - Given cell_table, but with new
            %       columns corresponding to tags of given tree with
            %       logical values of whether cell is within that gate of
            %       tree or not.
            
            if strcmp(obj.name, 'AllCells')
                for kid=1:numel(obj.kids)
                    recurse_gating(obj.kids.(kid), ones(size(cell_table, 1), 1));
                end
            else
                recurse_gating(obj, ones(size(cell_table, 1), 1));
            end

            new_cell_table = cell_table;
            function recurse_gating(obj, init_logic)
                %{
                    Params:
                        - obj - a child node
                        - init_logic - logic corresponding to gating of a
                                        parent
                    Modifies:
                        - cell_table - which is passed to gate_cells, and
                            then it's modified version is returned.
                %}
                switch obj.gate_type
                    case 'all'
                        Logic = init_logic;
                    case 'ellipsoid'
%                         params = [focis; center; [a, 0]] =  obj.gate_points
                        focis = obj.gate_points(1:(end-1), :);
% % %                         center = obj.gate_points(end-1, :);
                        a = obj.gate_points(end, 1);

                        % Points to gate
                        P = [table2array(cell_table(:, obj.gate_axes(1))), ...
                             table2array(cell_table(:, obj.gate_axes(2)))];

                        % Find the distance to F1 and F2
                        norm1 = sqrt(sum((P-focis(1,:)).^2, 2));
                        norm2 = sqrt(sum((P-focis(2,:)).^2, 2));
                        dist = norm1 + norm2;
                        % if the distance is more than 2*a the point is outside the ellipse
                        Logic = (dist <= 2*a);
% % %                         %Plot the resulting data
% % %                         figure(3)
% % %                         cla
% % %                         hold on
% % %                         plot(focis(:,1),focis(:,2),'mo');
% % %                         plot(P(~Logic, 1), P(~Logic, 2), 'b.');
% % %                         plot(P(Logic, 1), P(Logic, 2), 'r.');
% % %
% % %                         xlim([-10 200])
% % %                         ylim([-10 200])

%                         threshold = 0.05;
%                         a = max([obj.gate_points(1) obj.gate_points(2)]);
%                         b = min([obj.gate_points(1) obj.gate_points(2)]);
%                         sin_alpha = sin(obj.gate_points(3));
%                         cos_alpha = cos(obj.gate_points(3));
%                         x = obj.gate_points(4);
%                         y = obj.gate_points(5);
%
%                         x = table2array(cell_table(:, obj.gate_axes(1))) - x;
%                         y = table2array(cell_table(:, obj.gate_axes(2))) - y;
%
%                         % From the equation of the ellipsoid.
%                         if strcmp(obj.full_name, 'Lymphocytes')
%                             % DEBUG
%                             ((((x .* cos_alpha + y .* sin_alpha) .^ 2) ./ a^2) + ...
%                                 (((x .* sin_alpha - y .* cos_alpha) .^ 2) ./ b^2))
%                         end
%
%                         Logic = ((((x .* cos_alpha + y .* sin_alpha) .^ 2) ./ a^2) + ...
%                             (((x .* sin_alpha - y .* cos_alpha) .^ 2) ./ b^2)) ...
%                             <= 1 + threshold;
                    case 'polygon'
                        [Logic, ~] = inpolygon(table2array(cell_table(:, obj.gate_axes(1))), ...
                            table2array(cell_table(:, obj.gate_axes(2))), ...
                            obj.gate_points(:, 1), ...
                            obj.gate_points(:, 2));
                    case 'rectangle'
                        [Logic, ~] = inpolygon(table2array(cell_table(:, obj.gate_axes(1))), ...
                            table2array(cell_table(:, obj.gate_axes(2))), ...
                            obj.gate_points(:, 1), ...
                            obj.gate_points(:, 2));
                    case 'logic'
                        if ismember(obj.tag, cell_table.Properties.VariableNames)
                            Logic = table2array(cell_table(:, obj.tag));
                        elseif ~isempty(obj.gate_points)
                            Logic = table2array(obj.gate_points);
                            obj.gate_points = [];
                        else
                            Logic = 0 .* cell_table.X;
                        end

                    otherwise
                        error('Gate not recognized, or not yet implemented.')
                end
                Logic = Logic & init_logic; % So child doesn't contain more than parent.
                obj.n_gated = sum(Logic);
                cell_table.(obj.tag) = Logic;
                for kids_sub=1:numel(obj.kids)
                    recurse_gating(obj.kids.(kids_sub), Logic);
                end
            end
        end

        function tree_from_path = path2tree(obj, path)
            % PATH2TREE Given an absolute path, returns a tree that path is
            % pointing to.
            %
            % Input:
            %   - obj - Pointer to this tree
            %   - path - string, char, cell, table - Single path to the
            %       node which will be returned as a tree.
            %
            % Output:
            %   - tree_from_path - tree - Tree which given path was
            %       pointing to.
            %
            % Note:
            %   This function is basically identical to the get_subtree,
            %   but which one is faster is not determined, so both are
            %   kept. Additionally this function only allows for ONE path,
            %   while get_subtree can take care of multiple subtrees.
            

            % Format path correctly (to string)
            if isa(path, 'table')
                path = table2cell(path);
            end
            if iscell(path)
                path = path{1};
            end

            % If path starts with name of current node, then ignore it
            if startsWith(path, obj.name)
                path = path(numel(obj.name) + 1:end);
            end

            % Call recursively until there is no '/', and then remove:
            %  a child path is point to, or itself if path is empty
            if contains(path, '/')
                [child_name, new_path] = strtok(path, '/');
                new_path = new_path(2:end);
                try
                    tree_from_path = path2tree(obj.kids.(child_name), new_path);
                catch
                end
            elseif ~isempty(path)
                tree_from_path = obj.kids.(path);
            else
                tree_from_path = obj;
            end
        end

        function ui_tree = tree2uitree(obj, fig)
            % TREE2UITREE Returns a uitree representation of that tree.
            % It's nested under given fig. For more info on nodes look for
            % documentation of tree2uitreenode.
            
            ui_tree = uitree(fig);
            obj.tree2uitreenode(ui_tree);
        end

        function ui_node = tree2uitreenode(obj, ui_tree, parent_gated)
            % TREE2UITREENODE Returns a uitreenode representation of that tree.
            % Tag is included in NodeData property of uitreenode.
            
            if ~exist('parent_gated', 'var')
                parent_gated = obj.n_gated;
            end

            if parent_gated == 0
                ui_node = uitreenode(ui_tree, 'Text', obj.name, 'NodeData', [obj.n_gated 0.0 obj.tag]);
            else
                ui_node = uitreenode(ui_tree, 'Text', obj.name, 'NodeData', [obj.n_gated 1.0*obj.n_gated/parent_gated obj.tag]);
            end
            for i=obj.kids
                i{1, 1}.tree2uitreenode(ui_node, obj.n_gated);
            end
        end

        function names = get_full_name_kids(obj)
            % GET_FULL_NAME_KIDS Returns all full names of kids from all generations below of
            % a given tree (i.e. similar to get_kid_names with
            % recurse=true, and full names instead of valid ones). 
            
            names = {};
            for i=obj.kids
                names = horzcat(names, priv_recurse(i{1, 1}));
            end
            names = names';

            function names = priv_recurse(obj)
                names = {strcat(obj.full_name, ' (', obj.tag, ')')};
                for j=obj.kids
                    names = horzcat(names, priv_recurse(j{1, 1}));
                end
            end
        end
    end
end
