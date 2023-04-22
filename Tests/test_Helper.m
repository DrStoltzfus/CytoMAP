% Test Helper - Start Settings
classdef test_Helper
    methods (Static)
        function Import_Definitions_Func
            addpath('../CytoMAP/');
            addpath('../CytoMAP/3rdPartyFunctions');
            addpath('../Tests/');

            import Helper.*;
            import CytoMAP.*;
            import data_loader.*;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Other, more basic functions; or singular functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % TODO
        % function test_SurfCreate()
        %     %% Init cytomap, load data
        %     app = CytoMAP;
        %     data_loader.load_2d_dataset_mat(app);
        %
        %     %% 2D Dataset, normal name
        %
        %     %% 2D Dataset, weird name
        %
        %     %% Close the CytoMAP instance
        %     Helper.func_exit(app);
        %
        %     %% Init cytomap, load data
        %     app = CytoMAP;
        %     data_loader.load_2d_dataset_mat(app);
        %
        %     % 3D Dataset, normal name
        %
        %     % 3D Dataset, weird name
        %
        %     %% Close the CytoMAP instance
        %     Helper.func_exit(app);
        % end

        function test_setequal()
            a = {'a', 'b', 'c', 'd'};
            %% Test 1: Equal sets, different cases
            b = {'a', 'b', 'c', 'd'};
            assert(Helper.setequal(a, b));

            b = {'b', 'a', 'd', 'c'};
            assert(Helper.setequal(a, b));

            b = {"b", 'a', "d", 'c'};
            assert(Helper.setequal(a, b));

            assert(Helper.setequal({}, {}));

            %% Test 2: Set inequalities
            b = {'a', 'b', 'c', 'e'};
            assert(~Helper.setequal(a, b));

            b = {"e", 'a', "d", 'c'};
            assert(~Helper.setequal(a, b));

            assert(~Helper.setequal({}, a));
        end

        function test_logical2cell()
            % Empty cell
            out = Helper.logical2cell([]);
            assert(iscell(out) && isempty(out));

            % Simple case
            in = [1, 0, 0, 1, 0];
            out = Helper.logical2cell(in);
            assert(iscell(out));
            for idx=1:numel(out)
                assert(out{idx} == logical(in(idx)));
            end

            % Convert to logical case
            in = [1, 0, 2, 1, 0];
            out = Helper.logical2cell(in);
            assert(iscell(out));
            for idx=1:numel(out)
                assert(out{idx} == logical(in(idx)));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Processing gates, channels etc. to make it into valid MatLab vars.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function test_gate_processing()
            app = CytoMAP;

            data_loader.load_2d_dataset_mat(app);
            size = numel(app.data.all_1.GateTags{2, :});

            %% full2tag function
            % All Valid Tags, Sample given
            assert(Helper.setequal(...
                Helper.gate_full2tag(app, app.data.all_1.GateTags{2, :}, 'all_1'), ...
                app.data.all_1.GateTags.Properties.VariableNames ...
            ));

            % All Valid Tags, Sample default
            assert(Helper.setequal(...
                Helper.gate_full2tag(app, app.data.no_T4.GateTags{2, :}), ...
                app.data.no_T4.GateTags.Properties.VariableNames ...
            ));

            % Valid Tags + Others, Sample given
            input = [app.data.all_1.AllCells.Properties.VariableNames(1:end - size), app.data.all_1.GateTags{2, :}];
            assert(Helper.setequal(...
                Helper.gate_full2tag(app, input, 'all_1'), ...
                app.data.all_1.AllCells.Properties.VariableNames ...
            ));

            % Some valid Tags + Others, Sample given
            input = [app.data.all_1.AllCells.Properties.VariableNames(1:end - size), app.data.all_1.GateTags{2, 1:end-1}];
            assert(Helper.setequal(...
                Helper.gate_full2tag(app, input, 'all_1'), ...
                app.data.all_1.AllCells.Properties.VariableNames(1:end - 1) ...
            ));

            % No valid tags, Sample given
            assert(Helper.setequal(...
                Helper.gate_full2tag(app, app.data.all_1.AllCells.Properties.VariableNames(1:end - size), 'all_1'), ...
                app.data.all_1.AllCells.Properties.VariableNames(1:end - size) ...
            ));

            assert(Helper.setequal(...
                Helper.gate_full2tag(app, {}), ...
                {} ...
            ));

            %% tag2full function
            % All Valid Tags, Sample given
            assert(Helper.setequal(...
                Helper.gate_tag2full(app, app.data.all_1.GateTags.Properties.VariableNames, 'all_1'), ...
                app.data.all_1.GateTags{2, :} ...
            ));

            % All Valid Tags, Sample default
            assert(Helper.setequal(...
                Helper.gate_tag2full(app, app.data.no_T4.GateTags.Properties.VariableNames), ...
                app.data.no_T4.GateTags{2, :} ...
            ));

            % Valid Tags + Others, Sample given
            output = [app.data.all_1.AllCells.Properties.VariableNames(1:end - size), app.data.all_1.GateTags{2, :}];
            assert(Helper.setequal(...
                Helper.gate_tag2full(app, app.data.all_1.AllCells.Properties.VariableNames, 'all_1'), ...
                output ...
            ));

            % Some valid Tags + Others, Sample given
            output = [app.data.all_1.AllCells.Properties.VariableNames(1:end - size), app.data.all_1.GateTags{2, 1:end-1}];
            assert(Helper.setequal(...
                Helper.gate_tag2full(app, app.data.all_1.AllCells.Properties.VariableNames(1:end - 1), 'all_1'), ...
                output ...
            ));

            % No valid tags, Sample given
            assert(Helper.setequal(...
                Helper.gate_tag2full(app, app.data.all_1.AllCells.Properties.VariableNames(1:end - size), 'all_1'), ...
                app.data.all_1.AllCells.Properties.VariableNames(1:end - size) ...
            ));

            assert(Helper.setequal(...
                Helper.gate_tag2full(app, {}), ...
                {} ...
            ));

            %% full2path
            % All Valid Tags, Sample given
            assert(Helper.setequal(...
                Helper.gate_full2path(app, app.data.all_1.GateTags{2, :}, 'all_1'), ...
                app.data.all_1.GateTags{1, :} ...
            ));

            % All Valid Tags, Sample default
            assert(Helper.setequal(...
                Helper.gate_full2path(app, app.data.no_T4.GateTags{2, :}), ...
                app.data.no_T4.GateTags{1, :} ...
            ));

            % Valid Tags + Others, Sample given
            input = [app.data.all_1.AllCells.Properties.VariableNames(1:end - size), app.data.all_1.GateTags{2, :}];
            output = [app.data.all_1.AllCells.Properties.VariableNames(1:end - size), app.data.all_1.GateTags{1, :}];
            assert(Helper.setequal(...
                Helper.gate_full2path(app, input, 'all_1'), ...
                output ...
            ));

            % Some valid Tags + Others, Sample given
            input = [app.data.all_1.AllCells.Properties.VariableNames(1:end - size), app.data.all_1.GateTags{2, 1:end-1}];
            output = [app.data.all_1.AllCells.Properties.VariableNames(1:end - size), app.data.all_1.GateTags{1, 1:end-1}];
            assert(Helper.setequal(...
                Helper.gate_full2path(app, input, 'all_1'), ...
                output ...
            ));

            % No valid tags, Sample given
            assert(Helper.setequal(...
                Helper.gate_full2path(app, app.data.all_1.AllCells.Properties.VariableNames(1:end - size), 'all_1'), ...
                app.data.all_1.AllCells.Properties.VariableNames(1:end - size) ...
            ));

            assert(Helper.setequal(...
                Helper.gate_full2path(app, {}), ...
                {} ...
            ));

            %% get_tag
            % All Valid Tags, Sample given
            assert(Helper.setequal(...
                Helper.get_tag(app, app.data.all_1.GateTags{2, :}, 'all_1'), ...
                app.data.all_1.GateTags.Properties.VariableNames ...
            ));

            % All Valid Tags, Sample default
            assert(Helper.setequal(...
                Helper.get_tag(app, app.data.no_T4.GateTags{2, :}), ...
                app.data.no_T4.GateTags.Properties.VariableNames ...
            ));

            % Valid Tags + Others, Sample given
            input = [app.data.all_1.AllCells.Properties.VariableNames(1:end - size), app.data.all_1.GateTags{2, :}];
            output = cell(1, numel(app.data.all_1.AllCells.Properties.VariableNames(1:end - size)));
            for allcells_idx = 1:numel(app.data.all_1.AllCells.Properties.VariableNames(1:end - size))
                output(allcells_idx) = {strcat(Constants.gate_tag, num2str(allcells_idx + size))};
            end
            output = [output, app.data.all_1.GateTags.Properties.VariableNames];
            assert(Helper.setequal(...
                Helper.get_tag(app, input, 'all_1'), ...
                output ...
            ));

            % Some valid Tags + Others, Sample given
            input = [app.data.all_1.AllCells.Properties.VariableNames(1:end - size), app.data.all_1.GateTags{2, 1:end-1}];
            output = cell(1, numel(app.data.all_1.AllCells.Properties.VariableNames(1:end - size)));
            for allcells_idx = 1:numel(app.data.all_1.AllCells.Properties.VariableNames(1:end - size))
                output(allcells_idx) = {strcat(Constants.gate_tag, num2str(allcells_idx + size))};
            end
            output = [output, app.data.all_1.GateTags.Properties.VariableNames(1:end-1)];
            assert(Helper.setequal(...
                Helper.get_tag(app, input, 'all_1'), ...
                output ...
            ));

            % No valid tags, Sample given
            output = cell(1, numel(app.data.all_1.AllCells.Properties.VariableNames(1:end - size)));
            for allcells_idx = 1:numel(app.data.all_1.AllCells.Properties.VariableNames(1:end - size))
                output(allcells_idx) = {strcat(Constants.gate_tag, num2str(allcells_idx + size))};
            end
            assert(Helper.setequal(...
                Helper.get_tag(app, app.data.all_1.AllCells.Properties.VariableNames(1:end - size), 'all_1'), ...
                output ...
            ));

            assert(Helper.setequal(...
                Helper.get_tag(app, {}), ...
                {} ...
            ));

            %% valid_gate
            % All weird cases
            input = {"_A", '______A_+_b_', "", "B", 'C'};
            output = cellstr(strcat(Constants.gate_tag, ["A", "A_POS_b_", "", "B", "C"]));
            assert(Helper.setequal(...
                Helper.valid_gate(input), ...
                output ...
            ));

            input = Constants.other_names;
            output = input;
            assert(Helper.setequal(...
                Helper.valid_gate(input), ...
                output ...
            ));

            input = {};
            output = {};
            assert(Helper.setequal(...
                Helper.valid_gate(input), ...
                output ...
            ));

            %% full_gate
            input = {"Gate_A_B_C", 'Gate_vvPOSNEGs', 'not_gate', '', "", "Gate__POS"};
            output = {'A B C', 'vv+-s', 'not_gate', '', '', ' +'};
            assert(Helper.setequal(...
                Helper.full_gate(input), ...
                output ...
            ));

            input = {};
            output = {};
            assert(Helper.setequal(...
                Helper.full_gate(input), ...
                output ...
            ));

            %% Chane Sample loaded to multi-sample one.
            Helper.func_exit(app);
            app = CytoMAP;
            data_loader.load_2d_dataset_mat(app);

            %% get_gates
            % Single Sample Exisiting
            [gates, short_names] = Helper.get_gates(app, 'all_2');
            gates_out = {...
                'Gate_All/Gate_Lymphocytes_CD11cPOS', ...
                'Gate_All/Gate_T_Cells_OTNEGII_4GETPOS', ...
                'Gate_All/Gate_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                'Gate_All/Gate_T_Cells_OT_II_CXCR3POS', ...
                'Gate_All/Gate_T_cells_OTNEGII_T_Cells' ...
            };
            assert(...
                Helper.setequal(gates, gates_out) && ...
                Helper.setequal(short_names, app.data.all_2.GateTags{2, :}) ...
            );

            % Single Sample Non-Existing
            [gates, short_names] = Helper.get_gates(app, 'all_23');
            assert(isempty(gates) && isempty(short_names));

            % Multiple Samples, all existing (same gates)
            [gates, short_names] = Helper.get_gates(app, {'all_1', 'all_2'});
            assert(...
                Helper.setequal(gates, gates_out) && ...
                Helper.setequal(short_names, app.data.all_2.GateTags{2, :}) ...
            );

            % Multiple Sample, all exisiting (different gates)
            [gates, short_names] = Helper.get_gates(app, {'no_T4', 'all_1'});
            gates_out = { ...
                'Gate_All/Gate_Lymphocytes_CD11cPOS', ...
                'Gate_All/Gate_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                'Gate_All/Gate_T_Cells_OT_II_CXCR3POS', ...
                'Gate_All/Gate_T_cells_OTNEGII_T_Cells' ...
            };
            assert(...
                Helper.setequal(gates, gates_out) && ...
                Helper.setequal(short_names, app.data.no_T4.GateTags{2, :}) ...
            );

            % Multiple Sample, multiple exisiting, multiple non-existing (different gates)
            [gates, short_names] = Helper.get_gates(app, {'NE1', 'no_T4', 'all_1', 'Full_csv_23'});
            gates_out = { ...
                'Gate_All/Gate_Lymphocytes_CD11cPOS', ...
                'Gate_All/Gate_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                'Gate_All/Gate_T_Cells_OT_II_CXCR3POS', ...
                'Gate_All/Gate_T_cells_OTNEGII_T_Cells' ...
            };
            assert(...
                Helper.setequal(gates, gates_out) && ...
                Helper.setequal(short_names, app.data.no_T4.GateTags{2, :}) ...
            );

            % Multiple Samples, one existing, multiple non-existing
            [gates, short_names] = Helper.get_gates(app, {'NE', 'all_2', 'NE2'});
            gates_out = {...
                'Gate_All/Gate_Lymphocytes_CD11cPOS', ...
                'Gate_All/Gate_T_Cells_OTNEGII_4GETPOS', ...
                'Gate_All/Gate_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                'Gate_All/Gate_T_Cells_OT_II_CXCR3POS', ...
                'Gate_All/Gate_T_cells_OTNEGII_T_Cells' ...
            };
            assert(...
                Helper.setequal(gates, gates_out) && ...
                Helper.setequal(short_names, app.data.all_2.GateTags{2, :}) ...
            );

            %% get_gate_tags
            % Single Sample Exisiting
            assert(Helper.setequal(...
                Helper.get_gate_tags(app, 'all_2'), ...
                app.data.all_2.GateTags.Properties.VariableNames ...
            ));

            % Single Sample Non-Existing
            assert(isempty(Helper.get_gate_tags(app, 'all_23')));

            % Multiple Samples, all existing (same gates)
            assert(Helper.setequal(...
                Helper.get_gate_tags(app, {'all_2', 'all_1'}), ...
                app.data.all_2.GateTags.Properties.VariableNames ...
            ));


            % Multiple Samples, all exisiting (different gates)
            assert(Helper.setequal(...
                Helper.get_gate_tags(app, {'no_T4', 'all_1'}), ...
                app.data.no_T4.GateTags.Properties.VariableNames(:) ...
            ));
            assert(Helper.setequal(...
                Helper.get_gate_tags(app, {'all_1', 'no_T4'}), ...
                strcat(Constants.gate_tag, {'1', '3', '4', '5'}) ...
            ));

            % Multiple Samples, all existing (no common gates)
            assert(isempty(...
                Helper.get_gate_tags(app, {'no_lymph', 'other'}) ...
            ));

            % Multiple Samples, multiple exisiting, multiple non-existing (different gates)
            assert(Helper.setequal(...
                Helper.get_gate_tags(app, {'NE1', 'no_T4', 'all_2', 'all_23'}), ...
                app.data.no_T4.GateTags.Properties.VariableNames(:) ...
            ));

            % Multiple Samples, one existing, multiple non-existing
                assert(Helper.setequal(...
                Helper.get_gate_tags(app, {'NE', 'all_2', 'NE2'}), ...
                app.data.all_2.GateTags.Properties.VariableNames ...
            ));

            Helper.func_exit(app);
        end

        function test_channel_processing()
            %% valid_channel
            % All weird cases
            input = {"_A", '______A_+_b_', "", "B", 'C'};
            output = cellstr(strcat(Constants.channel_tag, ["A", "A_POS_b_", "", "B", "C"]));
            assert(Helper.setequal(...
                Helper.valid_channel(input), ...
                output ...
            ));

            input = Constants.other_names;
            output = input;
            assert(Helper.setequal(...
                Helper.valid_channel(input), ...
                output ...
            ));

            input = {};
            output = {};
            assert(Helper.setequal(...
                Helper.valid_channel(input), ...
                output ...
            ));

            %% full_channel
            input = [cellstr(strcat(Constants.channel_tag, ["A_B_C", 'vvPOSNEGs', '_POS'])), {'not_channel', ''}];
            output = {'A B C', 'vv+-s', 'not_channel', "", " +"};
            assert(Helper.setequal(...
                Helper.full_channel(input), ...
                output ...
            ));

            input = {};
            output = {};
            assert(Helper.setequal(...
                Helper.full_channel(input), ...
                output ...
            ));

            %% Create CytoMAP instance, and load data
            app = CytoMAP;
            data_loader.load_2d_dataset_mat(app);

            %% get_channels
            % Sample given
            output = app.data.all_2.AllCells.Properties.VariableNames;
            output = output(startsWith(output, Constants.channel_tag) | ismember(output, {'X', 'Y', 'Z'}));
            assert(Helper.setequal(...
                Helper.get_channels(app, 'all_2'), ...
                output ...
            ));

            % Default sample
            output = app.data.(app.DataN.Value).AllCells.Properties.VariableNames;
            output = output(startsWith(output, Constants.channel_tag) | ismember(output, {'X', 'Y', 'Z'}));
            assert(Helper.setequal(...
                Helper.get_channels(app), ...
                output ...
            ));

            % Multiple samples, all common channels
            output = app.data.all_2.AllCells.Properties.VariableNames;
            output = output(startsWith(output, Constants.channel_tag) | ismember(output, {'X', 'Y', 'Z'}));
            assert(Helper.setequal(...
                Helper.get_channels(app, {'all_1', 'all_2'}), ...
                output ...
            ));

            % Multiple samples, no common channels
            output = app.data.all_2.AllCells.Properties.VariableNames;
            output = intersect(output, app.data.other.AllCells.Properties.VariableNames);
            output = output(startsWith(output, Constants.channel_tag) | ismember(output, {'X', 'Y', 'Z'}));
            assert(Helper.setequal(...
                Helper.get_channels(app, {'all_1', 'other'}), ...
                output ...
            ));

            Helper.func_exit(app);
        end

        function test_neigh_processing()
            %% valid_neigh
            % All weird cases
            input = {"_A", '______A_+_b_', "", "B", 'C'};
            output = cellstr(strcat(Constants.neigh_tag, ["A", "A_POS_b_", "", "B", "C"]));
            assert(Helper.setequal(...
                Helper.valid_neigh(input), ...
                output ...
            ));

            input = Constants.other_names;
            output = input;
            assert(Helper.setequal(...
                Helper.valid_neigh(input), ...
                output ...
            ));

            input = {};
            output = {};
            assert(Helper.setequal(...
                Helper.valid_neigh(input), ...
                output ...
            ));

            %% full_neigh
            input = [cellstr(strcat(Constants.neigh_tag, ["A_B_C", 'vvPOSNEGs', '_POS'])), {'not_neigh', ''}];
            output = {'A B C', "vv+-s", "not_neigh", "", " +"};
            assert(Helper.setequal(...
                Helper.full_neigh(input), ...
                output ...
            ));

            input = {};
            output = {};
            assert(Helper.setequal(...
                Helper.full_neigh(input), ...
                output ...
            ));

            %% Load data
            app = CytoMAP;
            data_loader.load_3d_dataset_csv(app);
            data_loader.load_2d_dataset_mat(app, 'c');

            %% get_neighs
            % Sample given - no neighborhoods
            assert(isempty(...
                Helper.get_neighs(app, 'Sample_1', 'MFIRSN') ...
            ));

            % Sample given - same neighborhoods
            assert(Helper.setequal(...
                Helper.get_neighs(app, 'all_2', 'MFIRSN'), ...
                {'Nei_rsn_every_poly_x_y', 'Nei_rsn_all_rect_x_y_'} ...
            ));

            % Multiple samples, all common neighs
            assert(Helper.setequal(...
                Helper.get_neighs(app, {'all_1', 'all_2'}), ...
                {'Nei_rsn_every_poly_x_y', 'Nei_rsn_all_rect_x_y_'} ...
            ));

            % Test ordering
            assert(Helper.setequal(...
                Helper.get_neighs(app, {'all_2', 'all_1'}), ...
                {'Nei_rsn_every_poly_x_y', 'Nei_rsn_all_rect_x_y_'} ...
            ));

            % Multiple samples, some common neighs
            assert(Helper.setequal(...
                Helper.get_neighs(app, {'other', 'all_2'}), ...
                {'Nei_rsn_every_poly_x_y'} ...
            ));

            % Multiple samples, one valid, multiple invalid
            assert(Helper.setequal(...
                Helper.get_neighs(app, {'NE1', 'all_2', 'NE2'}, 'MFIRSN'), ...
                {'Nei_rsn_every_poly_x_y', 'Nei_rsn_all_rect_x_y_'} ...
            ));

            % Multiple samples, multiple valid, multiple invalid
            assert(Helper.setequal(...
                Helper.get_neighs(app, {'NE1', 'other', 'NE2', 'no_lymph', 'ne3'}, 'MFIRSN'), ...
                {'Nei_rsn_every_poly_x_y'} ...
            ));

            Helper.func_exit(app);
        end

        function test_other_processing()
            %% valid_other
            % All weird cases
            input = {"_A", '______A_+_b_', "", "B", 'C'};
            assert(Helper.setequal(...
                Helper.valid_other(input), ...
                input ...
            ));

            input = Constants.other_names;
            output = strcat(Constants.other_tag, Constants.other_names);
            assert(Helper.setequal(...
                Helper.valid_other(input), ...
                output ...
            ));

            input = {};
            output = {};
            assert(Helper.setequal(...
                Helper.valid_other(input), ...
                output ...
            ));

            %% full_other
            input = [cellstr(strcat(Constants.other_tag, ["A_B_C", 'vvPOSNEGs', '_POS'])), {'not_other', ''}];
            output = {'A B C', "vv+-s", "not_other", "", " +"};
            assert(Helper.setequal(...
                Helper.full_other(input), ...
                output ...
            ));

            input = {};
            output = {};
            assert(Helper.setequal(...
                Helper.full_other(input), ...
                output ...
            ));

            %% Initialize CytoMAP, and load dataset
            app = CytoMAP;
            data_loader.load_3d_dataset_csv(app);
            data_loader.load_2d_dataset_mat(app, 'e');

            %% get_others
            % Single sample given
            assert(Helper.setequal(...
                Helper.get_others(app, 'all_2'), ...
                { ...
                    'OtherLocalDensityOf_All_Lymphocytes_CD11cPOS', ...
                    'OtherDistTo_All_Lymphocytes_CD11cPOS', ...
                    'OtherLocalDensityOf_All_T_Cells_OTNEGII_4GETPOS', ...
                    'OtherDistTo_All_T_Cells_OTNEGII_4GETPOS', ...
                    'OtherLocalDensityOf_All_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                    'OtherDistTo_All_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                    'OtherLocalDensityOf_All_T_Cells_OT_II_CXCR3POS', ...
                    'OtherDistTo_All_T_Cells_OT_II_CXCR3POS', ...
                    'OtherLocalDensityOf_All_T_cells_OTNEGII_T_Cells', ...
                    'OtherDistTo_All_T_cells_OTNEGII_T_Cells', ...
                    'OtherLocalDensityOf_All_rect_x_y', ...
                    'OtherDistTo_All_rect_x_y', ...
                    'OtherLocalDensityOf_All_poly_all_x_y', ...
                    'OtherDistTo_All_poly_all_x_y', ...
                    'OtherDistTopoly_1_x_y_Polygon', ...
                    'OtherDistToUDS_surf_every_x_y', ...
                    'Otherall_7' ...
                } ...
            ));

            % Default sample
            assert(Helper.setequal(...
                Helper.get_others(app), ...
                { ...
                    'OtherLocalDensityOf_All_Lymphocytes_CD11cPOS', ...
                    'OtherDistTo_All_Lymphocytes_CD11cPOS', ...
                    'OtherLocalDensityOf_All_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                    'OtherDistTo_All_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                    'OtherLocalDensityOf_All_T_Cells_OT_II_CXCR3POS', ...
                    'OtherDistTo_All_T_Cells_OT_II_CXCR3POS', ...
                    'OtherLocalDensityOf_All_T_cells_OTNEGII_T_Cells', ...
                    'OtherDistTo_All_T_cells_OTNEGII_T_Cells', ...
                    'OtherLocalDensityOf_All_rect_x_y', ...
                    'OtherDistTo_All_rect_x_y', ...
                    'OtherLocalDensityOf_All_poly_all_x_y', ...
                    'OtherDistTo_All_poly_all_x_y', ...
                    'OtherLocalDensityOf_All_rect_x_y_nos', ...
                    'OtherDistTo_All_rect_x_y_nos', ...
                    'OtherLocalDensityOf_All_poly_x_y_nos', ...
                    'OtherDistTo_All_poly_x_y_nos', ...
                    'OtherDistTopoly_2_x_y_Polygon', ...
                    'OtherDistToUDS_surf_every_x_y', ...
                    'Otherall_7', ...
                    'Othernos_17' ...
                } ...
            ));

            % Single sample, no others
            assert(isempty(Helper.get_others(app, 'Sample_1')));

            % Multiple samples, all common channels
            assert(Helper.setequal(...
                Helper.get_others(app, {'all_1', 'all_2'}), ...
                { ...
                    'OtherDistToUDS_surf_every_x_y', ...
                    'OtherDistTo_All_Lymphocytes_CD11cPOS', ...
                    'OtherDistTo_All_T_Cells_OTNEGII_4GETPOS', ...
                    'OtherDistTo_All_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                    'OtherDistTo_All_T_Cells_OT_II_CXCR3POS', ...
                    'OtherDistTo_All_T_cells_OTNEGII_T_Cells', ...
                    'OtherDistTo_All_poly_all_x_y', ...
                    'OtherDistTo_All_rect_x_y', ...
                    'OtherDistTopoly_1_x_y_Polygon', ...
                    'OtherLocalDensityOf_All_Lymphocytes_CD11cPOS', ...
                    'OtherLocalDensityOf_All_T_Cells_OTNEGII_4GETPOS', ...
                    'OtherLocalDensityOf_All_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                    'OtherLocalDensityOf_All_T_Cells_OT_II_CXCR3POS', ...
                    'OtherLocalDensityOf_All_T_cells_OTNEGII_T_Cells', ...
                    'OtherLocalDensityOf_All_poly_all_x_y', ...
                    'OtherLocalDensityOf_All_rect_x_y', ...
                    'Otherall_7' ...
                } ...
            ));

            % Multiple samples, some common channels
            assert(Helper.setequal(...
                Helper.get_others(app, {'all_2', 'no_T4'}), ...
                { ...
                    'OtherDistToUDS_surf_every_x_y', ...
                    'OtherDistTo_All_Lymphocytes_CD11cPOS', ...
                    'OtherDistTo_All_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                    'OtherDistTo_All_T_Cells_OT_II_CXCR3POS', ...
                    'OtherDistTo_All_T_cells_OTNEGII_T_Cells', ...
                    'OtherDistTo_All_poly_all_x_y', ...
                    'OtherDistTo_All_rect_x_y', ...
                    'OtherLocalDensityOf_All_Lymphocytes_CD11cPOS', ...
                    'OtherLocalDensityOf_All_T_Cells_OTNEGII_CXCR3POS_4getNEG', ...
                    'OtherLocalDensityOf_All_T_Cells_OT_II_CXCR3POS', ...
                    'OtherLocalDensityOf_All_T_cells_OTNEGII_T_Cells', ...
                    'OtherLocalDensityOf_All_poly_all_x_y', ...
                    'OtherLocalDensityOf_All_rect_x_y', ...
                    'Otherall_7' ...
                } ...
            ));

            % Multiple samples, no common channels
            assert(isempty(...
                Helper.get_others(app, {'other', 'Sample_1'}) ...
            ));

            Helper.func_exit(app);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Data Preparation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function test_DataPrep()
            %% Initialize CytoMAP, load Data
            app = CytoMAP;
            data_loader.load_2d_dataset_mat(app, 'c');

            %% Tests
            % 1st option
            option = 'Cellularity: Number of Cells / Neighborhood';
            smpl = 'all_2';
            gates = app.data.(smpl).MFIRSN.Properties.VariableNames(startsWith(app.data.(smpl).MFIRSN.Properties.VariableNames, Constants.gate_tag));
            gates = table2array(app.data.(smpl).MFIRSN(:, gates));
            out = Helper.Func_DataPrep(app.data.(smpl).MFIRSN, gates, option);
            assert(~any(out - gates, 'all'));

            % 2nd option
            option = 'Density: Number of Cells / Volume (Area if 2D) Per Neighborhood';
            smpl = 'all_2';
            gates = app.data.(smpl).MFIRSN.Properties.VariableNames(startsWith(app.data.(smpl).MFIRSN.Properties.VariableNames, Constants.gate_tag));
            gates = table2array(app.data.(smpl).MFIRSN(:, gates));
            out = Helper.Func_DataPrep(app.data.(smpl).MFIRSN, gates, option);
            assert(~any(out - (gates ./ app.data.(smpl).MFIRSN.Neigh_Area), 'all'));
            smpl = 'other';
            gates = app.data.(smpl).MFIRSN.Properties.VariableNames(startsWith(app.data.(smpl).MFIRSN.Properties.VariableNames, Constants.gate_tag));
            gates = table2array(app.data.(smpl).MFIRSN(:, gates));
            out = Helper.Func_DataPrep(app.data.(smpl).MFIRSN, gates, option);
            assert(~any(out - (gates ./ app.data.(smpl).MFIRSN.Neigh_Volume), 'all'));

            % 3rd option
            option = 'Global Composition: Number of Cells / Max Cells in Tissue Neighborhoods';
            smpl = 'all_2';
            gates = app.data.(smpl).MFIRSN.Properties.VariableNames(startsWith(app.data.(smpl).MFIRSN.Properties.VariableNames, Constants.gate_tag));
            gates = table2array(app.data.(smpl).MFIRSN(:, gates));
            out = Helper.Func_DataPrep(app.data.(smpl).MFIRSN, gates, option);
            assert(all(out <= 1, 'all'));
            assert(all(out == (gates ./ max(gates)), 'all'));

            % 4th option
            option = 'Binary: If cell is in neighborhood';
            smpl = 'all_2';
            gates = app.data.(smpl).MFIRSN.Properties.VariableNames(startsWith(app.data.(smpl).MFIRSN.Properties.VariableNames, Constants.gate_tag));
            gates = table2array(app.data.(smpl).MFIRSN(:, gates));
            out = Helper.Func_DataPrep(app.data.(smpl).MFIRSN, gates, option);
            assert(all((out == 1 | out == 0) & out == (gates > 1), 'all'));

            % 5th option
            option = 'Standardize: subtract Mean, divide by standard deviation';
            smpl = 'all_2';
            gates = app.data.(smpl).MFIRSN.Properties.VariableNames(startsWith(app.data.(smpl).MFIRSN.Properties.VariableNames, Constants.gate_tag));
            gates = table2array(app.data.(smpl).MFIRSN(:, gates));
            out = Helper.Func_DataPrep(app.data.(smpl).MFIRSN, gates, option);
            assert(all(abs(mean(out)) < 1e-10, 'all'));
            assert(all(abs(std(out) - 1) < 1e-10, 'all'));

            % 6th option
            option = 'Corrected Density: Number of Cells / Volume (Area if 2D) of Neighborhood Inside Tissue';
            smpl = 'all_2';
            gates = app.data.(smpl).MFIRSN.Properties.VariableNames(startsWith(app.data.(smpl).MFIRSN.Properties.VariableNames, Constants.gate_tag));
            gates = table2array(app.data.(smpl).MFIRSN(:, gates));
            out = Helper.Func_DataPrep(app.data.(smpl).MFIRSN, gates, option);
            Vol = app.data.(smpl).MFIRSN.Effective_Neigh_Area;
            Vol(Vol==0) = 1;
            assert(~any(out - (gates ./ Vol), 'all'));
            option = 'Corrected Density: Number of Cells / Volume (Area if 2D) of Neighborhood Inside Tissue';
            smpl = 'other';
            gates = app.data.(smpl).MFIRSN.Properties.VariableNames(startsWith(app.data.(smpl).MFIRSN.Properties.VariableNames, Constants.gate_tag));
            gates = table2array(app.data.(smpl).MFIRSN(:, gates));
            out = Helper.Func_DataPrep(app.data.(smpl).MFIRSN, gates, option);
            Vol = app.data.(smpl).MFIRSN.Effective_Neigh_Volume;
            Vol(Vol==0) = 1;
            assert(~any(out - (gates ./ Vol), 'all'));

            % 7th option
            option = 'Composition: Number of Cells / Total cells in Neighborhood';
            smpl = 'all_2';
            gates = app.data.(smpl).MFIRSN.Properties.VariableNames(startsWith(app.data.(smpl).MFIRSN.Properties.VariableNames, Constants.gate_tag));
            gates = table2array(app.data.(smpl).MFIRSN(:, gates));
            out = Helper.Func_DataPrep(app.data.(smpl).MFIRSN, gates, option);
            NCell = app.data.(smpl).MFIRSN.NCells;
            NCell(NCell==0) = 1;
            assert(~any(out - (gates ./ NCell), 'all'));

            % 8th option
            option = 'Composition+Density: Number of Cells / Total cells in Neighborhood; Density of Neighborhood';
            smpl = 'all_2';
            gates = app.data.(smpl).MFIRSN.Properties.VariableNames(startsWith(app.data.(smpl).MFIRSN.Properties.VariableNames, Constants.gate_tag));
            gates = table2array(app.data.(smpl).MFIRSN(:, gates));
            out = Helper.Func_DataPrep(app.data.(smpl).MFIRSN, gates, option);
            assert(size(out, 2) == size(gates, 2) + 1);
            NCell = app.data.(smpl).MFIRSN.NCells;
            NCell(NCell==0) = 1;
            gates = gates ./ NCell;
            NCell = NCell./max(NCell);
            assert(all(out(:, 1:end-1) - gates == 0, 'all'));
            assert(all(out(:, end) - (0.5 * NCell) == 0, 'all'));

            %% Clean Workspace, close CytoMAP
            Helper.func_exit(app);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Networks, and any_sample checks.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TODO: Remake Helper.add_model to not use user interface (i.e. add backend function)
        % TAG: Version 2, it seems extremely hard, since it's interleaving with the returns, and processing a lot.
        % function test_add_model()
        % end

        function test_any_sample()
            % Create CytoMAP instance
            app = CytoMAP;

            % No Samples
            assert(~Helper.any_sample(app));

            % Single Sample
            data_loader.load_3d_dataset_mat(app);
            assert(Helper.any_sample(app));

            % Multiple Samples
            data_loader.load_2d_dataset_csv(app);
            assert(Helper.any_sample(app));
            Helper.func_exit(app);
        end

        function test_any_net()
            % Create CytoMAP instance
            app = CytoMAP;

            % No Samples/No Nets
            assert(~Helper.any_net(app));

            % Sample/No Nets Sample
            data_loader.load_2d_dataset_mat(app, 'a');
            assert(~Helper.any_net(app));

            % Samples/Nets
            data_loader.load_2d_dataset_mat(app, 'e');
            assert(Helper.any_net(app));
            Helper.func_exit(app);
        end

        function test_assert_valid_net()
            app = CytoMAP;

            % Not a struct test
            assert(~Helper.assert_valid_net(table));

            % Empty struct/Invalid struct
            assert(~Helper.assert_valid_net(struct));

            % Valid struct
            data_loader.load_2d_dataset_mat(app, 'e');
            assert(Helper.assert_valid_net(app.net.all_7));

            Helper.func_exit(app);
        end

        function test_get_net()
            %% Create CytoMAP instance
            app = CytoMAP;

            % No nets at all
            assert(isempty(...
                Helper.get_net(app) ...
            ));
            assert(isempty(...
                Helper.get_net(app, '') ...
            ));

            %% Load Data
            data_loader.load_2d_dataset_mat(app, 'e');

            %% Tests
            %% Purely correct net
            % Single as char
            output = Helper.get_net(app, 'all_7');
            assert(isequal(output, app.net.all_7));

            % Single as cell
            output = Helper.get_net(app, {'all_7'});
            assert(iscell(output) && isequal(output{1}, app.net.all_7));

            % Multiple in cells
            output = Helper.get_net(app, {'all_7', 'other_5', 'nos_17'});
            assert(...
                iscell(output) && ...
                numel(output) == 3 && ...
                isequal(output{1}, app.net.all_7) && ...
                isequal(output{2}, app.net.other_5) && ...
                isequal(output{3}, app.net.nos_17) ...
            );

            % With other notation
            output = Helper.get_net(app, Helper.valid_other({'Regionall_7', 'Modelnos_17'}));
            assert(...
                iscell(output) && ...
                numel(output) == 2 && ...
                isequal(output{1}, app.net.all_7) && ...
                isequal(output{2}, app.net.nos_17) ...
            );

            % Not valid variable
            output = Helper.get_net(app, Helper.valid_other({'Regionall____ 7', 'Modelnos _ 17'}));
            assert(...
                iscell(output) && ...
                numel(output) == 2 && ...
                isequal(output{1}, app.net.all_7) && ...
                isequal(output{2}, app.net.nos_17) ...
            );

            % Not valid variable + non-existent networks
            output = Helper.get_net(app, Helper.valid_other({'NE1', 'Regionall____ 7', 'NE2', 'Modelnos _ 17', 'NE3'}));
            assert(...
                iscell(output) && ...
                numel(output) == 5 && ...
                isempty(output{1}) && ...
                isequal(output{2}, app.net.all_7) && ...
                isempty(output{3}) && ...
                isequal(output{4}, app.net.nos_17) && ...
                isempty(output{5}) ...
            );
            Helper.func_exit(app);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % User Tables management (choosing samples/phenotype outside of plotting)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function test_reorder_cols()
            % Init CytoMAP, Load data
            app = CytoMAP;
            data_loader.load_2d_dataset_mat(app);

            % Test Empty Table
            assert(isempty(Helper.reorder_cols(table)));

            % Test Table with only X, Y cols
            out = Helper.reorder_cols(table('Size', [5, 2], 'VariableTypes', {'double', 'double'}, 'VariableNames', {'Y', 'X'}));
            assert(strcmp(out.Properties.VariableNames{1}, "X"));
            assert(strcmp(out.Properties.VariableNames{2}, "Y"));
            assert(strcmp(out.Properties.VariableNames{3}, "Z"));

            % Test Table with only X, Y, Z cols
            out = Helper.reorder_cols(table('Size', [5, 3], 'VariableTypes', {'double', 'double', 'double'}, 'VariableNames', {'Y', 'Z', 'X'}));
            assert(strcmp(out.Properties.VariableNames{1}, "X"));
            assert(strcmp(out.Properties.VariableNames{2}, "Y"));
            assert(strcmp(out.Properties.VariableNames{3}, "Z"));

            % Test real life data
            out = Helper.reorder_cols(app.data.all_1.AllCells);
            assert(strcmp(out.Properties.VariableNames{1}, "X"));
            assert(strcmp(out.Properties.VariableNames{2}, "Y"));
            assert(strcmp(out.Properties.VariableNames{3}, "Z"));

            % Close window etc.
            Helper.func_exit(app);
        end

        function test_populate_table()
            %% Simple Tests
            % Start app, and make it into the table
            app = CytoMAP;

            % Make sure it works properly for the empty data
            tab = Helper.populate_table(app);
            assert(isempty(tab));

            data_loader.load_3d_dataset_csv(app);

            % Simple, default settings
            tab = Helper.populate_table(app);
            assert(all(cell2mat(tab(:, 1))));
            assert(tab{strcmp(tab(:, 5), app.DataN.Value), 4});
            assert(~any(cell2mat(tab(~strcmp(tab(:, 5), app.DataN.Value), 4))));
            assert(size(tab, 2) == 5);

            % Not a proper MFI given
            tab2 = Helper.populate_table(app, 'mfi', 'MFICCN');
            assert(isempty(tab2));

            Helper.func_exit(app);

            %% Medium Dataset/Multiple samples chosen Tests
            % Start app, and load data
            app = CytoMAP;
            data_loader.load_3d_dataset_wsp(app);
            data_loader.load_2d_dataset_mat(app, 'e');
            % MFI - Single Dataset
            tab = Helper.populate_table(app, 'mfi', 'MFIRSN', 'smpls', {'no_lymph'});
            % Everything but the gates should be in
            assert(size(tab, 1) == size(app.data.no_lymph.MFIRSN, 2) - size(app.data.no_lymph.GateTags, 2));
            assert(~isempty(tab));
            assert(size(tab, 2) == 8);

            % MFI - Multiple Datasets, with only one network shared
            % (first one has 2 networks trained on it)
            tab2 = Helper.populate_table(app, 'mfi', 'MFIRSN', 'smpls', {'no_lymph', 'all_1'});
            assert(size(tab2, 1) == size(tab, 1) - 7);  % Some distances/models should dissapear.
            assert(~isempty(tab2));
            assert(size(tab, 2) == 8);

            % Multiple samples, default chosen
            tab = Helper.populate_table(app);
            idx = strcmp(tab(:, 5), app.DataN.Value);
            assert(tab{idx, 4});
            assert(~any(cell2mat(tab(~idx, 4))));
            assert(size(tab, 2) == 5);

            % All samples chosen at the same time; no common gates
            tab = Helper.populate_table(app, 'smpls', fieldnames(app.data));
            assert(all(cell2mat(tab(:, 4))));
            Helper.func_exit(app);

            %% Test prev table, and other key-word args; CytoMAP initData Loading
            app = CytoMAP;
            data_loader.load_2d_dataset_mat(app, 'e');

            %% Test prev table: Without MFIRSN
            tab_prev = Helper.populate_table(app, 'smpls', {'no_lymph', 'all_2'});

            % Do some changes so that we can see it prev_table works
            phns_to_mod = {'All/20170110 Day4-5 Alum2 T Cells OT-II 4GET+', 'All/20170110 Day4-5 Alum2 T cells OT-II T Cells'};
            idx_ph_ne = ~cellfun('isempty', tab_prev(:, 3));
            idx_ph = idx_ph_ne;
            idx_ph(idx_ph == 1) = ismember(tab_prev(idx_ph, 3), phns_to_mod);
            tab_prev(idx_ph, 2) = {false};
            tab_prev(idx_ph, 1) = {2};

            % Same samples check
            tab_new = Helper.populate_table(app, 'smpls', {'no_lymph', 'all_2'}, 'prev_table', tab_prev);
            assert(Helper.setequal(tab_new(idx_ph_ne, 1), tab_prev(idx_ph_ne, 1)));
            assert(Helper.setequal(tab_new(idx_ph_ne, 2), tab_prev(idx_ph_ne, 2)));
            assert(Helper.setequal(tab_new(idx_ph_ne, 3), tab_prev(idx_ph_ne, 3)));

            % Different samples check
            tab_new = Helper.populate_table(app, 'smpls', 'no_lymph', 'prev_table', tab_prev);
            idx_ph_ne = ~cellfun('isempty', tab_new(:, 3));
            idx_ph = idx_ph_ne;
            idx_ph(idx_ph == 1) = ismember(tab_new(idx_ph, 3), phns_to_mod);

            % Test Phenotype choices. In both columns. Both the ones changed, and the ones unchanged.
            assert(all(~cell2mat(tab_new(idx_ph, 2))));
            assert(all(cell2mat(tab_new(idx_ph, 1)) == 2));
            assert(all(cell2mat(tab_new(idx_ph_ne & ~idx_ph, 2))));
            assert(all(cell2mat(tab_new(idx_ph_ne & ~idx_ph, 1)) == 1));

            %% Test prev table: With MFIRSN
            tab_prev = Helper.populate_table(app, 'mfi', 'MFIRSN', 'smpls', {'no_lymph', 'all_2'});

            % Do some changes so that we can see it prev_table works
            phns_to_mod = {'All/20170110 Day4-5 Alum2 T Cells OT-II 4GET+', 'All/20170110 Day4-5 Alum2 T cells OT-II T Cells'};
            idx_ph_ne = ~cellfun('isempty', tab_prev(:, 3));
            idx_ph = idx_ph_ne;
            idx_ph(idx_ph == 1) = ismember(tab_prev(idx_ph, 3), phns_to_mod);
            tab_prev(idx_ph, 2) = {false};
            tab_prev(idx_ph, 1) = {2};

            mf_to_mod = {'Sphericity', 'Volume', 'NCells', 'Y', 'Z', '1900Para 11'};
            idx_mf_ne = ~cellfun('isempty', tab_prev(:, 6));
            idx_mf = idx_mf_ne;
            idx_mf(idx_mf == 1) = ismember(tab_prev(idx_mf, 6), mf_to_mod);
            tab_prev(idx_mf, 5) = {false};
            tab_prev(idx_mf, 4) = {5};

            % Same samples check
            tab_new = Helper.populate_table(app, 'mfi', 'MFIRSN', 'smpls', {'no_lymph', 'all_2'}, 'prev_table', tab_prev);
            assert(Helper.setequal(tab_new(idx_ph_ne, 1), tab_prev(idx_ph_ne, 1)));
            assert(Helper.setequal(tab_new(idx_ph_ne, 2), tab_prev(idx_ph_ne, 2)));
            assert(Helper.setequal(tab_new(idx_ph_ne, 3), tab_prev(idx_ph_ne, 3)));
            assert(Helper.setequal(tab_new(idx_mf_ne, 4), tab_prev(idx_mf_ne, 4)));
            assert(Helper.setequal(tab_new(idx_mf_ne, 5), tab_prev(idx_mf_ne, 5)));
            assert(Helper.setequal(tab_new(idx_mf_ne, 6), tab_prev(idx_mf_ne, 6)));

            % Different samples check
            tab_new = Helper.populate_table(app, 'mfi', 'MFIRSN', 'smpls', 'no_lymph', 'prev_table', tab_prev);
            idx_ph_ne = ~cellfun('isempty', tab_new(:, 3));
            idx_ph = idx_ph_ne;
            idx_ph(idx_ph == 1) = ismember(tab_new(idx_ph, 3), phns_to_mod);
            idx_mf_ne = ~cellfun('isempty', tab_new(:, 6));
            idx_mf = idx_mf_ne;
            idx_mf(idx_mf == 1) = ismember(tab_new(idx_mf, 6), mf_to_mod);

            % Test Phenotype choices. In both columns. Both the ones changed, and the ones unchanged.
            assert(all(~cell2mat(tab_new(idx_ph, 2))));
            assert(all(cell2mat(tab_new(idx_ph, 1)) == 2));
            assert(all(cell2mat(tab_new(idx_ph_ne & ~idx_ph, 2))));
            assert(all(cell2mat(tab_new(idx_ph_ne & ~idx_ph, 1)) == 1));
            % Test MFI choices. In both columns. Both the ones changed, and the ones unchanged.
            assert(all(~cell2mat(tab_new(idx_mf, 5))));
            assert(all(cell2mat(tab_new(idx_mf, 4)) == 5));
            assert(all(cell2mat(tab_new(idx_mf_ne & ~idx_mf, 5))));
            assert(all(cell2mat(tab_new(idx_mf_ne & ~idx_mf, 4)) == 1));

            %% Test fill args: With/WithoutMFI
            tab = Helper.populate_table(app, 'smpls', 'no_lymph', 'fill_weight', 5, 'fill_checkbox', false);
            idx_ph = ~cellfun('isempty', tab(:, 3));
            assert(all(~cell2mat(tab(idx_ph, 2))));
            assert(all(cell2mat(tab(idx_ph, 1)) == 5));

            tab = Helper.populate_table(app, 'mfi', 'MFIRSN', 'smpls', 'no_lymph', 'fill_weight', 5, 'fill_checkbox', false);
            idx_ph = ~cellfun('isempty', tab(:, 3));
            idx_mf = ~cellfun('isempty', tab(:, 6));
            assert(all(~cell2mat(tab(idx_ph, 2))));
            assert(all(cell2mat(tab(idx_ph, 1)) == 5));
            assert(all(~cell2mat(tab(idx_mf, 5))));
            assert(all(cell2mat(tab(idx_mf, 4)) == 5));

            %% Test fill args combined with prev table: MFI only
            tab_prev = Helper.populate_table(app, 'mfi', 'MFIRSN', 'smpls', {'no_lymph', 'no_T4'});

            % Do some changes so that we can see it prev_table works
            phns_to_mod = {'All/20170110 Day4-5 Alum2 T Cells OT II CXCR3+', 'All/20170110 Day4-5 Alum2 T cells OT-II T Cells'};
            idx_ph_ne = ~cellfun('isempty', tab_prev(:, 3));
            idx_ph = idx_ph_ne;
            idx_ph(idx_ph == 1) = ismember(tab_prev(idx_ph, 3), phns_to_mod);
            tab_prev(idx_ph, 2) = {false};
            tab_prev(idx_ph, 1) = {2};

            mf_to_mod = {'Sphericity', 'Volume', 'NCells', 'Y', 'Z', '1900Para 11'};
            idx_mf_ne = ~cellfun('isempty', tab_prev(:, 6));
            idx_mf = idx_mf_ne;
            idx_mf(idx_mf == 1) = ismember(tab_prev(idx_mf, 6), mf_to_mod);
            tab_prev(idx_mf, 5) = {false};
            tab_prev(idx_mf, 4) = {5};

            % Different samples check
            tab_new = Helper.populate_table(app, 'mfi', 'MFIRSN', 'smpls', 'no_lymph', 'prev_table', tab_prev, 'fill_weight', 25, 'fill_checkbox', false);
            idx_ph_ne = ~cellfun('isempty', tab_new(:, 3));
            phns_new = setdiff(tab_new(idx_ph_ne, 3), tab_prev(~cellfun('isempty', tab_prev(:, 3)), 3));
            idx_ph_keep = idx_ph_ne;
            idx_ph_keep(idx_ph_keep == 1) = ismember(tab_new(idx_ph_keep, 3), phns_to_mod);
            idx_ph_new = idx_ph_ne;
            idx_ph_new(idx_ph_new == 1) = ismember(tab_new(idx_ph_new, 3), phns_new);
            idx_mf_ne = ~cellfun('isempty', tab_new(:, 6));
            mf_new = setdiff(tab_new(idx_mf_ne, 6), tab_prev(~cellfun('isempty', tab_prev(:, 6)), 6));
            idx_mf_keep = idx_mf_ne;
            idx_mf_keep(idx_mf_keep == 1) = ismember(tab_new(idx_mf_keep, 6), mf_to_mod);
            idx_mf_new = idx_mf_ne;
            idx_mf_new(idx_mf_new == 1) = ismember(tab_new(idx_mf_new, 6), mf_new);

            % Test Phenotype choices. In both columns. Both the ones changed, new and unchanged.
            assert(all(~cell2mat(tab_new(idx_ph_keep, 2))));
            assert(all(cell2mat(tab_new(idx_ph_keep, 1)) == 2));
            assert(all(~cell2mat(tab_new(idx_ph_new, 2))));
            assert(all(cell2mat(tab_new(idx_ph_new, 1)) == 25));
            assert(all(cell2mat(tab_new(idx_ph_ne & ~idx_ph_new & ~idx_ph_keep, 2))));
            assert(all(cell2mat(tab_new(idx_ph_ne & ~idx_ph_new & ~idx_ph_keep, 1)) == 1));
            % Test MFI choices. In both columns. Both the ones changed, new and unchanged.
            assert(all(~cell2mat(tab_new(idx_mf_keep, 5))));
            assert(all(cell2mat(tab_new(idx_mf_keep, 4)) == 5));
            assert(all(~cell2mat(tab_new(idx_mf_new, 5))));
            assert(all(cell2mat(tab_new(idx_mf_new, 4)) == 25));
            assert(all(cell2mat(tab_new(idx_mf_ne & ~idx_mf_keep & ~idx_mf_new, 5))));
            assert(all(cell2mat(tab_new(idx_mf_ne & ~idx_mf_keep & ~idx_mf_new, 4)) == 1));

            Helper.func_exit(app);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Merging samples
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function test_merge_AllCells()
            %% Make CytoMAP instance, load datasets
            app = CytoMAP;
            data_loader.load_2d_dataset_mat(app);

            %% Single Sample
            [ch_oth_out, g_tab_out, g_tag_out] = Helper.merge_AllCells(app, {'all_2'});
            assert(...
                size(ch_oth_out, 1) == size(g_tab_out, 1) && ...
                size(ch_oth_out, 1) == size(app.data.all_2.AllCells, 1) && ...
                size(ch_oth_out, 2) + size(g_tab_out, 2) == size(app.data.all_2.AllCells, 2) && ...
                Helper.setequal(...
                    [g_tab_out.Properties.VariableNames, ch_oth_out.Properties.VariableNames], ...
                    app.data.all_2.AllCells.Properties.VariableNames ...
                ) && ...
                size(g_tab_out, 2) == size(g_tag_out, 2) && ...
                size(g_tag_out, 1) == 2 ...
            );

            %% Two same samples
            [ch_oth_out, g_tab_out, g_tag_out] = Helper.merge_AllCells(app, {'all_1', 'all_2'});
            assert(...
                size(ch_oth_out, 1) == size(g_tab_out, 1) && ...
                size(ch_oth_out, 1) == 2 * size(app.data.all_2.AllCells, 1) && ...
                size(unique(ch_oth_out, 'rows'), 1) == size(app.data.all_2.AllCells, 1) && ...
                size(ch_oth_out, 2) + size(g_tab_out, 2) == size(app.data.all_2.AllCells, 2) && ...
                Helper.setequal(...
                    [g_tab_out.Properties.VariableNames, ch_oth_out.Properties.VariableNames], ...
                    intersect(app.data.all_2.AllCells.Properties.VariableNames, app.data.all_1.AllCells.Properties.VariableNames) ...
                ) && ...
                size(g_tab_out, 2) == size(g_tag_out, 2) && ...
                size(g_tag_out, 1) == 2 ...
            );

            %% Two samples, no overlap (Apart from required/most typical channels, such as X, Y, Z etc.)
            [ch_oth_out, g_tab_out, g_tag_out] = Helper.merge_AllCells(app, {'all_2', 'other'});
            common_chcols = intersect(app.data.other.AllCells.Properties.VariableNames, app.data.all_2.AllCells.Properties.VariableNames);
            common_chcols = common_chcols(~startsWith(common_chcols, Constants.gate_tag));
            assert(...
                size(ch_oth_out, 1) == size(g_tab_out, 1) && ...
                size(ch_oth_out, 1) == size(app.data.other.AllCells, 1) + size(app.data.all_2.AllCells, 1) && ...
                size(unique(ch_oth_out(:, Helper.get_channels(app, {'all_2', 'other'}))), 1) == size(ch_oth_out, 1) && ...
                Helper.setequal(...
                    ch_oth_out.Properties.VariableNames, ...
                    common_chcols ...
                ) && ...
                size(ch_oth_out, 2) == numel(common_chcols) ...
            );

            %% Multi-samples case
            [ch_oth_out, g_tab_out, g_tag_out] = Helper.merge_AllCells(app, fieldnames(app.data));
            total_size = 0;
            for smpl_idx=1:numel(app.DataN.Items)
                total_size = total_size + size(app.data.(app.DataN.Items{smpl_idx}).AllCells, 1);
            end
            assert(...
                size(ch_oth_out, 1) == size(g_tab_out, 1) && ...
                size(ch_oth_out, 1) == total_size && ...
                size(unique(ch_oth_out(:, Helper.get_channels(app, fieldnames(app.data)))), 1) == ...
                    size(app.data.other.AllCells, 1) + size(app.data.all_2.AllCells, 1) && ...
                Helper.setequal(...
                    ch_oth_out.Properties.VariableNames, ...
                    common_chcols ...
                ) && ...
                size(ch_oth_out, 2) == numel(common_chcols) ...
            );

            %% Close
            Helper.func_exit(app);
        end

        function test_merge_MFI()
            %% Init CytoMAP, load data
            app = CytoMAP;
            data_loader.load_2d_dataset_mat(app, 'd');

            %% Tests
            % 2 samples - Exactly the same
            mfi_tab = Helper.merge_MFI(app, {'all_1', 'all_2'}, 'MFIRSN');
            assert(size(mfi_tab, 1) == size(app.data.all_2.MFIRSN, 1));
            assert(all(size(mfi_tab) - size(app.data.all_1.MFIRSN) == 0));
            assert(Helper.setequal( ...
                app.data.all_2.MFIRSN.Properties.VariableNames, ...
                mfi_tab.Properties.VariableNames ...
            ));
            assert(Helper.setequal( ...
                app.data.all_1.MFIRSN.Properties.VariableNames, ...
                mfi_tab.Properties.VariableNames ...
            ));

            % 2 samples - Overlapping
            mfi_tab = Helper.merge_MFI(app, {'no_lymph', 'no_T4'}, 'MFIRSN');
            assert(size(mfi_tab, 1) == size(app.data.no_lymph.MFIRSN, 1) + size(app.data.no_T4.MFIRSN, 1));
            assert(Helper.setequal( ...
                setdiff(app.data.no_lymph.MFIRSN.Properties.VariableNames, mfi_tab.Properties.VariableNames), ...
                {'Gate_1', 'Gate_2', 'Gate_3', 'Gate_4'} ...
            ));
            assert(Helper.setequal( ...
                setdiff(app.data.no_T4.MFIRSN.Properties.VariableNames, mfi_tab.Properties.VariableNames), ...
                {'Gate_1', 'Gate_2', 'Gate_3', 'Gate_4'} ...
            ));

            % 2 samples - No overlap
            mfi_tab = Helper.merge_MFI(app, {'all_2', 'other'}, 'MFIRSN');
            assert(size(mfi_tab, 1) == size(app.data.all_2.MFIRSN, 1) + size(app.data.other.MFIRSN, 1));
            correct = intersect(app.data.all_2.MFIRSN.Properties.VariableNames, app.data.other.MFIRSN.Properties.VariableNames);
            correct = correct(~startsWith(correct, Constants.gate_tag));
            assert(Helper.setequal( ...
                mfi_tab.Properties.VariableNames, ...
                correct ...
            ));

            % Multi-sample case
            smpls = {'all_2', 'all_1', 'no_lymph', 'no_T4'};
            mfi_tab = Helper.merge_MFI(app, smpls, 'MFIRSN');
            size_sum = 0;
            common_cols = app.data.(smpls{1}).MFIRSN.Properties.VariableNames;
            for s_idx=1:numel(smpls)
                common_cols = intersect(common_cols, app.data.(smpls{s_idx}).MFIRSN.Properties.VariableNames);
                size_sum = size_sum + size(app.data.(smpls{s_idx}).MFIRSN, 1);
                assert(size(mfi_tab, 1) >= size(app.data.(smpls{s_idx}).MFIRSN, 1));
            end
            common_cols = common_cols(~startsWith(common_cols, Constants.gate_tag));
            assert(size(mfi_tab, 1) <= size_sum);
            assert(Helper.setequal( ...
                setdiff(app.data.all_2.MFIRSN.Properties.VariableNames, mfi_tab.Properties.VariableNames), ...
                {'Gate_1', 'Gate_2', 'Gate_3', 'Gate_4', 'Gate_5'} ...
            ));
            assert(Helper.setequal( ...
                setdiff(app.data.all_1.MFIRSN.Properties.VariableNames, mfi_tab.Properties.VariableNames), ...
                {'Gate_1', 'Gate_2', 'Gate_3', 'Gate_4', 'Gate_5'} ...
            ));
            assert(Helper.setequal( ...
                setdiff(app.data.no_lymph.MFIRSN.Properties.VariableNames, mfi_tab.Properties.VariableNames), ...
                {'Gate_1', 'Gate_2', 'Gate_3', 'Gate_4'} ...
            ));
            assert(Helper.setequal( ...
                setdiff(app.data.no_T4.MFIRSN.Properties.VariableNames, mfi_tab.Properties.VariableNames), ...
                {'Gate_1', 'Gate_2', 'Gate_3', 'Gate_4'} ...
            ));
            assert(Helper.setequal( ...
                common_cols, ...
                mfi_tab.Properties.VariableNames ...
            ));

            %% Exit Test
            Helper.func_exit(app);
        end

        function test_intersect_smpls()
            %% Initialize CytoMAP instance, load data
            app = CytoMAP;
            data_loader.load_2d_dataset_mat(app);

            %% Tests
            % 2 samples - identical AllCells
            new_smpl = Helper.intersect_smpls(app, {'all_1', 'all_2'});
            assert(size(new_smpl.AllCells, 2) == size(app.data.all_1.AllCells, 2));
            assert(size(new_smpl.AllCells, 1) >= size(app.data.all_1.AllCells, 1));
            assert(Helper.setequal( ...
                new_smpl.AllCells.Properties.VariableNames, ...
                app.data.all_1.AllCells.Properties.VariableNames ...
            ));
            assert(Helper.setequal( ...
                new_smpl.AllCells.Properties.VariableNames(startsWith(new_smpl.AllCells.Properties.VariableNames, Constants.gate_tag)), ...
                new_smpl.GateTags.Properties.VariableNames ...
            ));
            assert(Helper.setequal( ...
                new_smpl.GateTags{2, :}, ...
                app.data.all_1.GateTags{2, :} ...
            ));
            assert(Helper.setequal( ...
                new_smpl.GateTags{1, :}, ...
                app.data.all_1.GateTags{1, :} ...
            ));
            assert(new_smpl.MetaData.Ver == Constants.CURR_VER);

            %% Clean workspace, Close CytoMAP
            Helper.func_exit(app);
        end
    end
end
