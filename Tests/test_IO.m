classdef test_IO
    methods (Static)
        function Import_Definitions_Func
            addpath(['..' filesep 'CytoMAP' filesep]);
            addpath(['..' filesep 'Tests' filesep]);
            addpath(['..' filesep 'CytoMAP' filesep '3rdPartyFunctions']);

            import Helper.*;
            import CytoMAP.*;
            import data_loader.*;
        end

        function test_load_wsp()
            app = CytoMAP;
            disp('  Initialized CytoMAP instance');
            assert(isempty(app.data));
            disp('  asserted that app.data isempty (CytoMAP Initialized correctly)');
            data_loader.load_3d_dataset_wsp(app);
            disp('  3D wsp dataset was loaded');
            assert(Helper.setequal(...
                {
                    'x20190113_auLN_Settings2_LN1_All', ...
                    'x20190113_auLN_Settings2_LN5_All', ...
                    'x20190113_auLN_Settings2_LN3_All', ...
                    'x20190113_auLN_Settings2_LN4_All', ...
                    'x20190113_auLN_Settings2_LN2_All', ...
                    'x20190113_bLN_Settings2_LN6_All', ...
                    'x20190113_auLN_Settings2_LN7_All', ...
                    'x20190113_bLN_Settings2_LN9_All', ...
                    'x20190113_bLN_Settings2_LN8_All' ...
                }, ...
                fieldnames(app.data) ...
            ));
            disp('  Checked that dataset was loaded correctly')
            Helper.func_exit(app);
        end

        function test_load_mat()
            app = CytoMAP;
            disp('  Initialized CytoMAP instance')
            assert(isempty(app.data));
            disp('  asserted that app.data isempty (CytoMAP Initialized correctly)')
            data_loader.load_2d_dataset_mat(app);
            disp('  2D mat dataset was loaded')
            assert(Helper.setequal(...
                {
                    'all_1', ...
                    'all_2', ...
                    'no_lymph', ...
                    'no_T4', ...
                    'other'
                }, ...
                fieldnames(app.data) ...
                ) ...
            );
            disp('  Checked that dataset was loaded correctly')
            Helper.func_exit(app);

            app = CytoMAP;
            disp('  Initialized CytoMAP instance')
            assert(isempty(app.data));
            disp('  asserted that app.data isempty (CytoMAP Initialized correctly)')
            data_loader.load_3d_dataset_mat(app);
            disp('  3D mat dataset was loaded')
            assert(Helper.setequal(...
                {
                    'All_1', ...
                    'All_2', ...
                    'Only_BVs', ...
                    'Only_Myleoids', ...
                    'Main_BV_Myleoid' ...
                }, ...
                fieldnames(app.data) ...
                ) ...
            );
            disp('  Checked that dataset was loaded correctly')
            Helper.func_exit(app);
        end

        function test_load_csv()
            app = CytoMAP;
            disp('  Initialized CytoMAP instance')
            assert(isempty(app.data));
            disp('  asserted that app.data isempty (CytoMAP Initialized correctly)')
            data_loader.load_2d_dataset_csv(app);
            disp('  2D csv dataset was loaded')
            assert(Helper.setequal(...
                    fieldnames(app.data), ...
                    {'Sample_1'}...
                )...
            );
            disp('  Checked that dataset was loaded correctly')
            Helper.func_exit(app);
        end

        function test_combined()
            app = CytoMAP;
            assert(isempty(app.data));

            data_loader.load_2d_dataset_csv(app);
            assert(Helper.setequal(...
                {
                    'Sample_1' ...
                }, ...
                fieldnames(app.data) ...
                ) ...
            );

            data_loader.load_3d_dataset_mat(app);
            assert(Helper.setequal(...
                {
                    'All_1', ...
                    'All_2', ...
                    'Only_BVs', ...
                    'Only_Myleoids', ...
                    'Main_BV_Myleoid', ...
                    'Sample_1' ...
                }, ...
                fieldnames(app.data) ...
                ) ...
            );

            % Check for adding new things to the dataset from mat
            data_loader.load_2d_dataset_mat(app);
            assert(Helper.setequal( ...
                {
                    'all_1', ...
                    'all_2', ...
                    'no_lymph', ...
                    'no_T4', ...
                    'other', ...
                    'All_1', ...
                    'All_2', ...
                    'Only_BVs', ...
                    'Only_Myleoids', ...
                    'Main_BV_Myleoid', ...
                    'Sample_1' ...
                }, ...
                fieldnames(app.data) ...
                ) ...
            );

            % Check for adding new things to the dataset from wsp
            data_loader.load_3d_dataset_wsp(app);
            assert(Helper.setequal(...
                {
                    'All_1', ...
                    'All_2', ...
                    'Only_BVs', ...
                    'Only_Myleoids', ...
                    'Main_BV_Myleoid', ...
                    'all_1', ...
                    'all_2', ...
                    'no_lymph', ...
                    'no_T4', ...
                    'other', ...
                    'Sample_1', ...
                    'x20190113_auLN_Settings2_LN1_All', ...
                    'x20190113_auLN_Settings2_LN5_All', ...
                    'x20190113_auLN_Settings2_LN3_All', ...
                    'x20190113_auLN_Settings2_LN4_All', ...
                    'x20190113_auLN_Settings2_LN2_All', ...
                    'x20190113_bLN_Settings2_LN6_All', ...
                    'x20190113_auLN_Settings2_LN7_All', ...
                    'x20190113_bLN_Settings2_LN9_All', ...
                    'x20190113_bLN_Settings2_LN8_All' ...
                }, ...
                fieldnames(app.data) ...
                ) ...
            );
            Helper.func_exit(app);
        end
    end
end
