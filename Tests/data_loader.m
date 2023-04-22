classdef data_loader
    %#ok<*UNRCH>

    methods (Static)
        function load_2d_dataset_mat(app, mode)
            if ~exist('mode', 'var')
                mode = "a";
            end
            if ~isstring(mode)
                mode = string(mode);
            end

            current_path = fileparts(mfilename('fullpath'));
            path_to_data = ['test_data' filesep '2D' filesep ];
            cd([current_path filesep path_to_data]);
            IO.func_load_backend(app, ...
                {{char("2D_dataset_" + mode + ".mat")}}, ...
                {path_to_data}, ...
                current_path ...
            );
            cd(current_path);
        end

        function load_2d_dataset_csv(app)
            current_path = fileparts(mfilename('fullpath'));
            path_to_data = ['test_data' filesep '2D' filesep 'csv' filesep];
            files = {{ ...
                '20170110 Lymphocytes_CD11c+.csv', ...
                '20170110 T Cells OT-II_4GET+.csv', ...
                '20170110 T Cells OT_II_CXCR3+.csv', ...
                '20170110 T Cells OT-II_CXCR3+_4get-.csv', ...
                '20170110 T cells_OT-II_T_Cells.csv'
            }};
            cd([current_path filesep path_to_data]);
            SampleNames = IO.get_sample_names(app, {path_to_data}, files);
            dat = IO.get_dat({path_to_data}, files, SampleNames);
            tab_data = IO.get_table_data(files, dat, SampleNames);
            IO.func_load_execute(...
                [], ...
                app, dat, files, {path_to_data}, SampleNames{end}, tab_data ...
            );
            cd(current_path);
        end

        function load_3d_dataset_mat(app, mode)
            if ~exist('mode', 'var')
                mode = "a";
            end
            if ~isstring(mode)
                mode = string(mode);
            end

            current_path = fileparts(mfilename('fullpath'));
            path_to_data = ['test_data' filesep '3D' filesep];
            cd([current_path filesep path_to_data]);
            IO.func_load_backend(app, ...
                {{char("3D_multi_thick_" + mode + ".mat")}}, ...
                {path_to_data}, ...
                current_path ...
            );
            cd(current_path)
        end
        
        function load_3d_dataset_csv(app)
            
            current_path = fileparts(mfilename('fullpath'));
            path_to_data = ['test_data' filesep '3D' filesep 'csv' filesep];
            files = {{ ...
                '0.csv', ...
                '1.csv', ...
                '2.csv', ...
                '3.csv', ...
                '4.csv' ...
            }};
            cd([current_path filesep path_to_data]);
            SampleNames = IO.get_sample_names(app, {path_to_data},files);
            dat = IO.get_dat({path_to_data}, files, SampleNames);
            tab_data = IO.get_table_data(files, dat, SampleNames);
            IO.func_load_execute(...
                [], ...
                app, dat, files, {path_to_data}, SampleNames{end}, tab_data ...
            );
            cd(current_path);
        end

        function load_3d_dataset_wsp(app)
            current_path = fileparts(mfilename('fullpath'));
            path_to_data = ['test_data' filesep '3D' filesep 'wsp' filesep];

            pnm = [current_path filesep path_to_data];
            fnm = '20190617_Settings2_LN1_9.wsp';
            
            IO_wsp.func_WSPload_backend(...
                app, ...
                pnm, ...
                fnm...
            )
        end
    end
end

