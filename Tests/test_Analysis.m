classdef test_Analysis
    methods (Static)
        function test_dist_empty(app)
            %% Load data in
            data_loader.load_2d_dataset_mat(app, 'b');
            Analysis.dist_func(app);

            %% Default test, only phenotypes
            % 1 sample
            Analysis.dist_func_backend(app, {'all_1'}, {'All/T Cells OT II CXCR3+', 'All/T Cells OT-II 4GET+'});
            assert(all(ismember( ...
                {[Constants.other_tag 'DistTo_All_T_Cells_OT_II_CXCR3POS'], [Constants.other_tag 'DistTo_All_T_Cells_OTNEGII_4GETPOS']}, ...
                app.data.all_1.AllCells.Properties.VariableNames ...
            )));
            assert(all(~ismember( ...
                {[Constants.other_tag 'DistTo_T_cells_OTNEGII_T_Cells'], [Constants.other_tag 'DistTo_All_T_Cells_OTNEGII_CXCR3PO']}, ...
                app.data.all_1.AllCells.Properties.VariableNames ...
            )));

            % Multiple samples/consistency-overwriting
            dist_cxcr3_before = app.data.all_1.AllCells{:, {[Constants.other_tag 'DistTo_All_T_Cells_OT_II_CXCR3POS']}};
            Analysis.dist_func_backend(app, {'all_1', 'all_2'}, {'All/T Cells OT II CXCR3+', 'All/T cells OT-II T Cells', 'All/T Cells OT-II CXCR3+ 4get-'});
            dist_cxcr3_after  = app.data.all_1.AllCells{:, {[Constants.other_tag 'DistTo_All_T_Cells_OT_II_CXCR3POS']}};
            assert(all(dist_cxcr3_before == dist_cxcr3_after));
            assert(all(ismember( ...
                { ...
                    [Constants.other_tag 'DistTo_All_T_Cells_OT_II_CXCR3POS'], ...
                    [Constants.other_tag 'DistTo_All_T_Cells_OTNEGII_4GETPOS'], ...
                    [Constants.other_tag 'DistTo_All_T_cells_OTNEGII_T_Cells'], ...
                    [Constants.other_tag 'DistTo_All_T_Cells_OTNEGII_CXCR3POS_4getNEG'] ...
                }, ...
                app.data.all_1.AllCells.Properties.VariableNames ...
            )));
            assert(all(ismember( ...
                {[Constants.other_tag 'DistTo_All_T_cells_OTNEGII_T_Cells'], [Constants.other_tag 'DistTo_All_T_Cells_OTNEGII_CXCR3POS_4getNEG']}, ...
                app.data.all_2.AllCells.Properties.VariableNames ...
            )));

            %% Test Spots, Regions, Polygons, Surface etc.
        end
    end
end
