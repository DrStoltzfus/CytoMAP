function generate_datasets_2()
    addpath('../../CytoMAP/')
    helper();
%     exit(0);
end

function helper()
    app = CytoMAP;

    current_path = pwd;
    path_to_data = 'Data3D/dataset_a/';
    path_out = 'Data3D';
    files = {strtrim(split(ls(path_to_data), '  '))'};

    cd(path_to_data);
    SampleNames = IO.get_sample_names(app, {path_to_data});
    dat = IO.get_dat({path_to_data}, files, SampleNames);
    tab_data = IO.get_table_data(files, dat, SampleNames);
    IO.func_load_execute(...
        [], ...
        app, dat, files, {path_to_data}, SampleNames{end}, tab_data ...
    );
    cd(current_path);

    % Save dataset B
    IO.func_save(app, [path_out '/dataset_b.mat']);

    % Create dataset C
    Analysis.raster_scan_backend(app, {'Sample_1'}, false);
    Analysis.cell_ctr_backend(app, {'Sample_1'}, app.data.Sample_1.GateTags{2, :}, 50);

    % Save dataset C
    IO.func_save(app, [path_out '/dataset_c.mat']);

    %% Create dataset D
    %% Make gates rectangular and polygon
    Plotting.func_newfig(app);
    
    Plt_Gating.new_gate(app, 'fig1', 'RectGate', ...
                            'name', 'rect_gate_x_y', ...
                            'position', [210, 241, 556, 549] ... 
    );
    Plt_Gating.new_gate(app, 'fig1', 'PolyGate', ...
                            'name', 'poly_gate_x_y', ...
                            'position', [292,  91;
                                          97, 251;
                                         292, 411;
                                         746, 302;
                                         746, 200] ... 
    );

    %% Make Surfaces on default axes
    [dat, ~, ~] = Plt_Helper.UpdateFigMenu(app, 'fig1');
    for smpl_idx=1:numel(app.figOpts.smpls.fig1)
        smpl = app.figOpts.smpls.fig1{smpl_idx};
        surf_name = "all_x_y_z";
        Plt_Helper.surf_create( ...
            app, smpl, app.figOpts.phenos.fig1, dat, ...
            app.figOpts.xaxIF.fig1.String{app.figOpts.xaxIF.fig1.Value}, ...
            app.figOpts.yaxIF.fig1.String{app.figOpts.yaxIF.fig1.Value}, ...
            app.figOpts.zaxIF.fig1.String{app.figOpts.zaxIF.fig1.Value}, ...
            100, 'UDS', surf_name, 50, 50 ...
        );
    end

    %% Switch to multiple phenotypes and different axes
    % TODO: It would be nice to have a programative way of changing those things in Plotting (but that's like Ver 3 stuff)
    app.figOpts.phenos.fig1 = app.data.Sample_1.GateTags{2, :};
    Plotting.func_plot(app, 'fig1');
    app.figOpts.phnTMP.fig1 = app.figOpts.phenos.fig1;

    app.figOpts.xaxIF.fig1.Value = find(strcmp(app.figOpts.xaxIF.fig1.String, 'A 0'));
    Plotting.func_plot(app, 'fig1');
    app.figOpts.xaxIFTMP.fig1 = app.figOpts.xaxIF.fig1.String{app.figOpts.xaxIF.fig1.Value};

    app.figOpts.yaxIF.fig1.Value = find(strcmp(app.figOpts.yaxIF.fig1.String, 'A 9'));
    app.figOpts.yaxIFTMP.fig1 = 'Y';
    Plotting.func_plot(app, 'fig1');
    app.figOpts.yaxIFTMP.fig1 = app.figOpts.yaxIF.fig1.String{app.figOpts.yaxIF.fig1.Value};

    app.figOpts.zaxIF.fig1.Value = find(strcmp(app.figOpts.yaxIF.fig1.String, 'A 5'));
    app.figOpts.zaxIFTMP.fig1 = 'Z';
    Plotting.func_plot(app, 'fig1');
    app.figOpts.zaxIFTMP.fig1 = app.figOpts.yaxIF.fig1.String{app.figOpts.yaxIF.fig1.Value};

    Plt_Gating.new_gate(app, 'fig1', 'RectGate', ...
                            'name', 'rect_gate_a0_a9', ...
                            'position', [114, 0.0001, 286, 0.0007] ... 
    );
    Plt_Gating.new_gate(app, 'fig1', 'PolyGate', ...
                            'name', 'poly_gate_a0_a9', ...
                            'position', [904, 0.0003
                                         659, 0.0009
                                          92, 0.0001] ... 
    );
    Plt_Gating.save_gate(app, 'fig1');

    %% Make Surfaces on non-default axes
    [dat, ~, ~] = Plt_Helper.UpdateFigMenu(app, 'fig1');
    for smpl_idx=1:numel(app.figOpts.smpls.fig1)
        smpl = app.figOpts.smpls.fig1{smpl_idx};
        surf_name = "multi_a0_a9_a5";
            
        Plt_Helper.surf_create( ...
            app, smpl, app.figOpts.phenos.fig1, dat, ...
            app.figOpts.xaxIF.fig1.String{app.figOpts.xaxIF.fig1.Value}, ...
            app.figOpts.yaxIF.fig1.String{app.figOpts.yaxIF.fig1.Value}, ...
            app.figOpts.zaxIF.fig1.String{app.figOpts.zaxIF.fig1.Value}, ...
            100, 'UDS', surf_name, 50, 50 ...
        );
    end

    %% Make neighborhoods
    app.figOpts.phenos.fig1 = {'Density/MFI RSN'};
    Plotting.func_plot(app, 'fig1');
    app.figOpts.phnTMP.fig1 = app.figOpts.phenos.fig1;

    Plt_Gating.new_gate(app, 'fig1', 'RectGate', ...
                            'name', 'rect_x_y', ...
                            'position', [422, 325, 350, 500] ... 
    );
    Plt_Gating.new_gate(app, 'fig1', 'PolyGate', ...
                            'name', 'poly_x_y', ...
                            'position', [136, 173;
                                         837, 187;
                                         674, 524;
                                         176, 620] ...
    );

    %% Make Surfaces on MFIRSN
    [dat, ~, ~] = Plt_Helper.UpdateFigMenu(app, 'fig1');
    for smpl_idx=1:numel(app.figOpts.smpls.fig1)
        smpl = app.figOpts.smpls.fig1{smpl_idx};
        surf_name = "mfirsn_x_y_z";
        Plt_Helper.surf_create( ...
            app, smpl, app.figOpts.phenos.fig1, dat, ...
            app.figOpts.xaxIF.fig1.String{app.figOpts.xaxIF.fig1.Value}, ...
            app.figOpts.yaxIF.fig1.String{app.figOpts.yaxIF.fig1.Value}, ...
            app.figOpts.zaxIF.fig1.String{app.figOpts.zaxIF.fig1.Value}, ...
            100, 'UDS', surf_name, 50, 50 ...
        );
    end

    %% Switch to multiple phenotypes and different axes
    app.figOpts.phenos.fig1 = {'Density/MFI CCN'};
    Plotting.func_plot(app, 'fig1');
    app.figOpts.phnTMP.fig1 = app.figOpts.phenos.fig1;

    app.figOpts.xaxIF.fig1.Value = find(strcmp(app.figOpts.xaxIF.fig1.String, 'A 0'));
    Plotting.func_plot(app, 'fig1');
    app.figOpts.xaxIFTMP.fig1 = app.figOpts.xaxIF.fig1.String{app.figOpts.xaxIF.fig1.Value};

    app.figOpts.yaxIF.fig1.Value = find(strcmp(app.figOpts.yaxIF.fig1.String, 'A 9'));
    app.figOpts.yaxIFTMP.fig1 = 'Y';
    Plotting.func_plot(app, 'fig1');
    app.figOpts.yaxIFTMP.fig1 = app.figOpts.yaxIF.fig1.String{app.figOpts.yaxIF.fig1.Value};
    Plt_Gating.new_gate(app, 'fig1', 'RectGate', ...
                            'name', 'rect_a0_a9', ...
                            'position', [320, 0.00025, 350, 0.00045] ... 
    );
    Plt_Gating.new_gate(app, 'fig1', 'PolyGate', ...
                            'name', 'poly_a0_a9', ...
                            'position', [700, 0.00015
                                         500, 0.00085
                                         150, 0.00025] ... 
    );
    Plt_Gating.save_gate(app, 'fig1');

    %% Make Spots
    Plt_Helper.spots_save([250, 350, 450], [100, 782, 123], [751, 456, 545], app, []);

    %% Make Surfaces on MFICCN
    [dat, ~, ~] = Plt_Helper.UpdateFigMenu(app, 'fig1');
    for smpl_idx=1:numel(app.figOpts.smpls.fig1)
        smpl = app.figOpts.smpls.fig1{smpl_idx};
        surf_name = "mficcn_a0_a9_a5";
        Plt_Helper.surf_create( ...
            app, smpl, app.figOpts.phenos.fig1, dat, ...
            app.figOpts.xaxIF.fig1.String{app.figOpts.xaxIF.fig1.Value}, ...
            app.figOpts.yaxIF.fig1.String{app.figOpts.yaxIF.fig1.Value}, ...
            app.figOpts.zaxIF.fig1.String{app.figOpts.zaxIF.fig1.Value}, ...
            100, 'UDS', surf_name, 50, 50 ...
        );
    end

    %% Save dataset D
    IO.func_save(app, [path_out '/dataset_d.mat']);

    Helper.func_exit(app);
end