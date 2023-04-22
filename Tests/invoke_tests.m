function invoke_tests(fnm, functions)
    disp('inside of invoke tests')
    warning('off','all');
    addpath(['..' filesep 'CytoMAP' filesep]);
    addpath(['..' filesep 'CytoMAP' filesep '3rdPartyFunctions']);


    if ~exist('functions', 'var')
        functions = [];
    end
    
    disp("Test " + string(fnm) + " started.");
    if ~invoke_tests_helper(fnm, functions)
        disp('Test Failed on test');
        exit(1);
    else
        disp('Test Completed');
        exit(0);
    end
end

function success = invoke_tests_helper(fnm, functions)
    disp('inside of invoke test helper')
    import_name = split(fnm, '.');
    import_name = import_name{1};
    import([import_name '.*']);

    success = true;

    fns = textread(fnm, '%s', 'delimiter', '\n', 'whitespace', '');
    fns = strtrim(fns);
    fns = fns(startsWith(fns, 'function'));
    fns = split(fns, 'function');
    if size(fns, 2) == 1
        fns = fns';
    end
    fns = strtrim(fns(:, 2));
    fns = fns(contains(fns, '('));

    if ~isempty(functions) && ~strcmp(functions, 'all')
        fns = fns(contains(fns, functions));
    end
    
    is_ret = contains(fns, '=');
    fns_ret = split(fns(is_ret), '=');
    if ~isempty(fns_ret)
        fns(is_ret) = strtrim(fns_ret(:, 2));
    end
    % FOR NOW NOT ACCEPT is_ret functions
    fns = fns(~is_ret);
    is_ret = contains(fns, '=');
    
    test_idx = startsWith(fns, 'test_');
    fns = fns(test_idx);
    is_ret = is_ret(test_idx);
    is_app = endsWith(fns, '(app)');

    fns = split(fns, '(');
    if size(fns, 2) == 1
        fns = fns';
    end
    fns = strtrim(fns(:, 1));
    for f_idx = 1:numel(fns)
        fn = fns{f_idx};
        disp("Running test " + string(fn))
        if is_app(f_idx) && is_ret(f_idx)
            app = CytoMAP;
            try
                out = eval(strcat(import_name, '.', fn, '(app);'));
            catch ME
                disp(getReport(ME, 'extended'));
                success = false;
            end
            Helper.func_exit(app);
            if ~success || ~out
                success = false;
                disp(strcat('Failed on test', fn));
                break;
            end
        elseif is_app(f_idx)
            app = CytoMAP;
            try
                eval(strcat(import_name, '.', fn, '(app);'));
            catch ME
                disp(getReport(ME, 'extended'));
                success = false;
            end
            Helper.func_exit(app);
            if ~success
                disp(strcat('Failed on test', fn));
                break;
            end
        elseif is_ret(f_idx)
            try
                out = eval(strcat(import_name, '.', fn, '();'));
            catch ME
                disp(getReport(ME, 'extended'));
                success = false;
            end
            if ~success || ~out
                success = false;
                disp(strcat('Failed on test', fn));
                break;
            end
        else
            try
                eval(strcat(import_name, '.', fn, '();'));
            catch ME
                disp(getReport(ME, 'extended'));
                success = false;
            end
            if ~success
                success = false;
                disp(strcat('Failed on test', fn));
                break;
            end
        end
        disp("Test " + string(fn) + " succeeded"); 
    end
end