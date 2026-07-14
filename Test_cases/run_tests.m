function run_tests()
    fprintf('=== Running 1D magma simulation test cases===\n\n');
    
    % Define all test cases - now store the folder name only
    test_cases = {
        'Pure_compaction'
    };
    
    passed = 0;
    failed = 0;
    
    % Get the current directory where run_tests.m is located
    current_dir = fileparts(mfilename('fullpath'));
    
    for i = 1:length(test_cases)
        test_name = test_cases{i};
        test_folder = fullfile(current_dir, test_name);
        
        fprintf('Test: %s... ', test_name);
        
        % Check if the folder exists
        if ~exist(test_folder, 'dir')
            error('Test folder not found: %s', test_folder);
        end
        
        % Add the test folder to path
        addpath(test_folder);
        
        % Run the test function (function name matches folder name)
        try
            func = str2func(test_name);
            [Results, err] = func();
        catch ME
            fprintf('ERROR: %s\n', ME.message);
            failed = failed + 1;
            rmpath(test_folder);
            continue;
        end
        
        % Remove from path
        rmpath(test_folder);
        
        % Check results
        if all(Results) && max(err) < 0.03
            passed = passed + 1;
            fprintf('PASSED\n');
        else
            fprintf('FAILED\n');
            failed = failed + 1;
        end
    end
    
    fprintf('\n=== Results: %d/%d tests passed ===\n', passed, passed + failed);
    
    % Exit with error if any test failed
    if failed > 0
        error('%d test(s) failed!', failed);
    end
end