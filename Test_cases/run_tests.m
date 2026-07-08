% tests/run_tests.m
function run_tests()
    fprintf('=== Running 1D magma simulation test cases===\n\n');
    
    % Define all test cases
    test_cases = {
        'Pure_compaction', 'Pure_compaction.m'
    };
    
    passed = 0;
    failed = 0;
    
    for i = 1:size(test_cases, 1)
        test_name = test_cases{i, 1};
        expected_file = test_cases{i, 2};
        
        fprintf('Test: %s... ', test_name);
        
        
        addpath(test_cases{i});
        % Run the settings file to set parameters
        func=str2func (test_cases{i});
        [Results, err]=func();
        rmpath(test_cases{i});

        if all(Results) && max(err)<0.03
            passed=passed+1;
        else
            failed=failed+1;
        end    
    end
    
    fprintf('\n=== Results: %d/%d tests passed ===\n', passed, passed + failed);
    
    % Exit with error if any test failed
    if failed > 0
        error('%d test(s) failed!', failed);
    end
end