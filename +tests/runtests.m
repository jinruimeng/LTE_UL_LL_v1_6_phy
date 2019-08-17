function [ test_results ] = runtests( scripts )
%RUNTESTS runs all config scripts provided in "scripts" with multisim, 
% catches the errors and add saves it into test_results
%
% (c) Lukas Nagel, ITC, 2016

test_results = struct;
test_results.tests = cell(1, length(scripts));

for i = 1:length(scripts)
    test_result = struct;
    test_result.success = true;
    test_result.name = scripts{i};
    
    try
        %eval(['tests.', scripts{i}])
        test_result.result = multisim(['test_scenarios/', scripts{i}]);
    catch tmp_error
        test_result.error = tmp_error;
        test_result.success = false;
    end
    
    test_results.tests{i} = test_result;
end

end

