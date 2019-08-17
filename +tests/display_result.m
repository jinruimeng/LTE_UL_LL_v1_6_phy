function  display_result( result )
%DISPLAY_RESULT gives an overview of the tests in result
% (c) Lukas Nagel, ITC, 2016

tests = result.tests;
N = length(tests);
passed = 0;

for i = 1:N
    if tests{i}.success 
        passed_phrase = 'passed';
        passed = passed + 1;
    else 
        passed_phrase = 'failed';
    end
    
    if length(tests{i}.name) > 15
        name = tests{i}.name(1:15);
    else
        name = tests{i}.name;
    end
    
    fprintf('%10s: %s\n', name, passed_phrase);
end

% summary

fprintf('-----------------------------------\n');
fprintf('%2d/%2d tests passed\n', passed, N);

