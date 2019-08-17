%% runs a couple of tests 
% (c) Lukas Nagel, ITC, 2016

clear; clc;

test_scripts = { ...
    'test_channels1', ...
    'test_2x2', ...
    'test_2x4', ...
    'test_4x4', ...
    'test_TR_36_873'
    };

result = tests.runtests(test_scripts);

clc;
tests.display_result(result);
