%% config script for subsimulations

simulation.name = 'iammmse_test';

simulation.parameter_scripts = {'LTE_UL_load_parameters'};

% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 1;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(0,40,8);  % SNR points
simulation.N_subframes = 30;  
simulation.show_plot = false;


% different subsimulations
s1 = multisim.Config('2 BS');
s1.add('nBS', 2);
simulation.add(s1);


s2 = multisim.Config('3 BS');
s2.add('nBS', 3);
simulation.add(s2);

s4 = multisim.Config('4 BS');
s4.add('nBS', 4);
simulation.add(s4);


