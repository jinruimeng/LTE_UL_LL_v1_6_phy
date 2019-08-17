%% config script for subsimulations

simulation.name = 'block_vs_fastfading';

simulation.parameter_scripts = {'LTE_UL_load_parameters'};

% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 1;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(-5,35,8);  % SNR points
simulation.N_subframes = 50;  
simulation.show_plot = true;


% different subsimulations
s1 = multisim.Config('BlockFading');
s1.add('ChanMod_config.filtering ', 'BlockFading');
simulation.add(s1);


s2 = multisim.Config('FastFading');
s2.add('ChanMod_config.filtering ', 'FastFading');
simulation.add(s2);



