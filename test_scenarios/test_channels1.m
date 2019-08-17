%% config script for subsimulations

simulation.name = 'channels1';

simulation.parameter_scripts = {'test_scenarios/test_load_parameters'};



% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 15;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(-5,35,4);  % SNR points
simulation.N_subframes = 3;  
simulation.show_plot = false;


% different subsimulations
s1 = multisim.Config('PedA');
s1.add('ChanMod_config.type', 'PedA'); 

s2 = multisim.Config('PedB');
s2.add('ChanMod_config.type', 'PedB');   







