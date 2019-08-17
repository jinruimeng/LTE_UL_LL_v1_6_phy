%% config script for subsimulations

simulation.name = '2x2';

simulation.parameter_scripts = {'test_scenarios/test_load_parameters'};



% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 15;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(-5,35,4);  % SNR points
simulation.N_subframes = 3;  
simulation.show_plot = false;


% different subsimulations
s1 = multisim.Config('2x2');
s1.add('UE_config.mode', 4); 
s1.add('UE_config.nTX', 2); 
s1.add('BS_config.nRX', 2); 






