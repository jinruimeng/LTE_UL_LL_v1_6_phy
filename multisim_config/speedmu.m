%% config script for subsimulations

simulation.name = 'speedmu';

simulation.parameter_scripts = {'LTE_UL_load_parameters'};

% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 1;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(0,40,8);  % SNR points
simulation.N_subframes = 20;  



% different subsimulations
s1 = multisim.Config('0 kmh');
s1.add('UE_config.user_speed', 0/3.6);
simulation.add(s1);


s2 = multisim.Config('50 kmh');
s1.add('UE_config.user_speed', 50/3.6);
simulation.add(s2);



