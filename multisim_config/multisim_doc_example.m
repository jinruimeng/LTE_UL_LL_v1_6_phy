%% config script for subsimulations
simulation.name = 'multisim_example';
simulation.parameter_scripts = {'multisim_config/load_parameters_multisim_doc'};

% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 1;   
simulation.SNR_vec = linspace(-5,35,8);  
simulation.N_subframes = 400;
simulation.show_plot = true;


%subsimulations
s1 = multisim.Config('PERFECT/ZF');
s1.add('BS_config.channel_estimation_method', 'PERFECT');
s1.add('BS_config.receiver', 'ZF');
simulation.add(s1);

s2 = multisim.Config('MMSE/MMSE');
s2.add('BS_config.channel_estimation_method', 'MMSE');
s2.add('BS_config.receiver', 'MMSE');
simulation.add(s2);




