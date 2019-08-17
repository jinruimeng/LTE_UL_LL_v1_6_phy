%% config script for subsimulations

simulation.name = 'doc mimo';

simulation.parameter_scripts = {'examples/load_params_doc_antenna'};


% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 15;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(0,40,8);  % SNR points
simulation.N_subframes = 20;
simulation.show_plot = false;

s1 = multisim.Config('CLSM');
s1.add('UE_config.nTX', 4); 
s1.add('BS_config.nRX', 4); 
s1.add('UE_config.mode', 4); 
s1.add('nUE', 4);
simulation.add(s1);

s2 = multisim.Config('MUMIMO');
s2.add('UE_config.nTX', 1); 
s2.add('BS_config.nRX', 4); 
s2.add('UE_config.mode', 5);
s2.add('nUE', 4);
s2.add('scheduler.type', 'random MU MIMO'); %'greedy MU MIMO');
simulation.add(s2);

