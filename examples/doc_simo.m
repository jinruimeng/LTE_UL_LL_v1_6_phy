%% config script for subsimulations

simulation.name = 'doc simo';

simulation.parameter_scripts = {'examples/load_params_doc_antenna'};


% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 15;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(-10,40,12);  % SNR points
simulation.N_subframes = 50;
simulation.show_plot = false;

s1 = multisim.Config('1x1');
s1.add('UE_config.nTX', 1); 
s1.add('BS_config.nRX', 1); 
simulation.add(s1);

s2 = multisim.Config('1x2');
s2.add('UE_config.nTX', 1); 
s2.add('BS_config.nRX', 2); 
simulation.add(s2);


s3 = multisim.Config('1x4');
s3.add('UE_config.nTX', 1); 
s3.add('BS_config.nRX',  4); 
simulation.add(s3);


