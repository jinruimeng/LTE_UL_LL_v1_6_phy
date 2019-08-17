%% config script for subsimulations

simulation.name = 'rank_adaptation';

simulation.parameter_scripts = {'LTE_UL_load_parameters'};

% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 1;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(-5,35,12);  % SNR points
simulation.N_subframes = 2000;  



% different subsimulations
s1 = multisim.Config('rank 1');
s1.add('UE_config.RI', '1');
s1.add('UE_config.PMI_fb', true);
s1.add('UE_config.RI_fb', false);
s1.add('UE_config.CQI_fb', true);

s2 = multisim.Config('rank 2');
s2.add('UE_config.RI', '2');
s2.add('UE_config.PMI_fb', true);
s2.add('UE_config.RI_fb', false);
s2.add('UE_config.CQI_fb', true);

s3 = multisim.Config('rank 3');
s3.add('UE_config.RI', '3');
s3.add('UE_config.PMI_fb', true);
s3.add('UE_config.RI_fb', false);
s3.add('UE_config.CQI_fb', true);

s4 = multisim.Config('rank 4');
s4.add('UE_config.RI', '4');
s4.add('UE_config.PMI_fb', true);
s4.add('UE_config.RI_fb', false);
s4.add('UE_config.CQI_fb', true);

s5 = multisim.Config('rank adaptive');
s5.add('UE_config.PMI_fb', true);
s5.add('UE_config.RI_fb', true);
s5.add('UE_config.CQI_fb', true);


