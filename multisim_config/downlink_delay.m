%% config script for subsimulations

simulation.name = 'downlinkdelayforstefan';

simulation.parameter_scripts = {'LTE_UL_load_parameters'};

% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 15;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(0,40,12);  % SNR points
simulation.N_subframes = 1000;  

% % channel estimation
% s1 = multisim.Config('1 delay, perfect CSI');
% s1.add('downlink_delay', 1);
% s1.add('BS_config.channel_estimation_method', 'PERFECT');
% 
% s2 = multisim.Config('1 delay, LS, ignore MSE');
% s2.add('downlink_delay', 1);
% s2.add('BS_config.channel_estimation_method', 'LS_SAV');
% s2.add('UE_config.ignore_channel_estimation', true);
% 
% s3 = multisim.Config('1 delay, LS');
% s3.add('downlink_delay', 1);
% s3.add('BS_config.channel_estimation_method', 'LS_SAV');
% s3.add('UE_config.ignore_channel_estimation', false);
% 
% s4 = multisim.Config('1 delay, LS, perfect CSI FB');
% s4.add('downlink_delay', 1);
% s4.add('BS_config.channel_estimation_method', 'LS_SAV');
% s4.add('UE_config.ignore_channel_estimation', false);
% s4.add('UE_config.MCS_and_scheduling_CSI', 'perfect');

% perfect channel knowledge
s1 = multisim.Config('0 delay');
s1.add('downlink_delay', 0);
simulation.add(s1);

s2 = multisim.Config('1 delay');
s2.add('downlink_delay', 1);
simulation.add(s2);


s3 = multisim.Config('2 delay');
s3.add('downlink_delay', 2);
simulation.add(s3);


s4 = multisim.Config('3 delay');
s4.add('downlink_delay', 3);
simulation.add(s4);

s5 = multisim.Config('4 delay');
s5.add('downlink_delay', 4);
simulation.add(s5);

s6 = multisim.Config('5 delay');
s6.add('downlink_delay', 5);
simulation.add(s6);





