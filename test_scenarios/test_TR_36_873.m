%% config script for subsimulations
simulation.name = 'TR 36.873 - 3D Model';
simulation.parameter_scripts = {'test_scenarios/test_TR_36_873_load_parameters'};

%% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 15;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = [-20, 0, 20];  % SNR points
simulation.N_subframes = 5;
simulation.show_plot = false;

%% Model specific
s1 = multisim.Config('Urban Macro');
s1.add('TR_36_873.environment', 'UMa');
simulation.add(s1);

s2 = multisim.Config('Urban Micro');
s2.add('TR_36_873.environment', 'UMi');
simulation.add(s2);

s3 = multisim.Config('LOS, indoor state and heights random');
s3.add('TR_36_873.LOS_according_model', true);
s3.add('TR_36_873.indoor_according_model', true);
s3.add('TR_36_873.UE_heights_according_to_model', true);
s3.add('TR_36_873.eNodeB_heights_according_to_model', true);
s3.add('nUE', 50);
s3.add('UE_config.mode', 1);
s3.add('UE_config.nTX', 1);
s3.add('BS_config.nRX', 1);
simulation.add(s3);

s4 = multisim.Config('LOS');
s4.add('TR_36_873.UE_is_LOS', true);
simulation.add(s4);

s5 = multisim.Config('NLOS');
s5.add('TR_36_873.UE_is_LOS', false);
simulation.add(s5);

s6 = multisim.Config('Indoor');
s6.add('TR_36_873.UE_is_indoor', true);
simulation.add(s6);

s7 = multisim.Config('Outdoor');
s7.add('TR_36_873.UE_is_indoor', false);
simulation.add(s7);

s8 = multisim.Config('XPOL antennas with 3D model antenna');
s8.add('TR_36_873.UE_antenna_polarization', 'XPOL');
s8.add('TR_36_873.antenna.antenna_polarization', 'XPOL');
s8.add('TR_36_873.antenna.antenna_gain_pattern', 'TR36.873 3D antenna');
simulation.add(s8);

%% Simulator specific
s9 = multisim.Config('FastFading');
s9.add('ChanMod_config.filtering', 'FastFading');
simulation.add(s9);

s10 = multisim.Config('Multiple eNodeBs and multiple users');
s10.add('nUE', 4);
s10.add('nBS', 2);
s10.add('TR_36_873.LOS_according_model', true);
s10.add('TR_36_873.indoor_according_model', true);
s10.add('TR_36_873.UE_heights_according_to_model', true);
s10.add('TR_36_873.eNodeB_heights_according_to_model', true);
s10.add('TR_36_873.antenna_azimuth_offset', zeros(1, 2));
simulation.add(s10);

s11 = multisim.Config('Large bandwidth');
s11.add('Bandwidth', 20e6);
s11.add('nUE', 1);
s11.add('nBS', 1);
s11.add('UE_config.nTX', 1);
s11.add('BS_config.nRX', 1);
s11.add('UE_config.mode', 1);
simulation.add(s11);




