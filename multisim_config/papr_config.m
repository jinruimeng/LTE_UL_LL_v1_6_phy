%%

%% config script for subsimulations

simulation.name = 'PAPR';

simulation.parameter_scripts = {'LTE_UL_load_parameters_papr'};

global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 7;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = 10;  % SNR points
simulation.N_subframes = 25000;  



% different subsimulations
s1 = multisim.Config('OFDM');
s1.add('DFT_spreading_off', true');
s1.add('cqi_i', 6);

s2 = multisim.Config('SC-FDM 4-QAM');
s2.add('DFT_spreading_off', false);
s2.add('cqi_i', 6);

s3 = multisim.Config('SC-FDM 16-QAM'); 
s3.add('DFT_spreading_off', false);
s3.add('cqi_i', 9);

s4 = multisim.Config('SC-FDM 64-QAM'); 
s4.add('DFT_spreading_off', false);
s4.add('cqi_i', 15);

s5 = multisim.Config('OFDM 10MHz');
s5.add('DFT_spreading_off', true');
s5.add('cqi_i', 6);
s5.add('Bandwidth', 10e6);

s6 = multisim.Config('SC-FDM 16-QAM 10MHz'); 
s6.add('DFT_spreading_off', false);
s6.add('cqi_i', 9);
s6.add('Bandwidth', 10e6);



