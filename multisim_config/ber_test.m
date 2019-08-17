%% config script for subsimulations

simulation.name = 'BER PAPER FINAL';

simulation.parameter_scripts = {'LTE_UL_load_parameters_ber_erich'};

% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 4;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(-5,35,16);  % SNR points
simulation.N_subframes = 100;  



% different subsimulations
s1 = multisim.Config('OFDM');
s1.add('ChanMod_config.type', 'PedB');
s1.add('DFT_spreading_off', true');

s2 = multisim.Config('SC-FDMA');
s2.add('ChanMod_config.type', 'PedB');  
s2.add('DFT_spreading_off', false);




