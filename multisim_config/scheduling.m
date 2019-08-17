%% config script for subsimulations

simulation.name = 'scheduling';

simulation.parameter_scripts = {'LTE_UL_load_parameters'};

% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 1;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(0,40,12);  % SNR points
simulation.N_subframes = 100;  



% different subsimulations
s1 = multisim.Config('Round Robin');
s1.add('scheduler.type', 'round robin');
simulation.add(s1);


s2 = multisim.Config('opt. max. thorughput');
s2.add('scheduler.type', 'opt max throughtput');
simulation.add(s2);


s3 = multisim.Config('approx. max. thorughput');
s3.add('scheduler.type', 'approx max throughtput');
simulation.add(s3);




