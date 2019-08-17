%% shows how to combine stuff with combine_to_multisim

sims = containers.Map();
sims('TU') = 'test_TU~all.mat';
sims('VehB') = 'test_VehB~all.mat';

simulation = multisim.combine_to_multisim( 'two channels', sims );

multisim.plot_simulation2(simulation)
