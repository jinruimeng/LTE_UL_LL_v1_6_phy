function [ new_simulation ] = combine_to_multisim( name, sims )
% combines simulation_results of single simulations to a multisim
% simulation


labels = sims.keys;
paths = sims.values;

N_sims = length(labels);

new_simulation = multisim.Simulation;
new_simulation.name = name;





for i__ = 1:N_sims
    
    load(paths{i__});
    
    % maybe only take one 
    new_simulation.SNR_vec = simulation_results.SNR_vector;  
    new_simulation.N_subframes = simulation_results.N_subframes;  
    

    
     sr = multisim.Result;
     sr.LTE_params = LTE_params;
     sr.simulation_results = simulation_results;

     new_simulation.subsimulations(labels{i__}) = sr;

end


end












