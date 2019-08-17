classdef Result < handle
    % Result contains the simulation results and parameters
    
    properties
        simulation_results
        LTE_params
        sweep_name
        sweep_value
        confidence
    end
    
    methods
        function SR = Result()
            SR.simulation_results = struct;
            SR.LTE_params = struct;
            
            sweep_name = '';
            sweep_value = '';
        end
    end
    
end

