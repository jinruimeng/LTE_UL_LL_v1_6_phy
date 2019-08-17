classdef Simulation < handle
    %Simulation wrapper for config and simulation results
    
    properties
        name
        subsimulations
        parameter_scripts
        
        DEBUG_LEVEL
        cqi_i   
        SNR_vec
        N_subframes
        sweep_name
        sweep_values
        sweep_label
        show_plot
        labels_in_order
        sim_configs
        
    end
    
    methods
        function s = Simulation()
            s.name = 'default';
            s.subsimulations  = containers.Map;
            
            
             % default values
            s.parameter_scripts = {'LTE_UL_load_parameters'};
            s.DEBUG_LEVEL = 1;
            s.cqi_i = 15;   % irrelevant if scheduler.assignment = 'dynamic'
            s.SNR_vec = linspace(-5,35,8);  % SNR points
            s.N_subframes = 100;
            s.sweep_name = '';
            s.sweep_values = [];
            s.sweep_label = '';
            s.show_plot = true;
            s.labels_in_order = {};
            s.sim_configs = {};
        end
        
        function add(simulation, subsim)
             simulation.sim_configs{end+1} = subsim;
             simulation.labels_in_order{end+1} = subsim.label;
        end
    end
    
end

