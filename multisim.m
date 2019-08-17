%% multisim
function [simulation] = multisim(filename)

clc;

addpath(pwd); % needed for classes to be visible, atm it's necessary that you are the files path




simulation = multisim.Simulation;


run(filename); % define your parameters here


all_variables = whos;
var_names = {all_variables.name};
N_vars = length(var_names);

sim_configs = {};

for s = 1:N_vars
    tmp = eval(var_names{s});
    if isa(tmp, 'multisim.Config')
        sim_configs{end+1} = tmp; % maybe optimize this..
    end
end

if isempty(simulation.sim_configs)
    simulation.sim_configs = sim_configs;
    warning('please add subsimulations explicitly by using simulation.add()');
end
     
N_sim = length(simulation.sim_configs);

DEBUG_LEVEL = simulation.DEBUG_LEVEL;
cqi_i = simulation.cqi_i;          
SNR_vec = simulation.SNR_vec;  
N_subframes = simulation.N_subframes;  

if length(SNR_vec) > 1 && length(simulation.sweep_values) > 1
    error('only one SNR value is allowed in parameter sweep mode');
end

% check params
if ischar(simulation.parameter_scripts) % only one script given
        run(simulation.parameter_scripts)
else
    run(simulation.parameter_scripts{1});
end

for s = 1:N_sim
    sim_conf = simulation.sim_configs{s};
    params = sim_conf.parameters.keys;
    

    for p=1:length(params)
        try
            eval(['LTE_params.',params{p},';']);
        catch
            error('non existing parameter: %s', params{p});
        end
        
    end
end

for s = 1:N_sim
    sim_conf = simulation.sim_configs{s};
    sim_result = multisim.Result;
    
    
    fprintf('simulating "%s" \n', sim_conf.label);
                            

    % set cqi if specified
    if ~isnan(sim_conf.cqi)
        cqi_i = sim_conf.cqi;
        LTE_params.cqi = sim_conf.cqi;
    end
    
    % load the parameterscripts
    if ischar(simulation.parameter_scripts) % only one script given
        run(simulation.parameter_scripts)
    else
        for i=1:length(simulation.parameter_scripts)
            run(simulation.parameter_scripts{i});
        end
    end
    
    
    
    % overwrite the variable parameter
    params = sim_conf.parameters.keys;
    for p=1:length(params)
        mcell = strsplit(params{p}, '.');
        LTE_params = setfield(LTE_params, mcell{:}, sim_conf.parameters(params{p}));
    end
    
    % hack to get multisim running for more different UE sizes
    LTE_params.connection_table =  utils.generate_connection_table(LTE_params.nBS, LTE_params.nUE); 
    LTE_params.pathloss_matrix = utils.generate_pathloss_matrix(LTE_params.connection_table, inf);
    
    
    if isempty(simulation.sweep_name) % no sweep parameter
        % Load dependent Parameters and generate Elements
        LTE_UL_load_parameters_dependent;
        LTE_UL_load_parameters_generate_elements;
        LTE_UL_check_parameters;

        sim_result.LTE_params = LTE_params; % save sim parameters

        LTE_UL_sim_main; % run the simulator
        
        if LTE_params.confidence_interval_probability > 0
            LTE_UL_calculate_confidence_intervals(simulation_results, 0.95);
        end
        sr = multisim.Result;
        sr.LTE_params = LTE_params;
        sr.simulation_results = simulation_results;

        simulation.subsimulations(sim_conf.label) = sr; 
        
    else % sweep parameters given
        
        if isempty(simulation.sweep_values)
            error('sweep_name given without a list of sweep_values.');
        end
        
        fprintf('--- parametersweep: --- "%s"\n', simulation.sweep_name);
        % debug print parameters?
        
        simulation.subsimulations(sim_conf.label) = [];
        
        for sv_=1:length(simulation.sweep_values)
            fprintf('"%s = %s"\n', simulation.sweep_name, simulation.sweep_values(sv_));
            
            mcell = strsplit(simulation.sweep_name, '.');
            LTE_params = setfield(LTE_params, mcell{:}, simulation.sweep_values(sv_));
            
            LTE_UL_load_parameters_dependent;
            LTE_UL_load_parameters_generate_elements;
            LTE_UL_check_parameters
            
            
            sim_result.LTE_params = LTE_params; % save sim parameters

            LTE_UL_sim_main; % run the simulator
            
            if LTE_params.confidence_interval_probability > 0
                LTE_UL_calculate_confidence_intervals(simulation_results, 0.95);
            end
            sr = multisim.Result;
            sr.LTE_params = LTE_params;
            sr.simulation_results = simulation_results;
            sr.sweep_name = simulation.sweep_name;
            sr.sweep_value = simulation.sweep_values(sv_);
            
            simulation.subsimulations(sim_conf.label) = [simulation.subsimulations(sim_conf.label), sr];
        end
        
        
    end




end

% save everything
dt = datestr(now,'_dd-mm-yy_HH-MM-SS');
full_filename = ['multisim_results/', strrep(simulation.name, ' ', '_'), dt];
save(full_filename,'simulation');

if simulation.show_plot
    multisim.plot_simulation2(simulation);
end


end





