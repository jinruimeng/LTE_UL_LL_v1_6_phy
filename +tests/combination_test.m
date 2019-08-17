%%
clear;clc;

paramset = { ...
    {{'ChanMod_config.filtering', 'UE_config.user_speed'}, {{'BlockFading', 0}, {'FastFading', 20}}}, ...
    {'ChanMod_config.time_correlation', {'independent', 'correlated'}}, ...
    {{'UE_config.PMI_fb', 'UE_config.RI_fb', 'scheduler.type', 'UE_config.mode', 'UE_config.nTX', 'BS_config.nRX', 'nUE', 'nBS', 'connection_table', 'pathloss_matrix'}, ...
        {{false, false, 'round robin', 1, 1, 1, 1, 1, 1, 0}, ...
        {false, false, 'round robin', 4, 2, 2, 1, 1, 1, 0}, ...
        {true, true, 'round robin', 4, 2, 2, 1, 1, 1, 0}, ...
        {true, true, 'greedy MU MIMO', 5, 1, 4, 4, 1, ones(1, 4), zeros(1, 4)}, ...
        {true, true, 'round robin', 1, 1, 1, 2, 2, utils.generate_connection_table(2, 2), utils.generate_pathloss_matrix(utils.generate_connection_table(2, 2), inf)}}}, ...
    {'UE_config.RI_fb', {true, false}} ...
    {'UE_config.channel_averaging', {true, false}}, ...
    {'UE_config.ignore_ISI_ICI', {true, false}}, ...
};
loadparamscript = 'LTE_UL_load_parameters';

N_par = length(paramset);
size_vec = zeros(1,N_par); % row per parameter entries size of choices

for ii = 1:N_par
    size_vec(ii) = length( paramset{ii}{2});
end

N_paramsims = prod(size_vec);

% generate all the combinations
tmp_c = cell(1,N_par);
for row_index = 1:N_par
    tmp_c{row_index} = 1:size_vec(row_index);
end
param_mat = combvec(tmp_c{:})'; % needs the neuronal networks toolbox...

testresults = cell(2, N_paramsims); % 1 success, 2 result, 3 errorlog

fprintf('running %d simulations\n', N_paramsims);
for ii=1:N_paramsims
    %fprintf('%d ----\n', i);
    
    tmps = multisim.Config(num2str(ii));
    
    for p=1:N_par
        key = paramset{p}{1};
        value = paramset{p}{2}{param_mat(ii,p)};
        
        if iscell(key)
            for k = 1:length(key)
                tmps.add(key{k}, value{k});
            end
        else
            tmps.add(key, value);
        end
    end
    
    % "multisim" part
    fprintf('test %d/%d... \n', ii,N_paramsims);
                            
    try
        DEBUG_LEVEL = 0;
        cqi_i = 15; 
        SNR_vec = linspace(-5,35,4);  % SNR points
        N_subframes = 3;  
        show_plot = false;
        
        % load the parameterscript
        run(loadparamscript)

        % overwrite the variable parameter
        params = tmps.parameters.keys;
        for p=1:length(params)
            mcell = strsplit(params{p}, '.');
            LTE_params = setfield(LTE_params, mcell{:}, tmps.parameters(params{p}));
        end
        
        % This many plots would crash Matlab
        LTE_params.show_plots = false;

        % Load dependent Parameters and generate Elements
        run('LTE_UL_load_parameters_dependent');
        run('LTE_UL_load_parameters_generate_elements');
        run('LTE_UL_check_parameters');

        sim_result.LTE_params = LTE_params; % save sim parameters

        LTE_UL_sim_main; % run the simulator

        LTE_UL_calculate_confidence_intervals(simulation_results, 0.95);

        testresults{1, ii} = true; % success
        testresults{2, ii} = simulation_results;
     
    catch tmp_error
        % rethrow(tmp_error);
        testresults{1, ii} = false;
        testresults{2, ii} = tmp_error;
    end
end

clc;

passed = 0;
failed_tests = {};
 
for ii = 1:N_paramsims
    if testresults{1, ii}
        passed = passed + 1;
    else
        failed_tests{end+1} = ii;
    end
end

fprintf('%3d/%3d tests passed\n', passed, N_paramsims);

if passed < N_paramsims
    fprintf('failed tests:\n');
    disp(failed_tests);
end