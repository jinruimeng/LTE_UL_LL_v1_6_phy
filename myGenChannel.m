function [channlData_alice,channlData_bob] = myGenChannel(N_subframes,type,user_speed,cqi_vec,SNR_vec)

    %% DEBUG level
    global DEBUG_LEVEL
    DEBUG_LEVEL = 1;    % 0: no text output at all
    % 1: shows current simulation status
    % 2: shows additional variables and figures

    %% Actual simulations
    if nargin<4
        cqi_vec = 1;      % only used for static scheduling
    end

    if nargin<5
        SNR_vec= 10;
        %sim_result_temp.UE_specific.throughput_coded =[];
    end


    for cqi_i = cqi_vec

        %N_subframes = k;%number of subframe

        % load Simulation Parameters
        LTE_UL_load_parameters;
        % test data intialization

        % Load dependent Parameters and generate Elements
        LTE_UL_load_parameters_dependent;
        LTE_UL_load_parameters_generate_elements;
        LTE_UL_check_parameters;
        % actual Simulation
        tic
        LTE_UL_sim_main;
        toc

    end

    channlData_alice =  H_test_results_alice';
    channlData_bob =  H_test_results_bob';

end

