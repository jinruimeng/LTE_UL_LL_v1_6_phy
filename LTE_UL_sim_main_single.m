% LTE system simulator main simulator file. Check the LTE_sim_batch files
% to check how to launch the simulator.
% [] = LTE_sim_main()
% Author: Jan Prokopec, prokopec@feec.vutbr.cz, Jiri Blumenstein xblume00@phd.feec.vutbr
% (c) 2016 by ITC
% modification by DREL
% www.nt.tuwien.ac.at
% www.urel.feec.vutbr.cz

maxStreams = min(LTE_params.UE_config.nTX,LTE_params.BS_config.nRX);

simulation_results = results.simulationResultsUL( LTE_params );


%% SNR and frame loops
StartTime = clock;
% NOTE: changing this with a parfor does not change any functionality, but
% messes up the timer (it goes backwards), so if you really want parallel
% functionality, change the line below to a parfor

% SNR_vector = LTE_params.SNR_vec;
% BER_UNcoded = zeros( length(SNR_vec(1,:)), N_subframes);
% BER_coded = zeros( length(SNR_vec(1,:)), N_subframes);

nUE = LTE_params.nUE;
Ns = LTE_params.Ns;

for SNR_i=1:size(SNR_vec,2)
    % Initialize variables that will be reused, such as the BS_output
    BS_output = outputs.bsOutput(LTE_params.nUE, LTE_params.nBS, LTE_params.Nrb,maxStreams);
    subframe_i = 1;
    delay_counter = 0;
    
    delay = true;

    % Set network clock. It will tell every network element in which TTI we are right now.
    network_clock = network_elements.clock(LTE_params.FrameDur/10);


    % Initialize downlink channel
    downlinkChannel = channels.downlinkChannel(LTE_params.downlink_delay,LTE_params.nUE,LTE_params.nBS);
    
%     % Receive feedback from the previous subframe
%     UE_output = downlinkChannel.receive_feedback;
    
    % Attach the clock to each network element
    % Attach UE_output.genie
    % Initilize the traffic data models
    for u_=1:LTE_params.nUE*LTE_params.nBS
        UE(u_).clock = network_clock;
        UE_output(u_).UE_genie = outputs.genieInformation;
    end
    for b_=1:LTE_params.nBS
        BS(b_).clock = network_clock;
        BS(b_).realization_num = 0;
    end
    
    % Some further initialization
    SNR = SNR_vec(:,SNR_i);
    network_clock.reset;
    
    for bb = 1:LTE_params.nBS
        BS(bb).reset_HARQ_process_index;
    end
    
    ChanMod_output = cell(LTE_params.nBS,LTE_params.nUE*LTE_params.nBS);
    
    % ACK of the previous frame. If this is the first frame, set the
    % ACK to correct so that the HARQ handling generates new data
    if subframe_i==1 || (subframe_i <= LTE_params.downlink_delay && ~strcmp(LTE_params.ChanMod_config.time_correlation,'independent'))% (number of max HARQ processes that will be used)
        for uu = 1:LTE_params.nUE
            UE_output(uu).ACK    = true(1,2);
            UE_output(uu).rv_idx = zeros(1,2);
        end
    end
    

%     % This is for signalization only
%     downlinkChannel = channels.downlinkChannel(LTE_params.downlink_delay,LTE_params.nUE);
%
    if DEBUG_LEVEL > 0
        disp('');
        disp(['*************** SNR = ' num2str(SNR(1)) 'dB, value ' num2str(SNR_i) ' of ' num2str(size(SNR_vec,2)) ' ***************']);
    end
  H_test_results_alice=zeros(2400,N_subframes);
  H_test_results_bob=zeros(2400,N_subframes);
    %while subframe_i <= N_subframes
    for subframe_i = 1:N_subframes
        % First of all, advance the network clock
        network_clock.advance_1_TTI;
        if mod(subframe_i,50)==1 && ~delay
            if DEBUG_LEVEL > 0
                disp(['   processing subframe #' num2str(subframe_i) ' of ' num2str(N_subframes)])
                disp(['---> remaining simulation time: ' num2str(etime(clock,StartTime)/((SNR_i-1)*N_subframes+subframe_i-1)*((size(SNR_vec,2)-SNR_i)*N_subframes+N_subframes-subframe_i)/60,'%5.3f') 'min']);
                pause(0.05);
            end
        end
        
        [BS_output, UE_output,H_test_results] = LTE_UL_sim_main_subframe(LTE_params, ChanMod, SNR, UE, UE_output, BS, BS_output, subframe_i, out, downlinkChannel, cqi_i, channel); % channel is used for the winner model
       H_test_results_alice(:,subframe_i)=H_test_results{1,1};
       H_test_results_bob(:,subframe_i)=H_test_results{1,2};
        if delay_counter == LTE_params.downlink_delay
            delay = false;
        end       
        delay_counter = delay_counter + 1;
        
%         if ~delay
        if true
            simulation_results.process_TTI_results(BS_output,UE_output,subframe_i,SNR_i);
        end 
    end
   

end
simulation_results.calculate_sim_aggregates(0);
simulation_results.SNR_vector = SNR_vec;


% %% Postprocessing
% 
% if LTE_params.confidence_interval_probability > 0
%     simulation_results = LTE_UL_calculate_confidence_intervals(simulation_results, LTE_params.confidence_interval_probability);
% end
% sim_result_alldata{Bigloop_idx} =  simulation_results;
% 
% simulation_results.SNR_vector = SNR_all(1:Bigloop_idx);
% simulation_results.UE_specific.throughput_coded = [sim_result_temp.UE_specific.throughput_coded, simulation_results.UE_specific.throughput_coded];
% sim_result_temp.UE_specific.throughput_coded = simulation_results.UE_specific.throughput_coded;
% % simulation_results.UE_specific.throughput_uncoded = [sim_result_temp.UE_specific.throughput_uncoded;simulation_results.UE_specific.throughput_uncoded];
% %% Plotting
% 
% if LTE_params.show_plots
%     if isempty(LTE_params.to_plot)
%         LTE_UL_plot_results( simulation_results, LTE_params);
%     else
%         LTE_UL_plot_results( simulation_results, LTE_params, LTE_params.to_plot);
%     end
% end


