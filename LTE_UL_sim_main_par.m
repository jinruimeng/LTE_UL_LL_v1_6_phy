% LTE system simulator main simulator file. Check the LTE_sim_batch files
% to check how to launch the simulator.
% [] = LTE_sim_main()
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% By using this simulator, you agree to the license terms stated in the license agreement included with this work.
% If you are using the simulator for your scientific work, please reference:


clear tmp_results;
clear filename_suffix;
clear StartTime;
clear maxStreams;
clear maxi;
clear num;
% clear output_filename;
clear simulation_results;

maxStreams = min(LTE_params.UE_config.nTX,LTE_params.BS_config.nRX);
simulation_results = results.simulationResultsUL( LTE_params );

maxi = max([UE.nTX]);
Maximum_number_of_packets=500;

% temporary result variables for parfor
UE_res =  struct(... 
    'ACK',false(N_subframes,maxStreams),...
    'rv_idx',zeros(N_subframes,maxStreams,'uint8'),...
    'RBs_assigned',zeros(N_subframes,'uint8'),...
    'biterrors_coded',zeros(N_subframes,maxStreams,'uint32'),...
    'biterrors_uncoded',zeros(N_subframes,maxStreams,'uint32'),...
    'blocksize_coded',zeros(N_subframes,maxStreams,'uint32'),...
    'blocksize_uncoded',zeros(N_subframes,maxStreams,'uint32'),...
    'throughput_coded',zeros(N_subframes,maxStreams,'uint32'),...
    'throughput_uncoded',zeros(N_subframes,maxStreams,'uint32'),...
    'FER_coded',zeros(N_subframes,maxStreams,'uint16'),...
    'FER_uncoded',zeros(N_subframes,maxStreams,'uint16'),...
    'used_codewords',zeros(N_subframes,maxStreams,'uint16'),...
    'channel_error',zeros(N_subframes,1,'double'),...
    'channel_pred_error',zeros(N_subframes,1,'double'),...
    'papr',zeros(LTE_params.Nsub,N_subframes,'double'),...
    'type',cell(1),...
    'data_buffer_left',zeros(1,Maximum_number_of_packets,'uint32'),...
    'data_generated',zeros(1,Maximum_number_of_packets,'uint32'),...
    'TTI_origin',zeros(1,Maximum_number_of_packets,'uint32'),...
    'delay_TTI',zeros(1,Maximum_number_of_packets,'uint32'),...
    'ID_count_current',zeros(1,'uint32'),...
    'ID_count_next',zeros(1,'uint32'), ...
    'used_cqi',zeros(N_subframes,maxStreams,'uint32'), ...
    'used_RI',zeros(N_subframes,1, 'uint32'));

cell_res =  struct(...
    'biterrors_coded',zeros(N_subframes,maxStreams,'uint32'),...
    'biterrors_uncoded',zeros(N_subframes,maxStreams,'uint32'),...
    'blocksize_coded',zeros(N_subframes,maxStreams,'uint32'),...
    'blocksize_uncoded',zeros(N_subframes,maxStreams,'uint32'),...
    'throughput_coded',zeros(N_subframes,maxStreams,'uint32'),...
    'throughput_uncoded',zeros(N_subframes,maxStreams,'uint32'),...
    'FER_coded',zeros(N_subframes,maxStreams,'uint16'),...
    'FER_uncoded',zeros(N_subframes,maxStreams,'uint16'),...
    'used_codewords',zeros(N_subframes,maxStreams,'uint16'),...
    'SINR_SC_dB',zeros(LTE_params.Ntot,N_subframes));
     
tmp_results.UE_specific = repmat(UE_res,1,LTE_params.nUE*LTE_params.nBS);
tmp_results.cell_specific = repmat(cell_res, LTE_params.nBS);
tmp_results = repmat(tmp_results,1,size(SNR_vec,2));
tmp_channel = channel;
clear UE_res cell_res;

%% Parallel toolbox

if verLessThan('matlab','8.3')
    % for MATLAB older than 2014a (8.3)
    if ~matlabpool('size')
        matlabpool('open');
    end
    
else
    % for MATLAB 2014a (8.3) or newer

    %delete(gcp);
    %num = parpool('local');
    p = gcp('nocreate');

    if isempty(p) % from matlab doc
        poolsize = 0;
    else
        poolsize = p.NumWorkers;
    end

    if strcmp(LTE_params.simulation_method,'parallel') && poolsize == 0;
        parpool('local');
    end
end

%% SNR and frame loops
StartTime = clock;
% NOTE: changing this with a parfor does not change any functionality, but
% messes up the timer (it goes backwards), so if you really want parallel
% functionality, change the line below to a parfor

%% File to save the progress through the SNR loop (necessary for estimation of remaining time - works only if all workers access the same harddisk)
fid = fopen('counter.m','w');
fprintf(fid,'%c','');
fclose(fid);

nUE = LTE_params.nUE;
Ns = LTE_params.Ns;

parfor SNR_i=1:size(SNR_vec,2) % parfor
    UE_tmp = UE;
    BS_tmp = BS;
    out_tmp = out;
    LTE_params_tmp = LTE_params;
    task = getCurrentTask;
    task_ID = get(task,'ID');
    
    % temporary result variables for parfor
    tmp_UE_ACK                  = false(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS);
    tmp_UE_rv_idx               = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint8');
    tmp_UE_RBs_assigned         = zeros(N_subframes,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint8');
    tmp_UE_biterrors_coded      = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_biterrors_uncoded    = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_blocksize_coded      = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_blocksize_uncoded    = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_throughput_coded     = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_throughput_uncoded   = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_FER_coded            = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint16');
    tmp_UE_FER_uncoded          = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint16');
    tmp_UE_used_codewords       = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint16');
    tmp_UE_used_cqi             = zeros(N_subframes,maxStreams,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_used_RI              = zeros(N_subframes,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_channel_error        = zeros(N_subframes,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'double');
    tmp_UE_channel_pred_error   = zeros(N_subframes,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'double');
    tmp_UE_papr                 = zeros(N_subframes,LTE_params_tmp.Nsub,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'double');
    
    tmp_UE_type                 = cell(1,LTE_params_tmp.nUE*LTE_params_tmp.nBS);
    tmp_UE_TTI_origin           = zeros(Maximum_number_of_packets,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_delay_TTI            = zeros(Maximum_number_of_packets,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_ID_count_current     = zeros(1,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    tmp_UE_ID_count_next        = zeros(1,LTE_params_tmp.nUE*LTE_params_tmp.nBS,'uint32');
    
    tmp_cell_biterrors_coded    = zeros(N_subframes,maxStreams, LTE_params_tmp.nBS, 'uint32');
    tmp_cell_biterrors_uncoded  = zeros(N_subframes,maxStreams, LTE_params_tmp.nBS, 'uint32');
    tmp_cell_blocksize_coded    = zeros(N_subframes,maxStreams, LTE_params_tmp.nBS, 'uint32');
    tmp_cell_blocksize_uncoded  = zeros(N_subframes,maxStreams, LTE_params_tmp.nBS, 'uint32');
    tmp_cell_throughput_coded   = zeros(N_subframes,maxStreams, LTE_params_tmp.nBS, 'uint32');
    tmp_cell_throughput_uncoded = zeros(N_subframes,maxStreams, LTE_params_tmp.nBS, 'uint32');
    tmp_cell_FER_coded          = zeros(N_subframes,maxStreams, LTE_params_tmp.nBS, 'uint16');
    tmp_cell_FER_uncoded        = zeros(N_subframes,maxStreams, LTE_params_tmp.nBS, 'uint16');
    tmp_cell_used_codewords     = zeros(N_subframes,maxStreams, LTE_params_tmp.nBS, 'uint16');
    tmp_SINR_SC_dB              = zeros(LTE_params_tmp.Ntot,N_subframes,  LTE_params_tmp.nBS);

    % Initialize variables that will be reused, such as the BS_output
    BS_output = outputs.bsOutput(LTE_params.nUE*LTE_params_tmp.nBS,LTE_params.Nrb,maxStreams);
    subframe_i = 1;
    delay_counter = 0;   
    
    delay = true;
    % Set network clock. It will tell every network element in which TTI we are right now.
    network_clock = network_elements.clock(LTE_params.FrameDur/10);

    
    % PROKOPEC 2011_11_23
    % Initialize uplink channel
    downlinkChannel = channels.downlinkChannel(LTE_params.downlink_delay,LTE_params.nUE,LTE_params.nBS);
    % Receive feedback from the previous subframe
    UE_output = downlinkChannel.receive_feedback;
    % Attach the clock to each network element
    % Attach UE_output.genie
    for u_ = 1:(LTE_params.nUE*LTE_params.nBS)
        UE_tmp(u_).clock = network_clock;
        UE_output(u_).UE_genie = outputs.genieInformation;
    end

    for b_=1:LTE_params.nBS
        BS_tmp(b_).clock = network_clock;
        BS_tmp(b_).reset_HARQ_process_index;
        %reset number of channel realizations
        BS_tmp(b_).realization_num = 0;
    end
    
    % Some further initialization
    SNR = SNR_vec(:,SNR_i);   
    network_clock.reset;
    ChanMod_output = cell(LTE_params.nBS,LTE_params.nUE);
    
    
    % ACK of the previous frame. If this is the first frame, set the
    % ACK to correct so that the HARQ handling generates new data
    if subframe_i==1 || (subframe_i <= LTE_params.downlink_delay && ~strcmp(LTE_params.ChanMod_config.time_correlation,'independent'))% (number of max HARQ processes that will be used)
        for uu = 1:LTE_params.nUE
            UE_output(uu).ACK    = true(1,2);
            UE_output(uu).rv_idx = zeros(1,2);
        end
    end
    
    if DEBUG_LEVEL > 0
        disp('');
        %disp(['*************** SNR = ' num2str(SNR(1)) 'dB, value ' num2str(SNR_i) ' of ' num2str(size(SNR_vec,2)) ' ***************']);
        fprintf('*************** SNR = %0 6.2f db, value %d of %d  ***************\n', SNR(1), SNR_i, size(SNR_vec,2));
    end
    
    while subframe_i <= N_subframes  
        % First of all, advance the network clock
        network_clock.advance_1_TTI;
        if mod(subframe_i,50)==1 && ~delay
            if DEBUG_LEVEL > 0
                if task_ID == 1 
                    job = getCurrentJob;
                    task_nr = length(get(job,'Tasks'));
                    fid = fopen('counter.m','r+');
                    A = fscanf(fid,'%d');
                    done = numel(A);
                    fclose(fid);
                   if (task_nr > size(SNR_vec,2) - done)
                        task_nr = size(SNR_vec,2) - done;
                   end
                    disp(['   processing subframe #' num2str(subframe_i) ' of ' num2str(N_subframes)])
                    disp(['---> remaining simulation time: ' num2str(etime(clock,StartTime)/60*((N_subframes*size(SNR_vec,2))-(subframe_i-1)*task_nr-N_subframes*done)/(N_subframes*done+(subframe_i-1)*task_nr),'%5.3f') 'min']);
                end
            end
        end
        
        [BS_output, UE_output] = LTE_UL_sim_main_subframe(LTE_params_tmp, ChanMod, SNR, UE_tmp, UE_output, BS_tmp, BS_output, subframe_i, out_tmp, downlinkChannel, cqi_i, tmp_channel);
        
        if delay_counter == LTE_params_tmp.downlink_delay;
            delay = false;
        end
        delay_counter = delay_counter + 1;
                       
        % Process results for this TTI
        % simulation_results.process_TTI_results(BS_output,UE_output,subframe_i,SNR_i);
%         if ~delay
        if true    
            [tmp_UE_ACK(subframe_i,:,:),...
             tmp_UE_rv_idx(subframe_i,:,:),...
             tmp_UE_RBs_assigned(subframe_i,:),...
             tmp_UE_biterrors_coded(subframe_i,:,:),...
             tmp_UE_biterrors_uncoded(subframe_i,:,:),...
             tmp_UE_blocksize_coded(subframe_i,:,:),...
             tmp_UE_blocksize_uncoded(subframe_i,:,:),...
             tmp_UE_throughput_coded(subframe_i,:,:),...
             tmp_UE_FER_coded(subframe_i,:,:),...
             tmp_UE_throughput_uncoded(subframe_i,:,:),...
             tmp_UE_FER_uncoded(subframe_i,:,:),...
             tmp_UE_used_codewords(subframe_i,:,:),...
             tmp_UE_channel_error(subframe_i,:),...
             tmp_UE_channel_pred_error(subframe_i,:),...
             tmp_UE_papr(subframe_i,:,:),...
             tmp_cell_biterrors_coded(subframe_i,:, :),...
             tmp_cell_biterrors_uncoded(subframe_i,:, :),...
             tmp_cell_blocksize_coded(subframe_i,:, :),...
             tmp_cell_blocksize_uncoded(subframe_i,:, :),...
             tmp_cell_throughput_coded(subframe_i,:, :),...
             tmp_cell_throughput_uncoded(subframe_i,:, :),...
             tmp_cell_FER_coded(subframe_i,:, :),...
             tmp_cell_FER_uncoded(subframe_i,:, :),...
             tmp_cell_used_codewords(subframe_i,:, :), ...
             tmp_UE_used_cqi(subframe_i, :, :), ...
             tmp_UE_used_RI(subframe_i, :)] = LTE_UL_processTTI_results(BS_output,...
                                                                                UE_output,...
                                                                                subframe_i,...
                                                                                SNR_i,...
                                                                                LTE_params_tmp.nUE,...
                                                                                LTE_params_tmp.nBS, ...
                                                                                tmp_cell_biterrors_coded(subframe_i,:, :),...
                                                                                tmp_cell_biterrors_uncoded(subframe_i,:, :),...
                                                                                tmp_cell_blocksize_coded(subframe_i,:, :),...
                                                                                tmp_cell_blocksize_uncoded(subframe_i,:, :),...
                                                                                tmp_cell_throughput_coded(subframe_i,:, :),...
                                                                                tmp_cell_throughput_uncoded(subframe_i,:, :),...
                                                                                tmp_cell_FER_coded(subframe_i,:, :),...
                                                                                tmp_cell_FER_uncoded(subframe_i,:, :),...
                                                                                tmp_cell_used_codewords(subframe_i,:, :),...
                                                                                maxStreams,...
                                                                                LTE_params_tmp,...
                                                                                LTE_params.connection_table);
        end
        subframe_i = subframe_i+1;
                    
    end
    
    
%     tmp_UE_data_buffer_left,
%     tmp_UE_data_generated,
%     tmp_UE_TTI_origin=UE_output.TTI_origin
%     tmp_UE_delay_TTI=UE_output.delay_TTI;
%     tmp_UE_ID_count_current,




    % Set a single output variable
    tmp_results(SNR_i) = LTE_set_results(tmp_UE_ACK,...
                                         tmp_UE_rv_idx,...
                                         tmp_UE_RBs_assigned,...
                                         tmp_UE_biterrors_coded,...
                                         tmp_UE_biterrors_uncoded,...
                                         tmp_UE_blocksize_coded,...
                                         tmp_UE_blocksize_uncoded,...
                                         tmp_UE_throughput_coded,...
                                         tmp_UE_throughput_uncoded,...
                                         tmp_UE_FER_coded,...
                                         tmp_UE_FER_uncoded,...
                                         tmp_UE_used_codewords,...
                                         tmp_UE_channel_error,...
                                         tmp_UE_channel_pred_error,...
                                         tmp_UE_papr,...
                                         tmp_UE_type,...
                                         tmp_UE_TTI_origin,...
                                         tmp_UE_delay_TTI,...
                                         tmp_UE_ID_count_current,...
                                         tmp_UE_ID_count_next,...
                                         tmp_cell_biterrors_coded,...
                                         tmp_cell_biterrors_uncoded,...
                                         tmp_cell_blocksize_coded,...
                                         tmp_cell_blocksize_uncoded,...
                                         tmp_cell_throughput_coded,...
                                         tmp_cell_throughput_uncoded,...
                                         tmp_cell_FER_coded,...
                                         tmp_cell_FER_uncoded,...
                                         tmp_cell_used_codewords,...
                                         LTE_params_tmp.nUE,...
                                         LTE_params_tmp.nBS, ...
                                         tmp_SINR_SC_dB, ...
                                         tmp_UE_used_cqi, ...
                                         tmp_UE_used_RI);
    
    fid = fopen('counter.m','a');
    fprintf(fid,'%d\n',SNR_i);
    fclose(fid);
end
simulation_results.SNR_vector = SNR_vec;
simulation_results.set_TTI_results(tmp_results);
clear tmp_results;
clear U_temp;
clear BS_temp;
clear LTE_params_tmp;
clear out_tmp;



%% Calculate simulation aggregates
if(strcmp(BS(1).autocorrelation_matrix_type,'estimated') && strcmp(BS(1).channel_estimation_method,'MMSE'))
    simulation_results.calculate_sim_aggregates(UE(1).realization_num_total/(BS.nRX * ChanMod.nTX));
else
    simulation_results.calculate_sim_aggregates(0);
end

%% Postprocessing

if LTE_params.confidence_interval_probability > 0
    simulation_results = LTE_UL_calculate_confidence_intervals(simulation_results, LTE_params.confidence_interval_probability);
end

%% Show plots at the end of the simulation
if LTE_params.show_plots
    if isempty(LTE_params.to_plot)
        LTE_UL_plot_results( simulation_results, LTE_params);
    else
        LTE_UL_plot_results( simulation_results, LTE_params, LTE_params.to_plot);
    end
end
