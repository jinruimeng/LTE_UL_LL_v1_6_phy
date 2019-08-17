classdef ueSpecificTraces < handle
% Class that stores all of the UE-scpecific traces
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

properties
    % Fields updated after every TTI
    ACK
    RBs_assigned
    rv_idx
    biterrors_coded
    biterrors_uncoded
    blocksize_coded
    blocksize_uncoded
    FER_coded
    FER_uncoded
    throughput_coded
    throughput_uncoded
    BER_PCFICH
    BER_PCFICH_CFI
    BER_PHICH
    papr
    
    % Aggregates calculated after the simulation is finisched
    BER_coded           % coded BER (each stream)
    BER_uncoded         % uncoded BER (each stream)
    BER_coded_overall   % coded BER (overall)
    BER_uncoded_overall % uncoded BER (overall)
    BLER                % BLER (each stream)
    BLER_overall        % BLER (overall)
    channel_error       % channel estimation error for each UE
    channel_pred_error  % channel prediction error for each UE
    used_codewords      % Used to signal what entries are valid
    
    used_cqi
    used_RI
    
    % confidence
    confidence
    
    %Traffic
    type
    data_buffer_left
    data_generated
    TTI_origin
    delay_TTI
    ID_count_current
    ID_count_next
end

   methods
       % Class contructor. Data preallocation
       function obj = ueSpecificTraces(N_subframes,SNR_vector_length,maxStreams,nRx,nTx)
           obj.BER_PCFICH               = zeros(N_subframes,SNR_vector_length);
           obj.BER_PCFICH_CFI           = zeros(N_subframes,SNR_vector_length);
           obj.BER_PHICH                = zeros(N_subframes,SNR_vector_length);
           obj.ACK                      = false(N_subframes,SNR_vector_length,maxStreams);          % true/false -> ACK of the received subframes for BLER calculation
           obj.RBs_assigned             = zeros(N_subframes,SNR_vector_length,'uint8');  % 0-255      -> number of assigned RBs in every TTI
           obj.rv_idx                   = zeros(N_subframes,SNR_vector_length,maxStreams,'uint8');  % 0-255      -> redundancy version index of the received subframes
           obj.biterrors_coded          = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.biterrors_uncoded        = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.blocksize_coded          = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.blocksize_uncoded        = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.FER_coded                = false(N_subframes,SNR_vector_length,maxStreams);          % This is actually ~ACK
           obj.FER_uncoded              = false(N_subframes,SNR_vector_length,maxStreams);          %
           obj.throughput_coded         = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.throughput_uncoded       = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.used_codewords           = false(N_subframes,SNR_vector_length,maxStreams);          % What codewords were used
           obj.channel_error            = zeros(N_subframes,SNR_vector_length);
           obj.channel_pred_error       = zeros(N_subframes,SNR_vector_length);
           obj.used_cqi                 = zeros(N_subframes,SNR_vector_length,maxStreams);
           obj.used_RI                 = zeros(N_subframes,SNR_vector_length);
           
%            obj.data_traffic             = traffic_models.generic_tm;
           %traffic
           % % maximum number of packet to be simulated
           Maximum_number_of_packets=500;
           obj.type                 = cell(1,SNR_vector_length);
           obj.data_buffer_left     = zeros(Maximum_number_of_packets,SNR_vector_length,'uint32');
           obj.data_generated       = zeros(Maximum_number_of_packets,SNR_vector_length,'uint32');
           obj.TTI_origin           = zeros(Maximum_number_of_packets,SNR_vector_length,'uint32');
           obj.delay_TTI            = zeros(Maximum_number_of_packets,SNR_vector_length,'uint32');
           obj.ID_count_current     = zeros(1,SNR_vector_length,'uint32');
           obj.ID_count_next        = zeros(1,SNR_vector_length,'uint32');
           
           obj.confidence = struct;
       end
   end
end 
