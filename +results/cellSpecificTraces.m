classdef cellSpecificTraces < handle
    % Class that stores all of the cell-scpecific traces
    % Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at
    
    properties
        % Fields updated after every TTI
        FER_coded          % Count of incorrect frames every TTI (coded)
        FER_uncoded        % Count of incorrect frames every TTI (uncoded)
        throughput_coded
        throughput_uncoded
        biterrors_coded
        biterrors_uncoded
        blocksize_coded
        blocksize_uncoded
        
        % Aggregates calculated after the simulation is finisched
        BER_coded           % coded BER (each stream)
        BER_uncoded         % uncoded BER (each stream)
        BER_coded_overall   % coded BER (overall)
        BER_uncoded_overall % uncoded BER (overall)
        BLER                % BLER (each stream)
        BLER_overall        % BLER (overall)
        
        % Used to signal what entries are valid
        used_codewords
        
        % Stores the SINR per-subcarrier
        % Note that for AWGN simulations, it is just flat (constant)
        SINR_SC_dB
        
        %confidence
        confidence
        
    end
    
    methods
        
        % Class constructor
        function obj = cellSpecificTraces(N_subframes,SNR_vector_length,maxStreams,nRx,nTx,Ntot,N_users)
            
            obj.biterrors_coded    = zeros(N_subframes, SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
            obj.biterrors_uncoded  = zeros(N_subframes, SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
            obj.blocksize_coded    = zeros(N_subframes, SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
            obj.blocksize_uncoded  = zeros(N_subframes, SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
            obj.FER_coded          = zeros(N_subframes, SNR_vector_length,maxStreams,'uint16'); % 0-65535
            obj.FER_uncoded        = zeros(N_subframes, SNR_vector_length,maxStreams,'uint16'); % 0-65535
            obj.throughput_coded   = zeros(N_subframes, SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
            obj.throughput_uncoded = zeros(N_subframes, SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
            obj.used_codewords     = zeros(N_subframes, SNR_vector_length,maxStreams,'uint16'); % 0-65535 -> What codewords were used
            
            obj.confidence = struct;

        end
        
    end
end
