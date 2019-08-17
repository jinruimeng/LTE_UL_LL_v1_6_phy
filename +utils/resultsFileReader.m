classdef resultsFileReader < handle
    % Encapsulates methods used for reading and plotting data from the results files
    % (c) Josep Colom Ikuno, ITC, 2016
    
    properties
        % The results data
        data
    end
    
    methods
        function obj = resultsFileReader(filename)
            fprintf('Loading results file: %s, ',filename);
            if exist(filename,'file')
                fprintf('exists, ');
                file_data = dir(filename);
                if length(file_data)==1
                    fprintf('%.1f MB, ',file_data.bytes/1024/1024);
                else
                    error('More than one file found. Non-unique filename');
                end
            else
                error('File does not exist');
            end
            obj.data = load(filename,'LTE_params','simulation_results','SNR_vec','N_subframes');
            fprintf('\n');
        end
        function clear(obj)
            obj.data = [];
        end
        function SNR_vec = get_SNR_vector(obj)
            SNR_vec = obj.data.SNR_vec;
        end
        function [cell_throughput_Mbps UE_throughput_Mbps] = get_throughput_over_SNR(obj)
            % Right now works only for SU, but it should be extended
            LTE_params  = obj.data.LTE_params;
            SNR_vector  = obj.data.SNR_vec;
            sim_results = obj.data.simulation_results;
            nUEs        = LTE_params.nUE;
            N_subframes = obj.data.N_subframes;
            
            rx_bits_UE = zeros(length(SNR_vector),nUEs);
            for u_=1:nUEs
                rx_bits_UE(:,u_) = sum(sum(sim_results.UE_specific(u_).throughput_coded,3),1);
            end
            
            UE_throughput_Mbps   = rx_bits_UE ./ N_subframes ./ 1e-3 ./ 1e6;
            cell_throughput_Mbps = sum(rx_bits_UE,2) ./ N_subframes ./ 1e-3 ./ 1e6;
        end
        
        function [UE_BLER] = get_UE_BLER_over_SNR(obj)
            % Right now works only for SU, but it should be extended
            LTE_params  = obj.data.LTE_params;
            SNR_vector  = obj.data.SNR_vec;
            sim_results = obj.data.simulation_results;
            nUEs        = LTE_params.nUE;
            N_subframes = obj.data.N_subframes;
            
            UE_BLER = zeros(length(SNR_vector),nUEs);
            for u_=1:nUEs
                UE_BLER(:,u_) = sim_results.UE_specific(u_).BLER_overall;
            end
        end
        
        function LTE_config = get_simulation_config(obj)
            LTE_config = obj.data.LTE_params;
        end
        function simulation_basic_data = get_config_params_summary(obj)
            LTE_params = obj.data.LTE_params;
            simulation_basic_data.bandwidth = sprintf('%.1f MHz',LTE_params.Bandwidth/1e6);
            switch LTE_params.UE_config.mode
                case 1
                simulation_basic_data.tx_mode = 'SISO';
                case 2
                    simulation_basic_data.tx_mode = 'TxD';
                case 3
                    simulation_basic_data.tx_mode = 'OLSM';
                case 4
                    simulation_basic_data.tx_mode = 'CLSM';
            end
            simulation_basic_data.antenna_config    = sprintf('%.0fx%.0f',LTE_params.UE_config.nTX,LTE_params.BS_config.nRX);
            simulation_basic_data.channel           = LTE_params.ChanMod_config.type;
            simulation_basic_data.simulation_length = obj.data.N_subframes;
            simulation_basic_data.feedback_delay    = LTE_params.downlink_delay;
            simulation_basic_data.scheduler         = LTE_params.scheduler.type;
            simulation_basic_data.eNodeB_count      = LTE_params.nBS;
            simulation_basic_data.UE_count          = LTE_params.nUE;
        end
    end
    methods (Static)
        function delete_results_file_info_except_throughput_and_BLER(filename,new_filename)
            % Loads the results file in filename and deletes all of the
            % results except the ones related to the BLER and throughput
            % (makes the results file smaller
            
            % Load data
            the_data = load(filename);
            
            % Delete unwanted data
            the_data.simulation_results.cell_specific = [];
            
            for u_=1:length(the_data.simulation_results.UE_specific)
                the_data.simulation_results.UE_specific(u_).biterrors_coded       = [];
                the_data.simulation_results.UE_specific(u_).biterrors_uncoded     = [];
                the_data.simulation_results.UE_specific(u_).blocksize_coded       = [];
                the_data.simulation_results.UE_specific(u_).blocksize_uncoded     = [];
                the_data.simulation_results.UE_specific(u_).FER_coded             = [];
                the_data.simulation_results.UE_specific(u_).FER_uncoded           = [];
%                 the_data.simulation_results.UE_specific(u_).DM_chan_est_error     = [];
                the_data.simulation_results.UE_specific(u_).used_codewords        = [];
%                 the_data.simulation_results.UE_specific(u_).used_CQI              = [];
                the_data.simulation_results.UE_specific(u_).ACK                   = [];
                the_data.simulation_results.UE_specific(u_).RBs_assigned          = [];
                the_data.simulation_results.UE_specific(u_).rv_idx                = [];
                the_data.simulation_results.UE_specific(u_).BER_coded             = [];
                the_data.simulation_results.UE_specific(u_).BER_uncoded           = [];
                the_data.simulation_results.UE_specific(u_).BER_coded_overall     = [];
                the_data.simulation_results.UE_specific(u_).BER_uncoded_overall   = [];
                the_data.simulation_results.UE_specific(u_).MSE_freq_offset       = [];
%                 the_data.simulation_results.UE_specific(u_).DM_chan_est_MSE       = [];
            end
            
            % Delete other stuff
             the_data.BS               = [];
             the_data.BS_output        = [];
             the_data.downlinkChannel  = [];
             the_data.mapping_data     = [];
             the_data.UE               = [];
             the_data.ChanMod          = [];
             the_data.UE_output        = [];
             the_data.complete_results = [];
             the_data.res              = [];
             the_data.averager         = [];
             the_data.Wn               = [];
             the_data.out              = [];
             the_data.channel          = [];
            
            % Re-save data
            save(new_filename,'-struct','the_data');
        end
    end
end

