classdef simulationResultsUL < handle
    % Class that stores all of the simulation results. Both Cell specific and
    % UE specific.
    % Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at

    properties
        cell_specific  % Cell specific traces
        UE_specific    % UE specific traces
        nUE            % number of users per basestation
        nBS            % number of basestations 
        N_subframes
        maxStreams
        SNR_vector     % SNR vector used for this simulation
        connection_table
        nRx
        nTx
        Ntot
        conf_interval_probability % (1-alpha) in bootci
        
    end

    methods
        % Class constructor. Preallocate all necessary space
        % BS_nRx num of eNodeB receive antennas, Prokopec
        function obj = simulationResultsUL( LTE_params )
            % initialize
            obj.nUE = LTE_params.nUE; 
            obj.nBS = LTE_params.nBS;
            obj.N_subframes = LTE_params.N_subframes;
            obj.maxStreams  = min(LTE_params.UE_config.nTX,LTE_params.BS_config.nRX);
            obj.SNR_vector  = LTE_params.SNR_vec;
            obj.connection_table = LTE_params.connection_table;
            obj.nTx = LTE_params.UE_config.nTX;
            obj.nRx = LTE_params.BS_config.nRX;
            obj.Ntot = LTE_params.Ntot;
            obj.conf_interval_probability = LTE_params.confidence_interval_probability;
            
            SNR_vector_length = size(obj.SNR_vector,2);

            % Preallocate cell specific traces
            obj.cell_specific = results.cellSpecificTraces(1,1,1,1,1,obj.Ntot,obj.nUE);
            for b_ = 1:obj.nBS
                obj.cell_specific(b_) = results.cellSpecificTraces(obj.N_subframes,SNR_vector_length,obj.maxStreams,obj.nRx,obj.nTx,obj.Ntot,obj.nUE);
            end

            % Preallocate UE-specific traces
            obj.UE_specific = results.ueSpecificTraces(1,1,1,1,1);
            for u_= 1:(obj.nUE * obj.nBS)
                obj.UE_specific(u_) = results.ueSpecificTraces(obj.N_subframes,SNR_vector_length,obj.maxStreams,obj.nRx,obj.nTx);
            end
            
        end

        function process_TTI_results(obj, BS_output, UE_output, subframe_i,SNR_i)
            % Loop over all UEs and BS! and streams
            
            for uu = 1:obj.nUE * obj.nBS
                
                bb = utils.findBS(obj.connection_table, uu);
                local_uu = utils.globalToLocalUser( obj.connection_table, bb, uu ); 
                
                
                for stream_i = 1:BS_output.UE_signaling_UL(bb, local_uu).MCS_and_scheduling_UL.nCodewords
                    
                    if BS_output.UE_signaling_UL(bb, local_uu).MCS_and_scheduling_UL.assigned_RBs
                        % Update ACK (UE traces) and biterrors (cell traces)
                        obj.UE_specific(uu).ACK(subframe_i,SNR_i,stream_i)             = UE_output(uu).ACK(stream_i);
                        obj.UE_specific(uu).rv_idx(subframe_i,SNR_i,stream_i)          = UE_output(uu).rv_idx(stream_i);
                        obj.UE_specific(uu).RBs_assigned(subframe_i,SNR_i)             = BS_output.UE_signaling_UL(bb, local_uu).MCS_and_scheduling_UL.assigned_RBs;
                        obj.UE_specific(uu).used_cqi(subframe_i, SNR_i, stream_i)      = BS_output.UE_signaling_UL(bb, local_uu).MCS_and_scheduling_UL.cqi(stream_i);
                        obj.UE_specific(uu).used_RI(subframe_i, SNR_i)                 = BS_output.UE_signaling_UL(bb, local_uu).MCS_and_scheduling_UL.nLayers;
                    else
                         % Update ACK (UE traces) and biterrors (cell traces)
                        obj.UE_specific(uu).ACK(subframe_i,SNR_i,stream_i)             = UE_output(uu).ACK(stream_i);
                        if subframe_i ~= 1
                            obj.UE_specific(uu).rv_idx(subframe_i,SNR_i,stream_i)       = obj.UE_specific(uu).rv_idx(subframe_i-1,SNR_i,stream_i);
                        else
                            obj.UE_specific(uu).rv_idx(subframe_i,SNR_i,stream_i)       = 0;
                        end
                        obj.UE_specific(uu).RBs_assigned(subframe_i,SNR_i)             = BS_output.UE_signaling_UL(bb, local_uu).MCS_and_scheduling_UL.assigned_RBs;
                    end
                    
                    % Update the stats that are only meaningful when the UE has been scheduled
                    if UE_output(uu).UE_scheduled

                        % Update biterrors (UE traces)
                        if ~isempty(UE_output(uu).UE_genie.data_bits{stream_i}) 
                            obj.UE_specific(uu).biterrors_coded(subframe_i,SNR_i,stream_i)   = sum(abs(UE_output(uu).rx_data_bits{stream_i}  - UE_output(uu).UE_genie.data_bits{stream_i}));
                        else
                            obj.UE_specific(uu).biterrors_coded(subframe_i,SNR_i,stream_i)   = 0;
                        end
                        %obj.UE_specific(uu).biterrors_uncoded(subframe_i,SNR_i,stream_i) = sum(xor(double(UE_output(uu).rx_coded_bits{stream_i}), double(BS_output.genie(uu).sent_bits{stream_i})));
                        obj.UE_specific(uu).biterrors_uncoded(subframe_i,SNR_i,stream_i) = sum(abs(UE_output(uu).rx_coded_bits{stream_i} - BS_output.genie(bb,local_uu).sent_bits{stream_i}));
                        
                        % Update blocksize (UE traces)
                        obj.UE_specific(uu).blocksize_coded(subframe_i,SNR_i,stream_i)   = length(UE_output(uu).rx_data_bits{stream_i});
                        obj.UE_specific(uu).blocksize_uncoded(subframe_i,SNR_i,stream_i) = length(UE_output(uu).rx_coded_bits{stream_i});
                        
                        % Update PAPR
                        obj.UE_specific(uu).papr(subframe_i,SNR_i,:) = UE_output(uu).papr;
                        
                        % Update FER and throughput (UE traces)
                        
                        % Coded
                        if UE_output(uu).ACK(stream_i)
                            obj.UE_specific(uu).throughput_coded(subframe_i,SNR_i,stream_i) = length(UE_output(uu).rx_data_bits{stream_i});
                            obj.UE_specific(uu).FER_coded(subframe_i,SNR_i,stream_i) = 0;
                        else
                            obj.UE_specific(uu).FER_coded(subframe_i,SNR_i,stream_i) = 1;
                        end
                        

                        if (obj.UE_specific(uu).biterrors_uncoded(subframe_i,SNR_i,stream_i)==0)
                            obj.UE_specific(uu).throughput_uncoded(subframe_i,SNR_i,stream_i) = length(UE_output(uu).rx_coded_bits{stream_i});
                            obj.UE_specific(uu).FER_uncoded(subframe_i,SNR_i,stream_i)        = 0;
                        else
                            obj.UE_specific(uu).FER_uncoded(subframe_i,SNR_i,stream_i) = 1;
                        end
                        
                        % Update what codewords were used (valid positions in the traces, 0 if no RBs were allocated)
                        obj.UE_specific(uu).used_codewords(subframe_i,SNR_i,stream_i) = 1;
                        
                        %obj.UE_specific(uu)
                        
                    else
                        % Update what codewords were used (valid positions in the traces, 0 if no RBs were allocated)
                        obj.UE_specific(uu).used_codewords(subframe_i,SNR_i,stream_i) = 0;
                    end
                    
                end
                
                % update channel estimation error
                if UE_output(uu).UE_scheduled
                    % channel error is only nonzero if the UE was scheduled --> take the mean only of scheduled TTIs
                    obj.UE_specific(uu).channel_error(subframe_i,SNR_i)         = obj.UE_specific(uu).channel_error(subframe_i,SNR_i)           + UE_output(uu).channel_estimation_error;
                    obj.UE_specific(uu).channel_pred_error(subframe_i,SNR_i)    = obj.UE_specific(uu).channel_pred_error(subframe_i,SNR_i)      + UE_output(uu).channel_prediction_error;
                end

                % Update cell coded and uncoded bit errors
                
                if utils.isUEattached( obj.connection_table, bb, uu )
                    obj.cell_specific(bb).biterrors_coded(subframe_i,SNR_i,:)   = obj.cell_specific(bb).biterrors_coded(subframe_i,SNR_i,:)     + obj.UE_specific(uu).biterrors_coded(subframe_i,SNR_i,:);
                    obj.cell_specific(bb).biterrors_uncoded(subframe_i,SNR_i,:) = obj.cell_specific(bb).biterrors_uncoded(subframe_i,SNR_i,:)   + obj.UE_specific(uu).biterrors_uncoded(subframe_i,SNR_i,:);

                    % Update blocksize (cell traces)
                    obj.cell_specific(bb).blocksize_coded(subframe_i,SNR_i,:)   = obj.cell_specific(bb).blocksize_coded(subframe_i,SNR_i,:)     + obj.UE_specific(uu).blocksize_coded(subframe_i,SNR_i,:);
                    obj.cell_specific(bb).blocksize_uncoded(subframe_i,SNR_i,:) = obj.cell_specific(bb).blocksize_uncoded(subframe_i,SNR_i,:)   + obj.UE_specific(uu).blocksize_uncoded(subframe_i,SNR_i,:);

                    % Update FER and throughput (cell traces)
                    obj.cell_specific(bb).throughput_coded(subframe_i,SNR_i,:)  = obj.cell_specific(bb).throughput_coded(subframe_i,SNR_i,:)    + obj.UE_specific(uu).throughput_coded(subframe_i,SNR_i,:);
                    obj.cell_specific(bb).throughput_uncoded(subframe_i,SNR_i,:)= obj.cell_specific(bb).throughput_uncoded(subframe_i,SNR_i,:)  + obj.UE_specific(uu).throughput_uncoded(subframe_i,SNR_i,:);
                    obj.cell_specific(bb).FER_coded(subframe_i,SNR_i,:)         = obj.cell_specific(bb).FER_coded(subframe_i,SNR_i,:)           + uint16(obj.UE_specific(uu).FER_coded(subframe_i,SNR_i,:));
                    obj.cell_specific(bb).FER_uncoded(subframe_i,SNR_i,:)       = obj.cell_specific(bb).FER_uncoded(subframe_i,SNR_i,:)         + uint16(obj.UE_specific(uu).FER_uncoded(subframe_i,SNR_i,:));

                    % The number of codewords is the maximum (bitwise OR) of the used codewords matrix
                    obj.cell_specific(bb).used_codewords(subframe_i,SNR_i,:) = obj.cell_specific(bb).used_codewords(subframe_i,SNR_i,:) + uint16(obj.UE_specific(uu).used_codewords(subframe_i,SNR_i,:));
                end               
            end
        end
        
        function save_traffic_result(obj,UE_output,SNR_i)
            for u_=1:obj.nUE*obj.nBS
                obj.UE_specific(u_).type{SNR_i} = UE_output(u_).type;
                obj.UE_specific(u_).data_buffer_left(:,SNR_i) = UE_output(u_).data_buffer_left;
                obj.UE_specific(u_).data_generated(:,SNR_i) = UE_output(u_).data_generated;
                obj.UE_specific(u_).TTI_origin(:,SNR_i) = UE_output(u_).TTI_origin;
                obj.UE_specific(u_).delay_TTI(:,SNR_i) = UE_output(u_).delay_TTI;
                obj.UE_specific(u_).ID_count_current(SNR_i) = UE_output(u_).ID_count_current;
                obj.UE_specific(u_).ID_count_next(SNR_i) = UE_output(u_).ID_count_next;
            end          
        end

        % Calculate simulations aggregates.
        function calculate_sim_aggregates(obj,elements_to_remove)
            %remove results, which havent used estimated autocorrelations matrix
            for bb = 1:obj.nBS
                obj.cell_specific(bb).biterrors_coded(1:elements_to_remove,:,:)    = [];
                obj.cell_specific(bb).biterrors_uncoded(1:elements_to_remove,:,:)  = [];
                obj.cell_specific(bb).FER_coded(1:elements_to_remove,:,:)          = [];
                obj.cell_specific(bb).blocksize_coded(1:elements_to_remove,:,:)    = [];
                obj.cell_specific(bb).blocksize_uncoded(1:elements_to_remove,:,:)  = [];
                obj.cell_specific(bb).used_codewords(1:elements_to_remove,:,:)     = [];
                obj.cell_specific(bb).throughput_coded(1:elements_to_remove,:,:)   = [];
                obj.cell_specific(bb).throughput_uncoded(1:elements_to_remove,:,:) = [];
                obj.cell_specific(bb).FER_uncoded(1:elements_to_remove,:,:)        = [];

                % Cell specific results
                obj.cell_specific(bb).BER_coded   = squeeze(sum(obj.cell_specific(bb).biterrors_coded,1)   ./ sum(obj.cell_specific(bb).blocksize_coded,1));
                obj.cell_specific(bb).BER_uncoded = squeeze(sum(obj.cell_specific(bb).biterrors_uncoded,1) ./ sum(obj.cell_specific(bb).blocksize_uncoded,1));
                obj.cell_specific(bb).BLER        = squeeze(sum(obj.cell_specific(bb).FER_coded,1) ./ sum(obj.cell_specific(bb).used_codewords,1));

                obj.cell_specific(bb).BER_coded_overall   = squeeze( sum(sum(obj.cell_specific(bb).biterrors_coded,1),3)   ./ sum(sum(obj.cell_specific(bb).blocksize_coded,1),3) );
                obj.cell_specific(bb).BER_uncoded_overall = squeeze( sum(sum(obj.cell_specific(bb).biterrors_uncoded,1),3) ./ sum(sum(obj.cell_specific(bb).blocksize_uncoded,1),3) );
                obj.cell_specific(bb).BLER_overall        = squeeze( sum(sum(obj.cell_specific(bb).FER_coded,1),3) ./ sum(sum(obj.cell_specific(bb).used_codewords,1),3));

            
            end
            %obj.cell_specific.confidence.ber_uncoded =[mean(ber_uncoded_cell); bootci(2000, {@mean, ber_uncoded_cell})];
            
            
            % UE-specific results
            for uu = 1:obj.nUE*obj.nBS
                obj.UE_specific(uu).biterrors_coded(1:elements_to_remove,:,:)   = [];
                obj.UE_specific(uu).biterrors_uncoded(1:elements_to_remove,:,:) = [];
                obj.UE_specific(uu).FER_coded(1:elements_to_remove,:,:)         = [];
                obj.UE_specific(uu).blocksize_coded(1:elements_to_remove,:,:)   = [];
                obj.UE_specific(uu).blocksize_uncoded(1:elements_to_remove,:,:) = [];
                obj.UE_specific(uu).used_codewords(1:elements_to_remove,:,:)    = [];
                
                obj.UE_specific(uu).BER_coded    = sum(obj.UE_specific(uu).biterrors_coded,1)   ./ sum(obj.UE_specific(uu).blocksize_coded,1);
                obj.UE_specific(uu).BER_uncoded  = sum(obj.UE_specific(uu).biterrors_uncoded,1) ./ sum(obj.UE_specific(uu).blocksize_uncoded,1);
                obj.UE_specific(uu).BLER         = squeeze(sum(obj.UE_specific(uu).FER_coded,1) ./ sum(obj.UE_specific(uu).used_codewords,1));
                
                obj.UE_specific(uu).BER_coded_overall   = squeeze(sum(sum(obj.UE_specific(uu).biterrors_coded,1),3)   ./ sum(sum(obj.UE_specific(uu).blocksize_coded,1),3));
                obj.UE_specific(uu).BER_uncoded_overall = squeeze(sum(sum(obj.UE_specific(uu).biterrors_uncoded,1),3) ./ sum(sum(obj.UE_specific(uu).blocksize_uncoded,1),3));
                obj.UE_specific(uu).BLER_overall        = squeeze(sum(sum(obj.UE_specific(uu).FER_coded,1),3) ./ sum(sum(obj.UE_specific(uu).used_codewords,1),3));
                
                
                
                % confidence intervals
                % (row 1: means, row 2: lower bound, row 3: upper bound) x SNR length 
                % consider only channel MSE when UE was actually scheduled
                for snr_i = 1:size(obj.SNR_vector,2)                                                                % number of SNR points
                    channel_error_tmp = obj.UE_specific(uu).channel_error(:,snr_i);                                 % channel error, for one UE and one SNR
                    channel_error_tmp_scheduled = channel_error_tmp(obj.UE_specific(uu).RBs_assigned(:,snr_i)>0);   % scheduled UE's channel estimation error
                    TTIs_scheduled = sum((obj.UE_specific(uu).RBs_assigned(:,snr_i)>0),1);                          % number of subframes the UE was scheduled                    
                    ber_func = @(bit_error,nr_of_bits) sum(bit_error)/sum(nr_of_bits);                              % function handle for bootci
                    if size(channel_error_tmp_scheduled,1) > 2
                        obj.UE_specific(uu).confidence.channel_error(:,snr_i) = [mean(channel_error_tmp_scheduled); bootci(2000,{ber_func, channel_error_tmp_scheduled,TTIs_scheduled},'alpha',0.05)];
                    else
                        obj.UE_specific(uu).confidence.channel_error(:,snr_i) = NaN(3,1);
                    end
                end
                %obj.UE_specific(uu).channel_error = obj.UE_specific(uu).confidence.channel_error(1,:);              % now overwrite the actual variables with the mean
            
            end
        end
        
		% Dumps the output struct from the parallel simulation mode into the results object
        function set_TTI_results(obj,tmp_results)
            for ii = 1:size(obj.SNR_vector,2)
                for uu = 1:(obj.nUE * obj.nBS)

                    obj.UE_specific(uu).ACK(:,ii,:)                 = tmp_results(ii).UE_specific(uu).ACK(:,:);
                    obj.UE_specific(uu).rv_idx(:,ii,:)              = tmp_results(ii).UE_specific(uu).rv_idx(:,:);
                    obj.UE_specific(uu).RBs_assigned(:,ii)          = tmp_results(ii).UE_specific(uu).RBs_assigned(:);
                    obj.UE_specific(uu).biterrors_coded(:,ii,:)     = tmp_results(ii).UE_specific(uu).biterrors_coded(:,:);
                    obj.UE_specific(uu).biterrors_uncoded(:,ii,:)   = tmp_results(ii).UE_specific(uu).biterrors_uncoded(:,:);
                    obj.UE_specific(uu).blocksize_coded(:,ii,:)     = tmp_results(ii).UE_specific(uu).blocksize_coded(:,:);
                    obj.UE_specific(uu).blocksize_uncoded(:,ii,:)   = tmp_results(ii).UE_specific(uu).blocksize_uncoded(:,:);
                    obj.UE_specific(uu).throughput_coded(:,ii,:)    = tmp_results(ii).UE_specific(uu).throughput_coded(:,:);
                    obj.UE_specific(uu).throughput_uncoded(:,ii,:)  = tmp_results(ii).UE_specific(uu).throughput_uncoded(:,:);
                    obj.UE_specific(uu).FER_coded(:,ii,:)           = tmp_results(ii).UE_specific(uu).FER_coded(:,:);
                    obj.UE_specific(uu).FER_uncoded(:,ii,:)         = tmp_results(ii).UE_specific(uu).FER_uncoded(:,:);
                    obj.UE_specific(uu).used_codewords(:,ii,:)      = tmp_results(ii).UE_specific(uu).used_codewords(:,:);
                    obj.UE_specific(uu).channel_error(:,ii)         = tmp_results(ii).UE_specific(uu).channel_error(:);
                    obj.UE_specific(uu).channel_pred_error(:,ii)    = tmp_results(ii).UE_specific(uu).channel_pred_error(:);
                    obj.UE_specific(uu).papr(:,ii,:)                = tmp_results(ii).UE_specific(uu).papr(:,:);
                    obj.UE_specific(uu).used_cqi(:,ii,:)            = tmp_results(ii).UE_specific(uu).used_cqi(:,:);
                    obj.UE_specific(uu).used_RI(:,ii)               = tmp_results(ii).UE_specific(uu).used_RI(:);
                    obj.UE_specific(uu).type{ii}                    = tmp_results(ii).UE_specific(uu).type;
                    obj.UE_specific(uu).TTI_origin(:,ii)            = tmp_results(ii).UE_specific(uu).TTI_origin;
                    obj.UE_specific(uu).delay_TTI(:,ii)             = tmp_results(ii).UE_specific(uu).delay_TTI;
                    obj.UE_specific(uu).ID_count_current(ii)        = tmp_results(ii).UE_specific(uu).ID_count_current;
                    obj.UE_specific(uu).ID_count_next(ii)           = tmp_results(ii).UE_specific(uu).ID_count_next;
                end
                
                for bb = 1:obj.nBS
                    obj.cell_specific(bb).biterrors_coded(:,ii,:)    = tmp_results(ii).cell_specific(bb).biterrors_coded(:,:);
                    obj.cell_specific(bb).biterrors_uncoded(:,ii,:)  = tmp_results(ii).cell_specific(bb).biterrors_uncoded(:,:);
                    obj.cell_specific(bb).blocksize_coded(:,ii,:)    = tmp_results(ii).cell_specific(bb).blocksize_coded(:,:);
                    obj.cell_specific(bb).blocksize_uncoded(:,ii,:)  = tmp_results(ii).cell_specific(bb).blocksize_uncoded(:,:);
                    obj.cell_specific(bb).throughput_coded(:,ii,:)   = tmp_results(ii).cell_specific(bb).throughput_coded(:,:);
                    obj.cell_specific(bb).throughput_uncoded(:,ii,:) = tmp_results(ii).cell_specific(bb).throughput_uncoded(:,:);
                    obj.cell_specific(bb).FER_coded(:,ii,:)          = tmp_results(ii).cell_specific(bb).FER_coded(:,:);
                    obj.cell_specific(bb).FER_uncoded(:,ii,:)        = tmp_results(ii).cell_specific(bb).FER_uncoded(:,:);
                    obj.cell_specific(bb).used_codewords(:,ii,:)     = tmp_results(ii).cell_specific(bb).used_codewords(:,:);
                end
            end
        end
        
        function plot_BLER_throughput(obj,varargin)
            if isempty(varargin)
                first_figure = figure;
                second_figure = figure;
            else
                first_figure  = varargin{1}(1);
                second_figure = varargin{1}(2);
            end
            
            first_axes = axes('Parent',first_figure);
            second_axes = axes('Parent',second_figure);
            
            N_subframes = size(obj.cell_specific.used_codewords,1); 
            plot(first_axes,obj.SNR_vector,obj.cell_specific.BLER_overall);
            plot(second_axes,obj.SNR_vector,squeeze(sum(sum(obj.cell_specific.throughput_coded,1),3))/(N_subframes*1e-3)/1e6);
            
            set(first_axes,'YScale','log');
            ylabel(first_axes,'BLER');
            xlabel(first_axes,'SNR [dB]');
            ylim(first_axes,[1e-3 1]);
            grid(first_axes,'on');
            
            
            ylabel(second_axes,'throughput [Mbps]');
            xlabel(second_axes,'SNR [dB]');
            grid(second_axes,'on');
            
        end
    end
end