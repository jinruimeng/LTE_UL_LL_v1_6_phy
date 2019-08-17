function [  UE_ACK,...
            UE_rv_idx,...
            UE_RBs_assigned,...
            UE_biterrors_coded,...
            UE_biterrors_uncoded,...
            UE_blocksize_coded,...
            UE_blocksize_uncoded,...
            UE_throughput_coded,...
            UE_FER_coded,...
            UE_throughput_uncoded,...
            UE_FER_uncoded,...
            UE_used_codewords,...
            UE_channel_error,...
            UE_channel_pred_error,...
            UE_papr,...
            cell_biterrors_coded,...
            cell_biterrors_uncoded,...
            cell_blocksize_coded,...
            cell_blocksize_uncoded,...
            cell_throughput_coded,...
            cell_throughput_uncoded,...
            cell_FER_coded,...
            cell_FER_uncoded,...
            cell_used_codewords, ...
            UE_used_cqi, ...
            UE_used_RI] = LTE_UL_processTTI_results(   BS_output,...
                                                                UE_output,...
                                                                subframe_i,...
                                                                SNR_i,...
                                                                nUE,...
                                                                nBS, ...
                                                                cell_biterrors_coded,...
                                                                cell_biterrors_uncoded,...
                                                                cell_blocksize_coded,...
                                                                cell_blocksize_uncoded,...
                                                                cell_throughput_coded,...
                                                                cell_throughput_uncoded,...
                                                                cell_FER_coded,...
                                                                cell_FER_uncoded,...
                                                                cell_used_codewords,...
                                                                maxstreams,...
                                                                LTE_params_tmp,...
                                                                connection_table)
% Process temporary results for parallel simulations.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

% Initialization
UE_ACK                  = false(maxstreams,nUE*nBS);
UE_rv_idx               = zeros(maxstreams,nUE*nBS,'uint8');
UE_RBs_assigned         = zeros(nUE*nBS,1,'uint8');
UE_biterrors_coded      = zeros(maxstreams,nUE*nBS,'uint32');
UE_biterrors_uncoded    = zeros(maxstreams,nUE*nBS,'uint32');
UE_blocksize_coded      = zeros(maxstreams,nUE*nBS,'uint32');
UE_blocksize_uncoded    = zeros(maxstreams,nUE*nBS,'uint16');
UE_throughput_coded     = zeros(maxstreams,nUE*nBS,'uint32');
UE_throughput_uncoded   = zeros(maxstreams,nUE*nBS,'uint32');
UE_FER_coded            = zeros(maxstreams,nUE*nBS,'uint16');
UE_FER_uncoded          = zeros(maxstreams,nUE*nBS,'uint16');
UE_used_codewords       = zeros(maxstreams,nUE*nBS,'uint16');
UE_used_cqi             = zeros(maxstreams,nUE*nBS,'uint32');
UE_used_RI            = zeros(nUE*nBS, 1,'uint32');
UE_channel_error        = zeros(nUE*nBS,1,'double');
UE_channel_pred_error   = zeros(nUE*nBS,1,'double');
UE_papr                 = zeros(LTE_params_tmp.Nsub,nUE*nBS,'double');

% Loop over all UEs and streams
for uu = 1:nUE*nBS
    bb = utils.findBS(connection_table, uu);
    local_uu = utils.globalToLocalUser(connection_table, bb, uu);
    
    for stream_i = 1:BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.nCodewords
        
        % Update ACK (UE traces) and biterrors (cell traces)
        UE_ACK(stream_i,uu)       = UE_output(uu).ACK(stream_i);
        UE_rv_idx(stream_i,uu)    = UE_output(uu).rv_idx(stream_i);
        UE_RBs_assigned(uu)       = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.assigned_RBs;
        UE_used_cqi(stream_i, uu) = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.cqi(stream_i);
        UE_used_RI(uu)            = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.nLayers;
        
       
        if UE_output(uu).UE_scheduled
           
            UE_biterrors_coded(stream_i,uu)   = uint32(sum(abs(int32(UE_output(uu).rx_data_bits{stream_i})  - int32(UE_output(uu).UE_genie.data_bits{stream_i}))));
            UE_biterrors_uncoded(stream_i,uu) = uint32(sum(abs(int32(UE_output(uu).rx_coded_bits{stream_i}) - int32(UE_output(uu).UE_genie.sent_bits{stream_i}))));

            % Update blocksize (UE traces)
            UE_blocksize_coded(stream_i,uu)   = length(UE_output(uu).rx_data_bits{stream_i});
            UE_blocksize_uncoded(stream_i,uu) = length(UE_output(uu).rx_coded_bits{stream_i});
            
            % Update FER and throughput (UE traces)
            
            % Coded
            if UE_output(uu).ACK(stream_i)
                UE_throughput_coded(stream_i,uu) = length(UE_output(uu).rx_data_bits{stream_i});
                UE_FER_coded(stream_i,uu) = 0;
            else
                UE_FER_coded(stream_i,uu) = 1;
                UE_throughput_coded(stream_i,uu) = 0;
            end
            
            % Uncoded
            if (UE_biterrors_uncoded(stream_i,uu)==0)%!
                UE_throughput_uncoded(stream_i,uu) = length(UE_output(uu).rx_coded_bits{stream_i});
                UE_FER_uncoded(stream_i,uu) = 0;
            else
                UE_FER_uncoded(stream_i,uu) = 1;
                UE_throughput_uncoded(stream_i,uu) = 0;
            end
            
            % update papr
            UE_papr(1:size(UE_output(uu).papr,1),uu) = UE_output(uu).papr;

        end
        
        % Update what codewords were used (valid positions in the traces, 0 if no RBs were allocated)
        if(BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.assigned_RBs)
            UE_used_codewords(stream_i,uu) = 1;
        else
            UE_used_codewords(stream_i,uu) = 0;
        end
        
    end
    
    
    if UE_output(uu).UE_scheduled
        % update channel estimation error
       UE_channel_error(uu) = UE_channel_error(uu) + UE_output(uu).channel_estimation_error;
        % update channel prediction error
       UE_channel_pred_error(uu) = UE_channel_pred_error(uu) + UE_output(uu).channel_prediction_error;
    end
    
    if utils.isUEattached( connection_table, bb, uu )
        % Update cell coded and uncoded bit errors
        cell_biterrors_coded(:, :, bb)     = cell_biterrors_coded(:, :, bb) + UE_biterrors_coded(:,uu).';
        cell_biterrors_uncoded(:, :, bb)   = cell_biterrors_uncoded(:, :, bb) + UE_biterrors_uncoded(:,uu).';

        % Update blocksize (cell traces)
        cell_blocksize_coded(:,:, bb)     = cell_blocksize_coded(:, :,  bb)  + UE_blocksize_coded(:,uu).';
        cell_blocksize_uncoded(:,:, bb)   = uint32(cell_blocksize_uncoded(:,:, bb)) + uint32(UE_blocksize_uncoded(:,uu)).';

        % Update FER and throughput (cell traces)
        cell_throughput_coded(:,:, bb)    = cell_throughput_coded(:,:, bb)      + UE_throughput_coded(:,uu).';
        cell_throughput_uncoded(:,:, bb)  = cell_throughput_uncoded(:,:, bb)    + UE_throughput_uncoded(:,uu).';
        cell_FER_coded(:,:, bb)           = uint16(cell_FER_coded(:,:, bb))     + uint16(UE_FER_coded(:,uu)).';
        cell_FER_uncoded(:,:, bb)         = uint16(cell_FER_uncoded(:,:, bb))   + uint16(UE_FER_uncoded(:,uu)).';

        % The number of codewords is the maximum (bitwise OR) of the used codewords matrix
        cell_used_codewords(:,:, bb) = uint16(cell_used_codewords(:,:, bb)) + uint16(UE_used_codewords(:,uu)).';
    end
end
        