function dk = LTE_UL_rx_turbo_rate_matcher(LTE_params,ek,BS_output,UE,stream_index,uplink_user)
% LTE Turbo Code Rate Matcher, as of TS 36.212, Section 5.1.4.1.
% [dk] = LTE_rx_turbo_rate_matcher(ek,UE_signaling,UE)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input :   input_bits        ... ek (soft bits)
%           UE_signaling      ... BS signaling
% output:   output_bits       ... rate de-matched bits

% Nomenclature scheme
%        _______________________        ________________        ___________________________
%  dk-> | sub-block interleaver | vk-> | bit collection | wk-> | bit selection and pruning | ek->
%       |_______________________|      |________________       |___________________________|
%
% date of creation: 2008/08/11
% last changes: 
%   2008/09/15  Bosanska    added input [1 x 1]struct UE
%   2009/04/22  Jcolom      Changed how N_IR is set according to the new version of the standard (8.6.0)
%               Blumenstein Upravil jsem LTE_rx_turbo_rate_matcher na LTE_rx_ul_turbo_rate_matcher.

if UE(1,uplink_user).mode==3 || UE(1,uplink_user).mode==4
    K_MIMO = 2;
else
    K_MIMO = 1;
end
M_limit = 8;

N_IR = floor(UE(1,uplink_user).N_soft / (K_MIMO*min(LTE_params.HARQ_processes,M_limit)));


uplink_user = 1;   

    [bit_selection_and_pruning_mapping k_0 null_positions] = LTE_common_turbo_rate_matcher_bit_selection_and_pruning_mapping(...
        BS_output.UE_signaling(1,uplink_user).turbo_rate_matcher(stream_index).G,...
        N_IR,...
        BS_output.UE_signaling(1,uplink_user).MCS_and_scheduling(1,uplink_user).CQI_params(stream_index).modulation_order,...
        BS_output.UE_signaling(1,uplink_user).turbo_rate_matcher(stream_index).N_l,...
        BS_output.UE_signaling(1,uplink_user).turbo_rate_matcher(stream_index).rv_idx,...
        BS_output.UE_signaling(1,uplink_user).turbo_rate_matcher(stream_index).subblock_interleaver.K_pi,...
        BS_output.UE_signaling(1,uplink_user).turbo_rate_matcher(stream_index).subblock_interleaver.R_tc,...
        BS_output.UE_signaling(1,uplink_user).TB_segmentation(stream_index).C,...
        BS_output.UE_signaling(1,uplink_user).TB_segmentation(stream_index).F,...
        1-1,...
        BS_output.UE_signaling(1,uplink_user).turbo_rate_matcher(stream_index).subblock_interleaver.Nd,...
        2,...
        BS_output.UE_signaling(1,uplink_user).turbo_rate_matcher(stream_index).null_positions);
    
    % Actual bit selection and pruning
    wk = LTE_common_turbo_rate_matching_bit_selection_and_pruning(double(ek),...
        bit_selection_and_pruning_mapping,...
        2,...
        BS_output.UE_signaling(1,uplink_user).turbo_rate_matcher(stream_index).subblock_interleaver.K_pi*3);
    
    % Circular buffer
    vk = LTE_common_turbo_rate_matcher_circular_buffer(wk,2);
    % vk{i} = LTE_common_turbo_rate_matcher_circular_buffer(ek{i},2);
    
    % Sub-block interleaver
    dk = LTE_common_subblock_interleaver(LTE_params,vk,2,BS_output.UE_signaling(1,uplink_user).turbo_rate_matcher(stream_index).subblock_interleaver.Nd);
    
BS_output.UE_signaling(1,uplink_user).turbo_rate_matcher(stream_index).ek_sizes(1) = length(bit_selection_and_pruning_mapping);
    % The actual bit selection and pruning (1=interleave, as interleaving/punturing and its inverse are
    % implemented in the same function)
    
% end
