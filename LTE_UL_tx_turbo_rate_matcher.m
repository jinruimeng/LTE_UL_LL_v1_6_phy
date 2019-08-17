function ek = LTE_UL_tx_turbo_rate_matcher(LTE_params,dk,UE_signaling,UE,stream_index, uplink_user)
% LTE Turbo Code Rate Matcher, as of TS 36.212, Section 5.1.4.1.
% [ek UE_signaling] = LTE_tx_turbo_rate_matcher(dk,UE_signaling,UE)
% Author: Josep Colom Ikuno, josep.colom@nt.tuwien.ac.at
% added uplink modification Prokopec
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input :   input_bits        ... Coded bits dk, as outputted by the turbo 
%                                 encoder This is a matrix of 3xN bits.
%           UE_signaling      ... BS signaling
% output:   output_bits       ... rate-matched bits
%           UE_signaling      ... UE_signaling with extra information, such as
%                                 number of padding bits used for
%                                 the interleaving (foe each sub-block).
%                                 Needed for the rate matching at the
%                                 receiver.

% Nomenclature scheme
%        _______________________        ________________        ___________________________
%  dk-> | sub-block interleaver | vk-> | bit collection | wk-> | bit selection and pruning | ek->
%       |_______________________|      |________________       |___________________________|
%
% date of creation: 2008/08/11
% last changes:
%   2009/04/22  Jcolom      Changed how N_IR is set according to the new version of the standard (8.6.0)

if UE.mode==3 || UE.mode==4
    K_MIMO = 2;
else
    K_MIMO = 1;
end
M_limit = 8;

N_IR = floor(UE.N_soft / (K_MIMO*min(LTE_params.HARQ_processes,M_limit)));

UE_signaling.turbo_rate_matcher(stream_index).total_bits_before_punturing = 0;
length_dk = length(dk);
vk = cell(1,length_dk);
wk = cell(1,length_dk);
ek = cell(1,length_dk);
for i=1:length(dk)
    % Sub-block interleaver
    [   vk{i},...
        UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).Nd,...
        UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).R_tc,...
        UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).K_pi ] = LTE_common_subblock_interleaver(LTE_params,dk{i},1);
    UE_signaling.turbo_rate_matcher(stream_index).total_bits_before_punturing = UE_signaling.turbo_rate_matcher(stream_index).total_bits_before_punturing + 3*UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).K_pi;
    % Circular buffer
    wk{i} = LTE_common_turbo_rate_matcher_circular_buffer(vk{i},1);
      
    % Bit selection and pruning (first get the mapping)
    % for uplink, multiple UE's need this change of indexing

    [bit_selection_and_pruning_mapping k_0 null_positions real_number_used] = LTE_common_turbo_rate_matcher_bit_selection_and_pruning_mapping(...
        UE_signaling.turbo_rate_matcher(stream_index).G,...
        N_IR,...
        UE_signaling.MCS_and_scheduling(stream_index).CQI_params.modulation_order,...
        UE_signaling.turbo_rate_matcher(stream_index).N_l,...
        UE_signaling.turbo_rate_matcher(stream_index).rv_idx,...
        UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).K_pi,...
        UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).R_tc,...
        UE_signaling.TB_segmentation(stream_index).C,...
        UE_signaling.TB_segmentation(stream_index).F,...
        i-1,...
        UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).Nd,...
        1,...
        wk{i});


    null_positions = null_positions(1:real_number_used); % This is needed as not all filler bits may be in the output
    
    UE_signaling.turbo_rate_matcher(stream_index).null_positions{i} = null_positions; % Please take note that this is 0_indexed!
    UE_signaling.turbo_rate_matcher(stream_index).ek_sizes(i) = length(bit_selection_and_pruning_mapping);
    % The actual bit selection and pruning (1=interleave, as interleaving/punturing and its inverse are
    % implemented in the same function)
    ek{i} = LTE_common_turbo_rate_matching_bit_selection_and_pruning(wk{i},bit_selection_and_pruning_mapping,1);
end

