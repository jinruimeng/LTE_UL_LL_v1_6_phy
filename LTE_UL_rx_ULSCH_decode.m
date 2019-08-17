function [a BER ACK] = LTE_UL_rx_ULSCH_decode(LTE_params,f,BS_signaling,BS,UE,BS_ue_specific,genie_data,stream_index)
% Performs the decoding of a TB for the LTE Uplink LL simulator. Basically an implementation
% of TS 36.212. Naming convention also according to TS 36.212.
% [a BER ACK ] = LTE_UL_rx_ULSCH_decode(f,BS_signaling,UE), 
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input:    f                   ... received LLRs (soft bits)
%           BS_signaling        ... struct containing parameters needed
%                                   for decoding (number of padding bits, 
%                                   TB size...)
%           UE                  ... struct that contains the number of 
%                                   turbo iterations that the user will
%                                   perform and some rate matching
%                                   params, such as the soft buffer size.
%
% output:   a                   ... decoded data bits
%           BER                 ... BER, as outputted by the turbo
%                                   decoder (using genie info). Actually
%                                   the number of erroneous bits per
%                                   turbo iteration
%           ACK                 ... ACK, calculated from the TB CRCs. A
%                                   true of false value stating if the TB
%                                   CRC check was correct (true) or
%                                   failed (false).
%           UE                  ... updated UE struct.
%
% Note: a BLER_b value containing the CRC checks of the individual code
% blocks is also computed, but neither used nor passed outside of the
% function.
%
% date of creation: 2008/08/11
% last changes: 2008/09/15 Bosanska     added input BS_ue_specific - [1 x 1]struct with tx user specific HARQ parameters

% Code block de-concatenation
e = LTE_UL_rx_code_block_concatenation(f,BS_signaling,stream_index);
% Rate de-matching
d = LTE_UL_rx_turbo_rate_matcher(LTE_params,e,BS_signaling,UE,stream_index);

% At the encoder the <NULL> bits are treated as zeros, thus we will now replace those values by negative LLRs (d^(0) and d^(1))
d{1}([1 2],1:BS_signaling.TB_segmentation_UL(stream_index).F) = -10;

% HARQ combining. Could have been embedded into the rate matcher, but it's less messy to put it separate
[ d_HARQ BS ] = LTE_rx_HARQ_combine(d,BS,BS_ue_specific,stream_index);

% Turbo decoding. Genie info is used to get BER values for each turbo
% iteration
[ c BER ] = LTE_rx_turbo_decode(LTE_params,d_HARQ,...
    BS.turbo_iterations,...
    BS_signaling.TB_segmentation_UL(stream_index).C,...
    genie_data.bits_to_turboencode{stream_index});
% end

% Code block desegmentation
[b ACK_b] = LTE_rx_code_block_desegmentation(c,BS_signaling.TB_segmentation_UL(stream_index).F);

% Transport block CRC check
ACK = LTE_rx_check_crc(b,'24a');

% Get the Transport Block data
a = b(1:(end-24));