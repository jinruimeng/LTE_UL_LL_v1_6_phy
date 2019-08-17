function s = LTE_common_scrambling(b, NIDcell, NIDue, subframe, codeword, mode)
% LTE_common_scrambling - to scramble and descramble the coded bits in each of the codewords for PUSCH
% DEFINED IN STANDARD 3GPP TS 36.211 V11.4.0 (2013-09) Section 5.3.1
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input :   b             ... [1 x # coded bits]logical - coded bits from one codeword
%           NIDcell       ... [1 x 1]double cell identity
%           NIDue         ... [1 x 1]double UE identity
%           subframe      ... [1 x 1]double number of the subframe transmitted
%           codeword      ... [1 x 1]double number of the codeword transmitted
%           mode          ... string 'scramble' or 'descramble', in  LTE_TX for scrambling, 
%                             in LTE_RX for descrambling of bits in case of uncoded and LLRs in case of coded performance
% output:   s             ... [1 x # coded bits]logical/double - scrambled bits (in LTE_TX) or descrambled bits/LLRs (LTE_RX)  from one codeword


%% Generation of Pseudo-random sequence
% DEFINED IN STANDARD 3GPP TS 36.211 V11.4.0 (2013-09) Section 7.2

pn_seq = LTE_common_gen_gold_sequence(length(b),NIDcell,NIDue,subframe,codeword);

%% Bit scrambling / descrambling
% DEFINED IN STANDARD 3GPP TS 36.211 V11.4.0 (2013-09) Section 5.3.1
% placeholders for ARQ-ACK and RI bits not yet included
switch mode
    case 'scramble'
        s = logical(mod(double(b) + double(pn_seq), 2));
    case 'descramble'
        pn_seq = +pn_seq; % to convert from logical to double because here LLRs (1,-1) are processed
        pn_seq(pn_seq == 0) = -1;
        s = -(b.*pn_seq');
end