function [BS_output,LTE_params] = LTE_tx_PCFICH(LTE_params,BS,b_,subframe_corr,UE,BS_output) 
%% PCFICH processing 
% TS 36.211 V8.9.0 Section 6.7
% TS 36.212 V8.8.0 Section 5.3.4 
% Author: Petr Kej�k, xkejik00@stud.feec.vutbr.cz
% 
% input : Nrb           ... [1 x 1] double number of resource blocks
%         b_            ... [1 x 1] double number - b-th NodeB
%         subframe_corr ... [1 x 1] double number - actual subframe number
%         Nsc           ... [1 x 1] double number - resource block size
%         Nsub
%         Ns            ... [1 x 1] double number - number of slots
% output: CFIcw         ... [1 x 32]double generated CFI codeword

 CFI = 2;

%% Channel coding
% TS 36.212 V8.8.0 Section 5.3.4
if LTE_params.Nrb <= 10
    CFI_i = CFI + 1;
else
    CFI_i = CFI;
end

switch CFI_i
    case 1
        CFIcw = [0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1];
    case 2
        CFIcw = [1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0];
    case 3
        CFIcw = [1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1];
    otherwise
        CFIcw = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
        error('This combination is reseved (LTE-tx-PCFICH)')
end
LTE_params.PCFICH(b_,subframe_corr).CFIcodeword = CFIcw;
LTE_params.PCFICH(b_,subframe_corr).txCFI = CFI_i;

%% Scrambling 
% TS 36.211 V8.9.0 Section 6.7.1
c_init = (floor((2*subframe_corr-2)/2)+1)*(2*BS.NIDcell+1)*(2^9)+BS.NIDcell; % initialization value
% (2*subframe_corr-2) = slot number - only even slots
pn_seq = LTE_common_gen_gold_sequence_CCH(length(CFIcw),c_init); % pn generation
CFIcwSC = logical(mod(CFIcw + pn_seq, 2)); % scrambling
LTE_params.PCFICH(b_,subframe_corr).CFIcwSC = CFIcwSC;

BS_output(b_).PCFICH_tx = CFI_i;
BS_output(b_).PCFICH_CFI = CFIcwSC;

%% Symbol mapping - QPSK modulation 
% TS 36.211 V8.9.0 Section 6.7.2
nibble = 2*LTE_params.PCFICH(b_,subframe_corr).CFIcwSC(1:2:end) + LTE_params.PCFICH(b_,subframe_corr).CFIcwSC(2:2:end);
tx_PCFICH_symbols = LTE_params.SymbolAlphabet{LTE_params.PCFICH_common_param.modulation_order}(nibble+1).';
LTE_params.PCFICH(b_,subframe_corr).CFIcwSym = tx_PCFICH_symbols;

%% Layer mapping 
% TS 36.211 V8.9.0 Section 6.7.3
nLayers = LTE_params.BS_config.nAtPort; % Section 6.3.3.3
switch nLayers % transmission mode 
    case 1  % single antenna transmission, TS 36.211, Section 6.3.3.1
        layer_xCFI = tx_PCFICH_symbols;
    case 2  % transmit diversity, TS 36.211, Section 6.3.3.3
        layer_xCFI(1,:) = tx_PCFICH_symbols(1:2:end);
        layer_xCFI(2,:) = tx_PCFICH_symbols(2:2:end);
    case 4
        if (mod(length(tx_PCFICH_symbols),4)~=0)
           tx_PCFICH_symbols = [tx_PCFICH_symbols,0,0];
        end  
        layer_xCFI(1,:) = tx_PCFICH_symbols(1:4:end);
        layer_xCFI(2,:) = tx_PCFICH_symbols(2:4:end);
        layer_xCFI(3,:) = tx_PCFICH_symbols(3:4:end);
        layer_xCFI(4,:) = tx_PCFICH_symbols(4:4:end);
    otherwise
        error('Number of layers not supported (LTE-tx-PCFICH)');
end
LTE_params.PCFICH(b_,subframe_corr).CFIcwSymLM = layer_xCFI;

%% Precoding 
% TS 36.211 V8.9.0 Section 6.3.4.1 and 6.3.4.3
AtPorts = LTE_params.BS_config.nAtPort;
layer_x = layer_xCFI;
codebook_index = 1; % no meaning for PCFICH - spatial multiplexing is not used for PCFICH
RI = 1;             % no meaning for PCFICH
CDD = 1;            % no meaning for PCFICH 
rb_numbers = 1;     % no meaning for PCFICH
if UE.mode > 2; % spatial multiplexing is not used for PCFICH
    error('This configuration is not implemented yet - spatial multiplexing is not used for PCFICH (LTE-tx-PCFICH)')
end

[precode_y, UE, D, W, U] = LTE_tx_precoding(LTE_params, layer_x, UE, AtPorts, codebook_index, RI, CDD,rb_numbers);

%% Preliminary symbol assembling - see LTE_common_gen_PCFICH.m
% NIDcell shifts PCFICH mapping - 1st quadruplet is not always mapped to
% the first resource elements group - I am not sure about that
k_sort = LTE_params.PCFICH(b_,subframe_corr).k_sort;
precode_y_temp = zeros(LTE_params.BS_config.nAtPort,16);
for i6 = 1:4 % four REs groups, four REs in each group
     i7 = k_sort(i6,1);
     precode_y_temp(:,(i6-1)*4+1:(i6-1)*4+4) = precode_y(:,(i7-1)*4+1:(i7-1)*4+4);
 end

LTE_params.PCFICH(b_,subframe_corr).CFIcwSymPrecode = precode_y_temp;

