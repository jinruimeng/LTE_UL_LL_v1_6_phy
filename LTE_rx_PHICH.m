function [LTE_params,UE_output] = LTE_rx_PHICH(UE_signaling,LTE_params,H_est,subframe_corr,uu,UE,BS,y_rx_assembled,...
         ChanMod,sigma_n2,UE_output,BS_output)
%% PHICH processing 
% (TS 36.211 V8.9.0 Section 6.7) 
% (TS 36.212 V8.8.0 Section 5.3.4) 
% Author: Petr Kej�k, xkejik00@stud.feec.vutbr.cz
% 
% input : Nrb           ... [1 x 1] double number of resource blocks
%         b_            ... [1 x 1] double number - bth NodeB
%         subframe_corr ... [1 x 1] double number - actual subframe number
%         Nsc           ... [1 x 1] double number - resource block size
%         Nsub
%         Ns            ... [1 x 1] double number - number of slots
% output: CFIcw         ... [1 x 32]double generated CFI codeword
 
% Vys�laj� se i nuly, QPSK je detekuje :-(
%% Warnings, comments, possible changes
% disassembe symbols part is not clear
b_ = 1; % Node Bs number

%% Modification of some parameters I
% UEs use scheduled modulation schemes but PCFICH uses BPSK
modul_order_temp = UE_signaling.MCS_and_scheduling.CQI_params.modulation_order;
UE_signaling.MCS_and_scheduling.CQI_params.modulation_order = 2; % QPSK !!! constellation 
% diagram changes due to spreading with complex sequences (+j and -j) 
 
%% Disassemble symbols
% resource elements used for PHICH
PHICH_ResEl = LTE_params.PHICH(b_,subframe_corr).Mapping(:,1); 
 
for nn = 1:UE(uu).nRX
    for tt = 1:BS.nTX 
        H_temp = (H_est(:,1,nn,tt));
        H_est_user_a = H_temp(logical(PHICH_ResEl),:,:);
        H_est_user_HI(:,nn,tt) = H_est_user_a;
    end
    y_rx_assembled_temp = y_rx_assembled(:,1,nn);
    rx_user_symbols_a = y_rx_assembled_temp(logical(PHICH_ResEl));
    rx_user_symbols_HI(:,nn) = rx_user_symbols_a;
end

%% PHICH detecting
LTE_params.PHICH_common_param.help_index = 1; % indicates PHICH mapping to detectors - see LTE_detect_SIxO.m 
[LLR_SD_HI,M_HI] = LTE_detecting(UE_signaling.MCS_and_scheduling,BS.nAtPort,rx_user_symbols_HI,UE(uu).nRX,H_est_user_HI,ChanMod.filtering,H_est,LTE_params,...
    UE(uu).receiver,sigma_n2);
LTE_params.PHICH_common_param.help_index = 0;
 
%% Undo layer mapping
switch LTE_params.BS_config.nAtPort
    case 1  % single antenna transmission
        LLR_SS_HI{1} = reshape(LLR_SD_HI(M_HI(1):-1:1,:),1,[]).';
    case 2
        LLR_SS_HI{1} = reshape([LLR_SD_HI(M_HI(1):-1:1,:);LLR_SD_HI(M_HI(1)+M_HI(2):-1:M_HI(1)+1,:)],1,[]).';
    case 4
        LLR_SS_HI{1} = reshape([LLR_SD_HI(M_HI(1):-1:1,:);LLR_SD_HI(M_HI(1)+M_HI(2):-1:M_HI(1)+1,:);LLR_SD_HI(M_HI(1)+M_HI(2)+M_HI(3):-1:M_HI(1)+M_HI(2)+1,:);LLR_SD_HI(M_HI(1)+M_HI(2)+M_HI(3)+M_HI(4):-1:M_HI(1)+M_HI(2)+M_HI(3)+1,:)],1,[]).';
    otherwise
        error('This combination is not valid for PCFICH (LTE-rx-PCFICH)')
end
rx_bits_HI = (1+sign(LLR_SS_HI{1}.'))/2;

%% Symbol mapping - for despreading and decsreambling puprose
nibble = 2 * rx_bits_HI(1:2:end) + rx_bits_HI(2:2:end);
HI_rx_symbols = LTE_params.SymbolAlphabet{2}(nibble+1).';

%%
N_g=1;
switch LTE_params.CyclicPrefix         
    case 'normal'
        N_g_phich = ceil(N_g*(LTE_params.Nrb/8));
        N_SF = 4;
        w_phich = [ +1,  +1,  +1,  +1;
                    +1,  -1,  +1,  -1;
                    +1,  +1,  -1,  -1;
                    +1,  -1,  -1,  +1;
                   +1j, +1j, +1j, +1j;
                   +1j, -1j, +1j, -1j;
                   +1j, +1j, -1j, -1j;
                   +1j, -1j, -1j, +1j];
    case 'extended'
        N_g_phich = 2*ceil(N_g*(LTE_params.Nrb/8));
        N_SF = 2;
        w_phich = [ +1,  +1;
                    +1,  -1;
                   +1j, +1j;
                   +1j, -1j];
    otherwise
        error('Wrong cyclic prefix (LTE-tx-PHICH)')
end

%% Despreading and descrambling
% (TS 36.211 V8.9.0 Section 6.9.1)
c_init = (floor((2*subframe_corr-2)/2)+1)*(2*BS.NIDcell+1)*(2^9)+BS.NIDcell; % initialization value
pn_seq = LTE_common_gen_gold_sequence_CCH(N_SF*3,c_init); % pn generation

HI_rx_despr = HI_rx_symbols.*(1-2*pn_seq);

for ii=2:2
    HI_rx_descr = HI_rx_despr.* [conj(w_phich(ii,:)) conj(w_phich(ii,:)) conj(w_phich(ii,:))];
end

%% HI detection
HI_rx_dec = zeros(1,12);
for iii=1:12
    if real(HI_rx_descr(1,iii)) < 0 && imag(HI_rx_descr(1,iii)) < 0 
        HI_rx_dec(1,iii) = 1;
    elseif real(HI_rx_descr(1,iii)) > 0 && imag(HI_rx_descr(1,iii)) > 0 
        HI_rx_dec(1,iii) = -1;
    else
        HI_rx_dec(1,iii) = 0;
    end
end
HI_rx_sum = sum(HI_rx_dec)/12;

if HI_rx_sum < 0; 
    HI_rx = 0;
else
    HI_rx = 1;
end
UE_output(uu).PHICH_HI_rx = HI_rx;

 
%% Modification of some parameters II
% modulation order has to be restored
UE_signaling.MCS_and_scheduling.CQI_params.modulation_order = modul_order_temp;

   
     