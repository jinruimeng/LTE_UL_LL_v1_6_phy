function [rx_scrambled_bits_F,CFI_i,LTE_params,UE_signaling]=LTE_rx_PCFICH(UE_signaling,LTE_params,H_est,subframe_corr,uu,UE,BS,y_rx_assembled,ChanMod,sigma_n2)
%% PCFICH processing 
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

%% Warnings, comments, possible changes
% disassembe symbols part is not clear
  b_ = 1; % Node Bs number

%% Modification of some parameters I
% UEs use scheduled modulation schemes but PCFICH uses QPSK
modul_order_temp = UE_signaling.MCS_and_scheduling.CQI_params.modulation_order;
UE_signaling.MCS_and_scheduling.CQI_params.modulation_order = LTE_params.PCFICH_common_param.modulation_order;

%% Disassemble symbols
% resource elements used for PCFICH
PCFICH_ResEl = LTE_params.PCFICH(b_,subframe_corr).Mapping(:,1); 

for nn = 1:UE(uu).nRX
    for tt = 1:BS.nTX 
        H_temp = (H_est(:,1,nn,tt));
        H_est_user_a = H_temp(logical(PCFICH_ResEl),:,:);
        H_est_user_F(:,nn,tt) = H_est_user_a;
    end
    y_rx_assembled_temp = y_rx_assembled(:,1,nn);
    rx_user_symbols_a = y_rx_assembled_temp(logical(PCFICH_ResEl));
    rx_user_symbols_F(:,nn) = rx_user_symbols_a;
end

%% Disassemble symbols (preliminary) - see LTE_tx_PCFICH.m
% k_sort=[2; 3; 4; 1];
k_sort = LTE_params.PCFICH(b_,subframe_corr).k_sort;
rx_user_symbols_F_temp = zeros(16,UE(uu).nRX);
H_est_user_F_temp = zeros(16,UE(uu).nRX,BS.nTX);
for i6 = 1:4 % four REs groups
     i7 = find(k_sort(:,1)==i6);
     rx_user_symbols_F_temp((i6-1)*4+1:(i6-1)*4+4,:) = rx_user_symbols_F((i7-1)*4+1:(i7-1)*4+4,:);
     H_est_user_F_temp((i6-1)*4+1:(i6-1)*4+4,:,:) = H_est_user_F((i7-1)*4+1:(i7-1)*4+4,:,:);
end
 rx_user_symbols_F = rx_user_symbols_F_temp;
 H_est_user_F = H_est_user_F_temp;

%% PCFICH detecting
LTE_params.PCFICH_common_param.help_index = 1; % indicates PCFICH mapping to detectors - see LTE_detect_SIxO.m 
[LLR_SD_F,M_F] = LTE_detecting(UE_signaling.MCS_and_scheduling,BS.nAtPort,rx_user_symbols_F,UE(uu).nRX,H_est_user_F,ChanMod.filtering,H_est,LTE_params,UE(uu).receiver,sigma_n2);
LTE_params.PCFICH_common_param.help_index = 0;

%% Undo layer mapping
switch LTE_params.BS_config.nAtPort
    case 1  % single antenna transmission
        LLR_SS_F{1} = reshape(LLR_SD_F(M_F(1):-1:1,:),1,[]).';
    case 2
        LLR_SS_F{1} = reshape([LLR_SD_F(M_F(1):-1:1,:);LLR_SD_F(M_F(1)+M_F(2):-1:M_F(1)+1,:)],1,[]).';
    case 4
        LLR_SS_F{1} = reshape([LLR_SD_F(M_F(1):-1:1,:);LLR_SD_F(M_F(1)+M_F(2):-1:M_F(1)+1,:);LLR_SD_F(M_F(1)+M_F(2)+M_F(3):-1:M_F(1)+M_F(2)+1,:);LLR_SD_F(M_F(1)+M_F(2)+M_F(3)+M_F(4):-1:M_F(1)+M_F(2)+M_F(3)+1,:)],1,[]).';
    otherwise
        error('This combination is not valid for PCFICH (LTE-rx-PCFICH)')
end
rx_scrambled_bits_F = (1+sign(LLR_SS_F{1}.'))/2;

%% Descrambling 
% TS 36.211, Section 6.7.1
% initialization value
c_init = (floor((2*subframe_corr-2)/2)+1)*(2*BS.NIDcell+1)*(2^9)+BS.NIDcell;
pn_seq = LTE_common_gen_gold_sequence_CCH(length(rx_scrambled_bits_F),c_init); % pn generation
CFIcwSC_RX = logical(mod(rx_scrambled_bits_F + pn_seq, 2)); % descrambling

%% CFI detection
CFIcw = zeros(4,32);
CFI_pos = zeros(4,1);
CFIcw(1,:) = [0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1];
CFIcw(2,:) = [1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0];
CFIcw(3,:) = [1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1];
CFIcw(4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

for i1 = 1:4 % "matched filter approach"
    CFI_pos(i1,1) = sum(CFIcwSC_RX == CFIcw(i1,:)); 
end
% final decision about CFI
CFI_i = find(CFI_pos(:,1) == max(CFI_pos));
LTE_params.PCFICH(b_,subframe_corr).rxCFI = CFI_i;

%% Modification of some parameters II
% modulation order has to be restored
UE_signaling.MCS_and_scheduling.CQI_params.modulation_order = modul_order_temp;


    
    
    
    
    