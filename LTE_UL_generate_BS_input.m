function [BS_input,n_rstream, signal_size] = LTE_UL_generate_BS_input(LTE_params, ChanMod_output, SNR, UE_MCS_and_scheduling_info)
% LTE_UL_generate_BS_input - combine ChanMod_output to BS_input
% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input  :   LTE_params       ... generic parameter struct
%            ChanMod_output   ... output of each UE to each BS 
%            SNR              ... current SNR (dB)
% 
% output :   BS_input         ... input signal struct for each BS
%       
nBS = LTE_params.nBS;
nUE = LTE_params.nUE;
n_rstream = LTE_params.noise_RandStream;
pathloss_matrix = LTE_params.pathloss_matrix;

maxStreams = min(LTE_params.UE_config.nTX,LTE_params.BS_config.nRX);

% obtain signal length - implicit assumtion here is that all signals have
% the same length
for bb=1:nBS
    for uu=1:(nBS*nUE)
        if utils.isScheduled(UE_MCS_and_scheduling_info, bb, uu, LTE_params.connection_table)
            signal_size = size(ChanMod_output{bb,uu}.y_rx);
            break
        end
    end
end




% the main struct: BS_input
BS_input = cell(nBS, maxStreams);



for bb = 1:nBS
%%
    % first add the noise
   %BS_input{bb}.input = ( randn(n_rstream, signal_size) + 1i*randn(n_rstream, signal_size) ) * 10^(-SNR/20)/sqrt(2);
    
%%
    % weighted addition of received signals from all UEs
    for uu = 1:(nBS * nUE)
        if isinf(pathloss_matrix(bb,uu)) == 0 % no infinite pathloss 
            if utils.isScheduled(UE_MCS_and_scheduling_info, 1:LTE_params.nBS, uu, LTE_params.connection_table)
            % BS_input{bb}.input = BS_input{bb}.input + ChanMod_output{bb,uu}.y_rx * 10^(-pathloss_matrix(bb,uu)/20);
            BS_input{bb}.input =  ChanMod_output{bb,uu}.y_rx * 10^(-pathloss_matrix(bb,uu)/20);
            end 
        end 
    end
end 

end % LTE_UL_generate_BS_input






