function [LLR_SD,M,H_back] = LTE_detecting(MCS_and_scheduling,BS_nAtPort,rx_user_symbols,UE_nRX,H_est,filtering,LTE_params,receiver,sigma_n2,UE_output,uu, bb)

% Detection.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

LLR_SD = 0;
M = 0;
H_back = 0;



switch LTE_params.UE_config.mode 
    case 1      % SISO 
        if (~isempty(MCS_and_scheduling.freq_indices))
        [LLR_SD,M] = LTE_UL_detect_SIxO(MCS_and_scheduling,filtering,rx_user_symbols,H_est,LTE_params,receiver,sigma_n2,UE_output,uu, bb);
        end
    case 4     % Closed Loop Spatial Multiplexing 
        if (~isempty(MCS_and_scheduling.freq_indices))
            [LLR_SD,M] = LTE_UL_detect_CLSM(MCS_and_scheduling,rx_user_symbols,H_est,LTE_params,filtering,receiver,sigma_n2,UE_output,uu, bb);
        end
    case 5
        [LLR_SD,M] = LTE_UL_detect_MUMIMO(MCS_and_scheduling,filtering,rx_user_symbols,H_est,LTE_params,receiver,sigma_n2,UE_output, bb);
    otherwise
        error('This transmission mode is not standardized in Uplink LTE');
 end
