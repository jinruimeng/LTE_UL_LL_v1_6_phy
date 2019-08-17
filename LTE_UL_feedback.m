function [rank_i,PMI,CQI] =  LTE_UL_feedback(nAtPort,sigma_n2,LTE_params,channel,UE,uu,modus,cqi_i,H_cal,BS)
% author Stefan Schwarz
% contact stefan.schwarz@nt.tuwien.ac.at
% calculates the PMI, RI and CQI feedback

switch modus    % transmission mode
    case 1      % Single-antenna port
        if UE.CQI_fb 
            [CQI] = LTE_UL_feedback_SISO(sigma_n2,channel,UE,LTE_params,H_cal,BS);
        else
            CQI = cqi_i*ones(LTE_params.Nrb,2);
        end
        rank_i = [];
        PMI = [];
    case {2,3}     
        error('Transmission modes 2 and 3 are not standardized in Uplink LTE');
    case 4      % Closed-loop spatial multiplexing
        if UE.CQI_fb || UE.RI_fb || UE.PMI_fb
            [rank_i,PMI,CQI] = LTE_UL_feedback_CLSM(nAtPort,sigma_n2,LTE_params,channel,UE,uu,H_cal,cqi_i);
        else    
            rank_i = LTE_params.UE_config.RI;
            PMI = (LTE_params.UE_config.PMI).*ones(LTE_params.Nrb,2);
            CQI = cqi_i*ones(LTE_params.Nrb,2,min(2,rank_i),rank_i);
        end
end

