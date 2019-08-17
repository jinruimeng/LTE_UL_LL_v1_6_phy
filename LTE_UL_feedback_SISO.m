function [CQI,CQI_bar] = LTE_UL_feedback_SISO(sigma_n2,channel,UE,LTE_params,H_cal,BS,ICI_power)
% author Stefan Schwarz
% contact stefan.schwarz@nt.tuwien.ac.at
% calculates the CQI feedback for SISO


H_est_complete = channel.genie.H_fft;

ISI_power = 2.*diag( squeeze(H_cal(:,:,1,1))'*squeeze(H_cal(:,:,1,1)) ) / LTE_params.Nfft;

M=[zeros((LTE_params.Nfft-LTE_params.Ntot)/2,LTE_params.Ntot);eye(LTE_params.Ntot);zeros((LTE_params.Nfft-LTE_params.Ntot)/2,LTE_params.Ntot)];
ISI_power=M'*ISI_power;

CQI_bar = [];
if LTE_params.UE_config.channel_averaging % use channel averaging
    RB_max = LTE_params.Nrb;
    SNR = zeros(RB_max,2);
    
    ISI_temp=nan(RB_max,1);             %average the ISI power
    ICI_temp=nan(RB_max,1);             %average the ICI power
    for RB_i=1:RB_max
        freq_band = (RB_i-1)*LTE_params.Nsc+1:min(RB_i*LTE_params.Nsc,size(H_est_complete,1));
        ISI_temp(RB_i)=mean(ISI_power(freq_band));
        ICI_temp(RB_i)=mean(LTE_params.ICI_power(freq_band));
    end
    ISI_power=ISI_temp;
    ICI_power=ICI_temp;
else
    RB_max = LTE_params.Ntot;
    ICI_power = LTE_params.ICI_power.';
    SNR = zeros(RB_max,2);
end

for slot_ind=1:2
    
    H_mean=mean(H_est_complete(:,(slot_ind-1)*LTE_params.Ns+1:slot_ind*LTE_params.Ns )  ,2);        %average over time
    if LTE_params.UE_config.channel_averaging
        H_temp=nan(RB_max,1);           %average over subcarrier
        for RB_i=1:RB_max
            freq_band = (RB_i-1)*LTE_params.Nsc+1:min(RB_i*LTE_params.Nsc,size(H_est_complete,1));
            H_temp(RB_i)=mean(H_mean(freq_band));
        end
        H_mean=H_temp;
    end
    
    switch BS.receiver
        case 'ZF'
            F_n = pinv(diag(H_mean));
        case 'MMSE'
            F_n = diag(conj(H_mean)./sum(H_mean.*conj(H_mean)+sigma_n2,2));
    end
    
    if LTE_params.DFT_spreading_off   %OFDMA SNR formula
        SNR_k=((F_n*H_mean).^2)./ ((ISI_power+ICI_power+sigma_n2).*abs(diag(F_n)).^2);
    else   %SC-FDMA formula
        %length(H_mean) corresponds to LTE_params.Ntot for the channel averaging case
        
        %SNR_k = (norm(F_n*H_mean)^2) ./ ( norm(F_n*H_mean)^2 -(1/RB_max)*abs(sum(     F_n*H_mean)     )^2+           (F_n'*F_n)*ISI_power+(sigma_n2+ICI_power)*norm(diag(F_n))^2);
        SNR_k = (1/RB_max)*abs(sum(     F_n*H_mean)     )^2 ./ ( norm(F_n*H_mean)^2 -(1/RB_max)*abs(sum(     F_n*H_mean)     )^2+           (F_n'*F_n)*ISI_power+(sigma_n2+ICI_power)*norm(diag(F_n))^2);
    end
    
    SNR(:,slot_ind)=SNR_k;
    
end


CQI = zeros(LTE_params.Nrb,2);
if LTE_params.UE_config.CQI_fb
    if UE.CQI_fb_gran ~= 1  % single CQI value for whole bandwidth
        
        CQIs = 0:15;
        SINReff = UE.SINR_averager.average(SNR(:),CQIs,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(CQIs(2:end)).modulation_order]);
        CQI_temp = LTE_common_CQI_mapping_table(LTE_params.CQI_mapping,SINReff(1),CQIs+1);
        CQI = CQI_temp*ones(LTE_params.Nrb,2);
        CQI_bar = CQI_temp;
        
    else
        for RB_ii = 1:LTE_params.Nrb
            if LTE_params.UE_config.channel_averaging
                rb_ii = RB_ii;
            else
                rb_ii = (RB_ii-1)*LTE_params.Nsc+1:RB_ii*LTE_params.Nsc;
            end
            for slot_ii = 1:2
                CQIs = (0:15);
                SNR_tmp = SNR(rb_ii,slot_ii);
                SINReff = UE.SINR_averager.average(SNR_tmp,CQIs,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(CQIs(2:end)).modulation_order]);
                CQI_temp = LTE_common_CQI_mapping_table(LTE_params.CQI_mapping,SINReff,CQIs+1);
                CQI(RB_ii,slot_ii) = CQI_temp;
            end
        end
    end
else
    CQI = LTE_params.cqi_i.*ones(size(CQI));
end

end
