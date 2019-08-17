function [rank_i,PMI,CQI] =  LTE_UL_link_adaption_CLSM(UE_config,BS_config,CQI_params,CQI_mapping,SINR_averager,DFT_spreading_off,MI_data,Ntot,Ns,Nfft,Nrb,ICI_power,sigma_n2,channel_MSE,H_est_complete,H_cal,cqi_i,ue_allocation)

% if there is no CQI FB, there is no FB at all
if UE_config.CQI_fb % wheter CQI feedback is used

    % init varibales here
    H       = reshape(mean(mean(H_est_complete,1),2),size(H_est_complete,3),size(H_est_complete,4)); % this is just used to extract the MIMO dimensions
    N_TX    = size(H,2);
    N_RX    = size(H,1);
    RB_max  = size(H_est_complete,1);

    %% set the codebook size
    if N_TX == 2
        i_max = [5,0];
    else
        i_max = [23,15,11,0];
    end

    %% channel averaging
    if UE_config.channel_averaging % use channel averaging
        RB_max = size(H_est_complete,1)/12;     % average over subcarrier
        H_temp = NaN(RB_max,2*Ns,N_RX,N_TX);    

        ICI_temp = NaN(1,RB_max);               % average the ICI power
        for RB_i=1:RB_max
            freq_band           = (RB_i-1)*12+1:min(RB_i*12,size(H_est_complete,1));
            H_temp(RB_i,:,:,:)  = mean(H_est_complete(freq_band,:,:,:),1);
            ICI_temp(RB_i)      = mean(ICI_power(freq_band));
        end
        H_est_complete = H_temp;
        if UE_config.ignore_ISI_ICI
            ICI_temp = zeros(1,RB_max);
        end
        ICI_power = ICI_temp;
    end

    %% some initializations
    max_rank = min(size(H));
    I_g     = zeros(max(i_max)+1,max_rank,RB_max,2);   % total sum rate over all resource blocks
    SNR_g   = zeros(max(i_max)+1,max_rank,RB_max,2,max_rank);

    if UE_config.RI_fb % wheter RI feedback is activated
        rank_loop = 1:max_rank;
    else % iterate only over the actual RI value if RI feedback is deactivated
        rank_loop = UE_config.RI;
    end
    I       = zeros(max(i_max)+1,max_rank,RB_max);
    SNR     = zeros(max(i_max)+1,max(rank_loop),RB_max,max(rank_loop));

    for slot_i = 1:2
        H_mean  = squeeze(mean(H_est_complete(:,(slot_i-1)*Ns+1:slot_i*Ns,:,:),2)); % average over time
        H_Rx    = [];
        H_diag  = [];

        % construct diagonal channel
        for ii=1:N_TX
            for jj=1:N_RX
                H_Rx = [H_Rx; diag(squeeze(H_mean(:,jj,ii)))];
            end
            H_diag  = [H_diag,H_Rx];
            H_Rx    = [];
        end

        for rr = rank_loop
            % exhaustive search over all layer
            if UE_config.PMI_fb
                i_start = 0;
                i_stop  = i_max(rr);
                i_range = 0:max(i_max(rank_loop));
            else    % static scheduling with fixed PMI
                i_start = UE_config.PMI;
                i_stop  = UE_config.PMI;
                i_range = i_start: i_stop;
            end

            for i = i_start:i_stop   % exhaustive search over precoders
                % exhaustive search over all precoders
                W = LTE_UL_get_precoding_matrix(4,N_TX,i,rr,NaN);

                % some inits
                ISI_power   = zeros(rr,RB_max);
                E           = eye(RB_max);
                idx_h       = zeros(size(H_diag));
                H_eff_k     = zeros(size(H*W));
                idx_f       = zeros(size(W'*H'));
                F           = zeros( size(W'*H')*RB_max );

                % equalizer for effective channel
                for kk=1:RB_max
                    idx_h   = logical(kron(ones(size(H)),diag(E(:,kk))));
                    H_eff_k = reshape(H_diag( idx_h ),size(H))*W;
                    idx_f   = logical(kron(ones(size(H_eff_k')),diag(E(:,kk))));
                    switch BS_config.receiver
                        case 'ZF'
                            F_k = pinv(H_eff_k);
                        case 'MMSE'
                            F_k = pinv(H_eff_k'*H_eff_k+sigma_n2*eye(size(H_eff_k'*H_eff_k)))*H_eff_k';
                    end
                    F(idx_f) = F_k;
                end
                H_eff   = H_diag*kron(W,eye(RB_max));
                F       = reshape(F,size(H_eff'));

                % for ISI calculation
                if ~UE_config.ignore_ISI_ICI
                    for layer=1:rr
                        F_block = zeros(RB_max,RB_max,rr,N_RX);
                        for rx=1:N_RX
                            for ll=1:rr
                                selector            = zeros([rr,N_RX]);
                                selector(ll,rx)     = 1;
                                F_block(:,:,ll,rx)  = reshape(F(logical(kron(selector,ones(RB_max)))),[RB_max RB_max]);
                            end
                        end

                        % ISI calculation
                        M = [zeros((Nfft-Ntot)/2,Ntot); eye(Ntot); zeros((Nfft-Ntot)/2,Ntot)];
                        for ri=1:N_RX
                            for rj=1:N_RX
                                H_cal_sum = zeros(Nfft,Nfft);
                                for tt=1:N_TX
                                    H_cal_sum = H_cal_sum + H_cal(:,:,ri,tt)'*H_cal(:,:,rj,tt);
                                end

                                cutted_H_cal=(M'*diag(diag(H_cal_sum))*M)/Nfft;
                                % get user ISI
                                cutted_H_cal = diag(cutted_H_cal(logical(diag(ue_allocation))));

                                if UE_config.channel_averaging
                                    cutted_temp = NaN(RB_max,1);           %average over subcarrier
                                    for RB_i=1:RB_max
                                        freq_band           = (RB_i-1)*12+1:min(RB_i*12,size( cutted_H_cal,1));
                                        cutted_H_cal_diag   = diag(cutted_H_cal);
                                        cutted_temp(RB_i)   = mean(cutted_H_cal_diag(freq_band));
                                    end
                                    cutted_H_cal = diag(cutted_temp);
                                end
                                ISI_power(layer,:) = ISI_power(layer,:) + 2.*transpose(diag(F_block(:,:,layer,rj)*cutted_H_cal*F_block(:,:,layer,ri)'));
                            end

                        end
                    end
                end

                if DFT_spreading_off   %OFDMA
                    SNR_k = nan(rr,RB_max);
                    for kk=1:RB_max
                        idx_h       = logical(kron(ones(size(H)),diag(E(:,kk))));
                        H_eff_k     = reshape(H_diag(idx_h),size(H))*W  ;
                        idx_f       = logical(kron(ones(size(H_eff_k')),diag(E(:,kk))));
                        F_k         = reshape(F(idx_f),size(H_eff_k,2),[]);
                        K           = F_k*H_eff_k;
                        SNR_k(:,kk) = abs(diag(K)).^2./(sum(abs(K-diag(diag(K))).^2,2) + ISI_power(:,kk) + (channel_MSE + sigma_n2+ICI_power(kk)).*sum(abs(F_k).^2,2));
                    end
                    SNR(i+1,rr,:,1:rr) = SNR_k.';

                else % SC-FDMA
                    noise_enhancement   = zeros(rr,1);
                    signal_power        = zeros(rr,1);
                    interference_power  = zeros(rr,1);
                    for rrr=1:rr
                        S_l                                     = zeros( RB_max,rr*RB_max);
                        S_l(:,((rrr-1)*RB_max)+1:rrr*RB_max)    = eye(RB_max);
                        noise_enhancement(rrr)                  = norm(S_l*F,'fro')^2;
    %                     signal_power(rrr)                       = norm(S_l*diag(F*H_eff),2)^2;
                        signal_power(rrr)                       = (1/RB_max)*abs(sum(S_l*diag(F*H_eff)))^2;
                        interference_power(rrr)                 = norm(S_l*F*H_eff,'fro')^2-(1/RB_max)*abs(sum(S_l*diag(F*H_eff)))^2;
                    end

                    ISI_power_rep           = repmat(ISI_power(rr,:),rr,1);
                    ICI_power_rep           = repmat(ICI_power,rr,1);
                    noise_enhancement_rep   = repmat(noise_enhancement,1,RB_max);
                    signal_power_rep        = repmat(signal_power,1,RB_max);
                    interference_power_rep  = repmat(interference_power,1,RB_max);

                    SINR                = signal_power_rep./(interference_power_rep + RB_max*ISI_power_rep + (channel_MSE + sigma_n2+ICI_power_rep).*noise_enhancement_rep);
                    SNR(i+1,rr,:,1:rr)  = SINR.';
                end % if DFT_spreading_off

                SNR_temp = squeeze(SNR(i+1,rr,:,:));

                dummy_struct = [];
                dummy_struct.MI_data = MI_data;
                I(i+1,rr,:) = LTE_UL_feedback_getBICM(dummy_struct,harmmean(SNR_temp,1));

            end  %end i loop
        end  %end rr-loop

        I_g(:,:,:,slot_i) = I;

        SNR_dB = 10*log10(SNR);
        SNR_dB(isnan(SNR_dB))= -inf;
        SNR_g(i_range+1,1:max(rank_loop),:,slot_i,1:max(rank_loop)) = SNR_dB(i_range+1,1:max(rank_loop),:,1:max(rank_loop));
    end % for slot_i

    I_ges = sum(sum(I_g,3),4);
    [~,I1] = max(max(I_ges,[],1));

    %====================================

    if UE_config.RI_fb
        rank_i      = I1;    % choose rank indicator to maximize mutual information   => pre RI
        rank_start  = max(1,rank_i-1);
        rank_stop   = min(rank_i+1, max_rank);
    else
        rank_i      = UE_config.RI;
        rank_start  = UE_config.RI;
        rank_stop   = UE_config.RI;
    end


    %====================================
    LTE_rate = -Inf*ones(RB_max,2,rank_stop);
% if UE_config.CQI_fb % wheter CQI feedback is used
 
    for rank_i = rank_start:rank_stop
        
        if UE_config.PMI_fb
            [~,I2] = max(I_ges(:,rank_i));
            PMI = (I2-1); % codebook index starts at 0 not at 1
        else
            PMI = UE_config.PMI;
        end
        
        LTE_rate(:,:,rank_i) = 0;
        for i2 = 1:min(2,rank_i)   % iterate over the number of streams
            if rank_i == 4
                nLayers = 2;
                AWGN_SNR = squeeze(SNR_g(PMI(1)+1,rank_i,:,:,(i2-1)*2+1:i2*2));
            elseif rank_i == 3 && i2 == 2
                nLayers = 2;
                AWGN_SNR = squeeze(SNR_g(PMI(1)+1,rank_i,:,:,i2:rank_i));
            else
                nLayers = 1;
                AWGN_SNR = squeeze(SNR_g(PMI(1)+1,rank_i,:,:,i2));
            end
            CQIs = (0:15);
            SINReff             = SINR_averager.average(10.^(AWGN_SNR(:)/10),CQIs,[CQI_params(20).modulation_order,CQI_params(CQIs(2:end)).modulation_order]);
            SINReff(SINReff<-30)= -30;
            CQI_temp            = LTE_common_CQI_mapping_table(CQI_mapping,SINReff,CQIs+1);
            CQI(i2,rank_i)      = CQI_temp;

            if CQI_temp == 20
                LTE_rate(:,:,rank_i) = -Inf;
            else
                LTE_rate(:,:,rank_i) = LTE_rate(:,:,rank_i) + CQI_params(CQI_temp).efficiency*nLayers;
            end
        end
    end
    
    sum_LTE_rate = sum(sum(LTE_rate,1),2);
    
    %adjusted rank_i
    if UE_config.RI_fb
        %         rank_i=find(meanCQI==max(meanCQI), 1, 'last' );
        [~,rank_i] = max(sum_LTE_rate);
    else
        rank_i = UE_config.RI;
    end
    
    
    %% NOTE: now we know the rank_i --> we can select the PMI
    if UE_config.PMI_fb
        [~,I2] = max(I_ges(:,rank_i));
        PMI = (I2-1); % codebook index starts at 0
    else
        PMI = UE_config.PMI;
    end
    CQI=squeeze(CQI(:,rank_i));
    
    
else
    CQI = cqi_i*ones(UE_config.RI,1);
    rank_i = UE_config.RI;
    PMI = UE_config.PMI;
end

end

