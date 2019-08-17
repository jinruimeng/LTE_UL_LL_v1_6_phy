function [rank_i,PMI,CQI] =  LTE_UL_feedback_CLSM(nAtPort,sigma_n2,LTE_params,channel,UE,uu,H_cal,cqi_i)

channel=channel.genie.H_fft;

%% channel prediction
% save channel matrix for channel prediction
UE.previous_channels = circshift(UE.previous_channels,[0,0,-1,0,0]);
UE.previous_channels(:,:,end,:,:)=channel(:,:,:,:);
H_est_complete = LTE_UL_channel_predictor(UE.previous_channels,LTE_params.downlink_delay,LTE_params.ChanMod_config.filtering,UE.predict);
H = reshape(mean(mean(H_est_complete,1),2),size(H_est_complete,3),size(H_est_complete,4));

N_TX = size(H,2);
N_RX = size(H,1);


%% set the codebook size for the exhaustive search
if nAtPort == 2
    i_max = [5,0];
else
    i_max = [23,15,11,0];
end

%% channel averaging
if LTE_params.UE_config.channel_averaging % use channel averaging
    RB_max = LTE_params.Nrb;             % leads to a avaraging of 7 RE (1 subcarrier 7 OFDM symbols)
    
    H_temp=nan(RB_max,2*LTE_params.Ns,N_RX,N_TX);           %average over subcarrier
    
    
    ICI_temp=nan(1,RB_max);             %average the ICI power
    for RB_i=1:RB_max
        freq_band = (RB_i-1)*LTE_params.Nsc+1:min(RB_i*LTE_params.Nsc,size(H_est_complete,1));
        %         freq_band = freq_band(ceil(length(freq_band)/2));
        H_temp(RB_i,:,:,:)=mean(H_est_complete(freq_band,:,:,:),1);
        ICI_temp(RB_i)=mean(LTE_params.ICI_power(freq_band));
    end
    H_est_complete=H_temp;
    if LTE_params.UE_config.ignore_ISI_ICI
        ICI_temp = zeros(1,RB_max);
    end
    ICI_power=ICI_temp;
    
else
    RB_max = LTE_params.Ntot;            % calculations are performed for each subcarrier separately
    ICI_power = LTE_params.ICI_power;
end

%% some initializations
max_rank = min(size(H));
I_g = zeros(max(i_max)+1,max_rank,RB_max,2);   % total sum rate over all resource blocks
SNR_g = zeros(max(i_max)+1,max_rank,RB_max,2,max_rank);

if UE.RI_fb % wheter RI feedback is activated
    rank_loop = 1:max_rank;
else    % iterate only over the actual RI value if RI feedback is deactivated
    rank_loop = LTE_params.UE_config.RI;
end

I = zeros(max(i_max)+1,max_rank,RB_max);
SNR = zeros(max(i_max)+1,max(rank_loop),RB_max,max(rank_loop));


for slot_i = 1:2
    
    H_mean=squeeze(mean(H_est_complete(:,(slot_i-1)*LTE_params.Ns+1:slot_i*LTE_params.Ns,:,: )  ,2));        %average over time
    
    H_Rx=[];
    H_diag=[];
    %produce diagonal channel
    for ii=1:N_TX
        for jj=1:N_RX
            H_Rx =[H_Rx; diag(squeeze(H_mean(:,jj,ii)))];
        end
        H_diag=[H_diag,H_Rx];
        H_Rx=[];
    end
    
    
    for rr = rank_loop
        
        % exhaustive search over all layer
        if UE.PMI_fb
            i_start = 0;
            i_stop = i_max(rr);
            
            i_range=0:max(i_max(rank_loop));
        else  %static scheduling with fixed PMI
            i_start = LTE_params.UE_config.PMI;
            i_stop = LTE_params.UE_config.PMI;
            i_range=i_start: i_stop;
        end
        
        
        
        for i = i_start:i_stop   %exhaustive search over precoder!
            
            % exhaustive search over all precoders
            W = LTE_UL_get_precoding_matrix(4,nAtPort,i,rr,LTE_params);
            
            % some inits
            ISI_power = zeros(rr,RB_max);
            E=eye(RB_max);
            idx_h=zeros(size(H_diag));
            H_eff_k=zeros(size(H*W));
            idx_f =zeros(size(W'*H'));
            F=zeros( size(W'*H')*RB_max );
            
            
            % equalizer for effective channel
            for kk=1:RB_max
                idx_h=logical(kron(ones(size(H)), diag(E(:,kk))  ));
                H_eff_k=reshape(H_diag( idx_h ),size(H))    *W  ;
                idx_f = logical(kron(ones(size(H_eff_k')),diag(E(:,kk)) ));
                switch LTE_params.BS_config.receiver
                    case 'ZF'
                        F_k=pinv(H_eff_k);
                    case 'MMSE'
                        F_k = pinv(H_eff_k'*H_eff_k+sigma_n2*eye(size(H_eff_k'*H_eff_k)))*H_eff_k';
                end
                F( idx_f  ) = F_k;
            end
            H_eff=H_diag*kron(W ,eye(RB_max));
            F=reshape(F,size(H_eff'));
            
            % for ISI calculation
            if ~LTE_params.UE_config.ignore_ISI_ICI
                for layer=1:rr
                    F_block=zeros(RB_max,RB_max,rr,N_RX);
                    for rx=1:N_RX
                        for ll=1:rr
                            selector = zeros([rr,N_RX]);
                            selector(ll,rx) = 1;
                            F_block(:,:,ll,rx) = reshape(F(logical(kron(selector,ones(RB_max)))),[RB_max RB_max]);
                        end
                    end
                    
                    % ISI calculation
                    
                    M = [zeros((LTE_params.Nfft-LTE_params.Ntot)/2,LTE_params.Ntot); eye(LTE_params.Ntot); zeros((LTE_params.Nfft-LTE_params.Ntot)/2,LTE_params.Ntot)];
                    for ri=1:N_RX
                        for rj=1:N_RX
                            H_cal_sum = zeros(LTE_params.Nfft,LTE_params.Nfft);
                            for tt=1:N_TX
                                H_cal_sum = H_cal_sum +  H_cal(:,:,ri,tt)'*H_cal(:,:,rj,tt) ;
                            end
                            
                            cutted_H_cal=(M'*diag(diag(H_cal_sum))*M)/LTE_params.Nfft;
                            
                            if LTE_params.UE_config.channel_averaging
                                cutted_temp=nan(RB_max,1);           %average over subcarrier
                                for RB_i=1:RB_max
                                    freq_band = (RB_i-1)*LTE_params.Nsc+1:min(RB_i*LTE_params.Nsc,size( cutted_H_cal,1));
                                    cutted_H_cal_diag=diag(cutted_H_cal);
                                    cutted_temp(RB_i)=mean(cutted_H_cal_diag(freq_band));
                                end
                                cutted_H_cal=diag(cutted_temp);
                            end
                            ISI_power(layer,:) = ISI_power(layer,:) + 2.* transpose(diag(F_block(:,:,layer,rj)*cutted_H_cal*F_block(:,:,layer,ri)'));
                        end
                        
                    end
                end
            end
            
            
            
            
            
            if LTE_params.DFT_spreading_off   %OFDMA
                SNR_k=nan(rr,RB_max);
                for kk=1:RB_max
                    idx_h=logical(kron(ones(size(H)), diag(E(:,kk))  ));
                    H_eff_k=reshape(H_diag( idx_h ),size(H))    *W  ;
                    idx_f = logical(kron(ones(size(H_eff_k')),diag(E(:,kk)) ));
                    F_k = reshape(F( idx_f  ),size(H_eff_k,2),[]);
                    K= F_k*H_eff_k;
                    SNR_k(:,kk) = abs(diag(K)).^2./(sum(abs(K-diag(diag(K))).^2,2)+ISI_power(:,kk)+(sigma_n2+ICI_power(kk)).*sum(abs(F_k).^2,2));
                end
                
                SNR(i+1,rr,:,1:rr)=SNR_k.';
                
            else % SC-FDMA
                
                noise_enhancement = zeros(rr,1);
                signal_power = zeros(rr,1);
                interference_power = zeros(rr,1);
                for rrr=1:rr
                    S_l = zeros( RB_max, rr*RB_max);
                    S_l(:,((rrr-1)*RB_max)+1:rrr*RB_max) = eye(RB_max);
                    noise_enhancement(rrr) = norm(S_l*F,'fro')^2;
                    %signal_power(rrr) = norm(S_l*diag(F*H_eff),2)^2;
                    signal_power(rrr) =(1/RB_max)*abs(sum(S_l*diag(F*H_eff)))^2;
                    interference_power(rrr) = norm(S_l*F*H_eff,'fro')^2-(1/RB_max)*abs(sum(S_l*diag(F*H_eff)))^2;
                end
                
                ISI_power_rep = repmat(ISI_power(rr,:),rr,1);
                ICI_power_rep = repmat(ICI_power,rr,1);
                noise_enhancement_rep = repmat(noise_enhancement,1,RB_max);
                signal_power_rep = repmat(signal_power,1,RB_max);
                interference_power_rep = repmat(interference_power,1,RB_max);
                
                SINR=signal_power_rep./(interference_power_rep +RB_max*ISI_power_rep+ (sigma_n2+ICI_power_rep).* noise_enhancement_rep);
                
                SNR(i+1,rr,:,1:rr)=SINR.';
            end
            
            %             for kk=1:RB_max
            %                 I(i+1,rr,kk) = LTE_UL_feedback_getBICM(LTE_params,squeeze(SNR(i+1,rr,kk,:)));
            %             end
            %% NOTE: actually this code is more appropriate for SC-FDMA (at least with ZF)
            SNR_temp = squeeze(SNR(i+1,rr,:,:));
            %             I(i+1,rr,:) = LTE_UL_feedback_getBICM(LTE_params,1./(mean(1./SNR_temp,1)));
            
            %ezoechma: harmonic mean if ISI or ICI occured
            I(i+1,rr,:) = LTE_UL_feedback_getBICM(LTE_params,harmmean(SNR_temp,1));
            
            
            
        end  %end i loop
    end  %end rr-loop
    
    
    %        for slot_i = 1:2
    I_g(:,:,:,slot_i) =I;
    
    SNR_dB = 10*log10(SNR);
    SNR_dB(isnan(SNR_dB))= -inf;
    SNR_g(i_range+1,1:max(rank_loop),:,slot_i,1:max(rank_loop)) = SNR_dB(i_range+1,1:max(rank_loop),:,1:max(rank_loop));
end


I_ges = sum(sum(I_g,3),4);
[~,I1] = max(max(I_ges,[],1));

%====================================

if UE.RI_fb
    rank_i = I1;    % choose rank indicator to maximize mutual information   => pre RI
    rank_start = max(1,rank_i-1);
    rank_stop  = min(rank_i+1, max_rank);
else
    rank_i = LTE_params.UE_config.RI;
    rank_start = LTE_params.UE_config.RI;
    rank_stop = LTE_params.UE_config.RI;
end


%====================================
LTE_rate = -Inf*ones(RB_max,2,rank_stop);
if LTE_params.UE_config.CQI_fb % wheter CQI feedback is used
    
    for rank_i = rank_start:rank_stop
        
        if UE.PMI_fb
            [~,I2] = max(I_ges(:,rank_i));
            PMI = (I2-1)*ones(LTE_params.Nrb,2); % codebook index starts at 0 not at 1
        else
            PMI = (LTE_params.UE_config.PMI).*ones(LTE_params.Nrb,2);
        end
        
        
        LTE_rate(:,:,rank_i) = 0;
        for i2 = 1:min(2,rank_i)   % iterate over the number of streams
            if UE.CQI_fb_gran ~= 1  % single CQI value for whole bandwidth
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
                SINReff = UE.SINR_averager.average(10.^(AWGN_SNR(:)/10),CQIs,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(CQIs(2:end)).modulation_order]);
                SINReff(SINReff<-30)=-30;
                CQI_temp = LTE_common_CQI_mapping_table(LTE_params.CQI_mapping,SINReff,CQIs+1);
                CQI(:,:,i2,rank_i) = CQI_temp*ones(LTE_params.Nrb,2);
                
                if CQI_temp == 20
                    LTE_rate(:,:,rank_i) = -Inf;
                else
                    LTE_rate(:,:,rank_i) = LTE_rate(:,:,rank_i) + LTE_params.CQI_params(CQI_temp).efficiency*nLayers;
                end
                
            else  % single CQI value for every resource block
                for RB_ii = 1:RB_max
                    for slot_ii = 1:2
                        if rank_i == 4
                            %                             AWGN_SNR = squeeze(SNR_g(PMI(1)+1,i2,RB_ii,slot_ii,(i2-1)*2+1:min(i2*2,rank_i)));
                            %                             AWGN_SNR = squeeze(SNR_g(PMI(1)+1,i2,RB_ii,slot_ii,(i2-1)*2+1:min(i2*2,rank_i)));
                            nLayers = 2;
                            AWGN_SNR = squeeze(SNR_g(PMI(1)+1,rank_i,:,:,(i2-1)*2+1:i2*2));
                        elseif rank_i == 3 && i2 == 2
                            nLayers = 2;
                            AWGN_SNR = squeeze(SNR_g(PMI(1)+1,rank_i,RB_ii,slot_ii,i2:rank_i));
                        else
                            nLayers = 1;
                            AWGN_SNR = squeeze(SNR_g(PMI(1)+1,rank_i,RB_ii,slot_ii,i2));
                        end
                        CQIs = (0:15);
                        SINReff = UE.SINR_averager.average(10.^(AWGN_SNR(:)/10),CQIs,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(CQIs(2:end)).modulation_order]);
                        SINReff(SINReff<-30)=-30;
                        CQI_temp = LTE_common_CQI_mapping_table(LTE_params.CQI_mapping,SINReff,CQIs+1);
                        CQI(RB_ii,slot_ii,i2,rank_i) = CQI_temp;
                        
                        if CQI_temp == 20
                            LTE_rate(:,:,rank_i) = -Inf;
                        else
                            LTE_rate(RB_ii,slot_ii,rank_i) = LTE_rate(RB_ii,slot_ii,rank_i) + LTE_params.CQI_params(CQI_temp).efficiency*nLayers;
                        end
                    end
                end
                
                
                
                if ~LTE_params.UE_config.channel_averaging %Average of PMI and CQI now, to have right feedback size
                    
                    CQI_mean_size=size(CQI);
                    CQI_mean_size(1)=LTE_params.Nrb;
                    
                    CQI_mean=zeros(CQI_mean_size);
                    PMI_mean=zeros(LTE_params.Nrb,2);
                    
                    gran=LTE_params.Ntot/LTE_params.Nrb;        %number of CQI to be avaraged
                    
                    for i=0:(LTE_params.Nrb-1)
                        CQI_mean(i+1,:,:)=round(mean(CQI(i*gran+1:(i+1)*gran,:,:)));
                        PMI_mean(i+1,:)=round(mean(PMI(i*gran+1:(i+1)*gran,:,:)));
                    end
                    
                    CQI=CQI_mean;
                end
                
            end
        end
    end
    
    sum_LTE_rate = sum(sum(LTE_rate,1),2);
    
    %adjusted rank_i
    if LTE_params.UE_config.RI_fb
        %         rank_i=find(meanCQI==max(meanCQI), 1, 'last' );
        [~,rank_i] = max(sum_LTE_rate);
    else
        rank_i = LTE_params.UE_config.RI;
    end
    
    
    %% NOTE: now we know the rank_i --> we can select the PMI
    if UE.PMI_fb
        [~,I2] = max(I_ges(:,rank_i));
        PMI = (I2-1)*ones(LTE_params.Nrb,2); % codebook index starts at 0 not at 1
    else
        PMI = (LTE_params.UE_config.PMI).*ones(LTE_params.Nrb,2);
    end
    
    CQI=squeeze(CQI(:,:,:,rank_i));
    
    
else
    CQI = cqi_i*ones(LTE_params.Nrb,2,min(2,rank_i),rank_i);
end



end

