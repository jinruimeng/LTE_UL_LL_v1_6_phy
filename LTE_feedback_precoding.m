function [rank_i,PMI] =  LTE_feedback_precoding(nAtPort,sigma_n2,LTE_params,alphabet,channel,UE,uu)
% author Stefan Schwarz
% contact stefan.schwarz@nt.tuwien.ac.at
% calculates the PMI to be used, RI feedback right now not included!

%% channel prediction
% save channel matrix for channel prediction 
UE.previous_channels = circshift(UE.previous_channels,[0,0,-1,0,0]);
UE.previous_channels(:,:,end,:,:)=channel(:,:,:,:);
H_est_complete = LTE_channel_predictor(UE.previous_channels,LTE_params.downlink_delay,LTE_params.ChanMod_config.filtering,UE.predict);

rank_i = min(size(H_est_complete,3),size(H_est_complete,4));    % NOTE: remove when rank_i feedback is supported!
H = reshape(mean(mean(H_est_complete,1),2),size(H_est_complete,3),size(H_est_complete,4)); 
max_rank = min(size(H));
PMI = zeros(LTE_params.Nrb,2);
PMI_temp = zeros(max_rank,LTE_params.Nrb,2);
I_g = zeros(max_rank,1);   % total sum rate over all resource blocks (necessary to choose the RI)

for PMI_i = 1:ceil(LTE_params.Nrb/UE.PMI_fb_gran)
    for slot_i = 1:2
        
        freq_band = (PMI_i-1)*UE.PMI_fb_gran*LTE_params.Nsc+1:min(PMI_i*UE.PMI_fb_gran*LTE_params.Nsc,size(H_est_complete,1));
        freq_band_size = (freq_band(end)-freq_band(1)+1);
        H_est = H_est_complete(freq_band,(slot_i-1)*LTE_params.Ns+1:slot_i*LTE_params.Ns,:,:);
        %% calculate number of frequency intervals needed
        const = 1.2;        % value empirically determined (make it large, e.g. 10^8, then no averaging is done)
        maxi = zeros(size(H_est,3),size(H_est,4));
        for i1 = 1:size(H_est,3)
            for i2 = 1:size(H_est,4)
                Htemp = abs(H_est(:,1,i1,i2));
                indi = logical([1;abs(diff(sign(diff(Htemp)))/2);1]);
                tempi = 1:LTE_params.Ntot;
                maxi(i1,i2) = mean(abs(diff(Htemp(indi))).'./diff(tempi(indi))*freq_band_size*const);
            end
        end
        interval_nr = ceil(max(max(maxi)));
        if interval_nr > LTE_params.Ntot
            interval_nr = LTE_params.Ntot;
        end
        if interval_nr == 0
            interval_nr = 1;
        end

        %% calculate feedback   
        switch nAtPort
            case 2
                switch UE.RIandPMI_fb
                    case true
                        rank_i = 1; % initialize RI for 2x1 
                        i_max = 3;  
                        I = zeros(i_max+2,interval_nr);
                        interval_size_h =  ceil(freq_band_size/interval_nr);
                        interval_size_l =  floor(freq_band_size/interval_nr);
                        change = (freq_band_size-interval_nr*interval_size_l);
                        for kk = 1:interval_nr
                            if kk <= change
                                H_t = reshape(mean(mean(H_est((kk-1)*interval_size_h+1:kk*interval_size_h,:,:,:),1),2),size(H_est,3),size(H_est,4));
                            else
                                temp_ind = change*interval_size_h;
                                H_t = reshape(mean(mean(H_est(temp_ind+(kk-change-1)*interval_size_l+1:temp_ind+(kk-change)*interval_size_l,:,:,:),1),2),size(H_est,3),size(H_est,4));
                            end
                            for i = 0:i_max
                              [Z W U D] = LTE_common_get_precoding_matrix(4,nAtPort,i,1,LTE_params);
                              I(i+1,kk) = I(i+1,kk)+log2(det(eye(1)+1/sigma_n2*W'*(H_t'*H_t)*W));
                            end
                            if size(H,1) == 2   % if there are 2 receive antennas ( NOTE: this will be necessary when the rank_i is also evaluated!)
                                [Z W U D] = LTE_common_get_precoding_matrix(4,nAtPort,0,2,LTE_params);
                                I(i+2,kk) = I(i+2,kk)+abs(log2(det(eye(2)+1/sigma_n2*W'*(H_t'*H_t)*W)));
                            end
                        end
                        I_tot = sum(I(:,1:change),2)*interval_size_h+sum(I(:,change+1:end),2)*interval_size_l;  % mutual information for 2x2 (the same for all precoding matrices)
                        [C,PMI((PMI_i-1)*UE.PMI_fb_gran+1:min(PMI_i*UE.PMI_fb_gran,LTE_params.Nrb),slot_i)] = max(I_tot(:));
                        I_g = I_g + [max(I_tot(1:i_max+1));I_tot(i_max+2)];
                        if PMI_i == ceil(LTE_params.Nrb/UE.PMI_fb_gran) && slot_i == 2 && I_g(2) > I_g(1)   % this is just a temporary solution to get the feedback to run! 
                            rank_i = 2;
                            i_max = 3;
                            BEP4 = zeros(i_max,1);
                            BEP16 = zeros(i_max,1);
                            BEP64 = zeros(i_max,1);
                            switch LTE_params.UE_config.receiver
                                case 'ZF'
                                    for i = 1:i_max-1
                                        [Z W U D] = LTE_common_get_precoding_matrix(4,nAtPort,i,rank_i,LTE_params);
                                        hk(:,:) = pinv(H*W);
                                        SNR(1) = 1/(sigma_n2*sum(abs(hk(1,:)).^2));
                                        SNR(2) = 1/(sigma_n2*sum(abs(hk(2,:)).^2));
                                        % 4QAM
                                        BEP4(i+1) = mean(2*qfunc(sqrt(SNR))-qfunc(sqrt(SNR)).^2); % NOTE: which biterror probability has to be taken depends on the CQI
                                        % 16QAM approximation (for high SNR)
                                        BEP16(i+1) = mean(3/4*qfunc(sqrt(1/5*SNR)));              % also NOTE: the two independent streams may have different modulation order
                                                                                                  % which would require to take all possible combinations into account
                                        % 64QAM approximation
                                        BEP64(i+1) = mean(7/23*qfunc(sqrt(3/63*SNR)));
                                    end
                                case 'MMSE'
                                    for i = 1:i_max-1
                                        [Z W U D] = LTE_common_get_precoding_matrix(4,nAtPort,i,rank_i,LTE_params);
                                        hk(:,:) = H*W;
                                        SINR(1) = abs(hk(:,1)'/(sigma_n2*eye(2)+hk(:,2)*hk(:,2)')*hk(:,1));
                                        SINR(2) = abs(hk(:,2)'/(sigma_n2*eye(2)+hk(:,1)*hk(:,1)')*hk(:,2));
                                        % 4QAM
                                        BEP4(i+1) = mean(2*qfunc(sqrt(SINR))-qfunc(sqrt(SINR)).^2); 
                                        % 16QAM approximation
                                        BEP16(i+1) = mean(3/4*qfunc(sqrt(1/5*SINR)));
                                        % 64QAM approximation
                                        BEP64(i+1) = mean(7/23*qfunc(sqrt(3/63*SINR)));
                                    end
                                case {'MMSE_SIC','SSD'}
                                    for i = 1:i_max-1   
                                        clear SINR
                                        [Z W U D] = LTE_common_get_precoding_matrix(4,nAtPort,i,rank_i,LTE_params);
                                        hk(:,:) = H*W;
                                        SINR(1) = abs(hk(:,1)'/(sigma_n2*eye(2)+hk(:,2)*hk(:,2)')*hk(:,1));
                                        SINR(2) = abs(hk(:,2)'/(sigma_n2*eye(2)+hk(:,1)*hk(:,1)')*hk(:,2));
                                        [C,I]=min(SINR);
                                        [C,I2]=max(SINR);
                                        SINR(3) = sum(abs(hk(:,I)).^2)/sigma_n2;
                                        SINR(I) = 1/(sigma_n2*sum(abs(pinv(hk(:,I))).^2)+abs((pinv(hk(:,I))*hk(:,I2)))^2*2);
                                        % 4QAM
                                        BEP_tmp4 = 2*qfunc(sqrt(SINR))-qfunc(sqrt(SINR)).^2;
                                        % 16QAM approximation
                                        BEP_tmp16 = 3/4*qfunc(sqrt(1/5*SINR));
                                        % 64QAM approximation
                                        BEP_tmp64 = 7/23*qfunc(sqrt(3/63*SINR));
                                        BEP4(i+1) = 1/2*(BEP_tmp4(I2)+(1-BEP_tmp4(I2))*BEP_tmp4(3)+BEP_tmp4(I2)*BEP_tmp4(I)); 
                                        BEP16(i+1) = 1/2*(BEP_tmp16(I2)+(1-BEP_tmp16(I2))*BEP_tmp16(3)+BEP_tmp16(I2)*BEP_tmp16(I)); 
                                        BEP64(i+1) = 1/2*(BEP_tmp64(I2)+(1-BEP_tmp64(I2))*BEP_tmp64(3)+BEP_tmp64(I2)*BEP_tmp64(I)); 
                                    end
                            end
                            switch alphabet     % NOTE: one should switch here between the QAM constellation needed depending on the BEP 
                                case 2
                                    [C,PMI]=min(BEP4); 
                                case 4
                                    [C,PMI]=min(BEP16); 
                                case 6
                                    [C,PMI]=min(BEP64); 
                            end
                            PMI = PMI*ones(LTE_params.Nrb,2);
                        end
                    case false
                        switch LTE_params.scheduler.nLayers(uu)
                            case 1
                                rank_i = 1; % initialize RI for 2x1 
                                i_max = 3;  
                                I = zeros(i_max+2,interval_nr);
                                interval_size_h =  ceil(freq_band_size/interval_nr);
                                interval_size_l =  floor(freq_band_size/interval_nr);
                                change = (freq_band_size-interval_nr*interval_size_l);
                                for kk = 1:interval_nr
                                    if kk <= change
                                        H_t = reshape(mean(mean(H_est((kk-1)*interval_size_h+1:kk*interval_size_h,:,:,:),1),2),size(H_est,3),size(H_est,4));
                                    else
                                        temp_ind = change*interval_size_h;
                                        H_t = reshape(mean(mean(H_est(temp_ind+(kk-change-1)*interval_size_l+1:temp_ind+(kk-change)*interval_size_l,:,:,:),1),2),size(H_est,3),size(H_est,4));
                                    end
                                    for i = 0:i_max
                                      [Z W U D] = LTE_common_get_precoding_matrix(4,nAtPort,i,1,LTE_params);
                                      I(i+1,kk) = I(i+1,kk)+log2(det(eye(1)+1/sigma_n2*W'*(H_t'*H_t)*W));
                                    end
                                end
                                I_tot = sum(I(:,1:change),2)*interval_size_h+sum(I(:,change+1:end),2)*interval_size_l;  % mutual information for 2x2 (the same for all precoding matrices)
                                [C,PMI((PMI_i-1)*UE.PMI_fb_gran+1:min(PMI_i*UE.PMI_fb_gran,LTE_params.Nrb),slot_i)] = max(I_tot(:));
                            case 2
                                rank_i = 2;
                                i_max = 3;
                                BEP4 = zeros(i_max,1);
                                BEP16 = zeros(i_max,1);
                                BEP64 = zeros(i_max,1);
                                switch LTE_params.UE_config.receiver
                                    case 'ZF'
                                        for i = 1:i_max-1
                                            [Z W U D] = LTE_common_get_precoding_matrix(4,nAtPort,i,rank_i,LTE_params);
                                            hk(:,:) = pinv(H*W);
                                            SNR(1) = 1/(sigma_n2*sum(abs(hk(1,:)).^2));
                                            SNR(2) = 1/(sigma_n2*sum(abs(hk(2,:)).^2));
                                            % 4QAM
                                            BEP4(i+1) = mean(2*qfunc(sqrt(SNR))-qfunc(sqrt(SNR)).^2); % NOTE: which biterror probability has to be taken depends on the CQI
                                            % 16QAM approximation (for high SNR)
                                            BEP16(i+1) = mean(3/4*qfunc(sqrt(1/5*SNR)));              % also NOTE: the two independent streams may have different modulation order
                                                                                                      % which would require to take all possible combinations into account
                                            % 64QAM approximation
                                            BEP64(i+1) = mean(7/23*qfunc(sqrt(3/63*SNR)));
                                        end
                                    case 'MMSE'
                                        for i = 1:i_max-1
                                            [Z W U D] = LTE_common_get_precoding_matrix(4,nAtPort,i,rank_i,LTE_params);
                                            hk(:,:) = H*W;
                                            SINR(1) = abs(hk(:,1)'/(sigma_n2*eye(2)+hk(:,2)*hk(:,2)')*hk(:,1));
                                            SINR(2) = abs(hk(:,2)'/(sigma_n2*eye(2)+hk(:,1)*hk(:,1)')*hk(:,2));
                                            % 4QAM
                                            BEP4(i+1) = mean(2*qfunc(sqrt(SINR))-qfunc(sqrt(SINR)).^2); 
                                            % 16QAM approximation
                                            BEP16(i+1) = mean(3/4*qfunc(sqrt(1/5*SINR)));
                                            % 64QAM approximation
                                            BEP64(i+1) = mean(7/23*qfunc(sqrt(3/63*SINR)));
                                        end
                                    case {'MMSE_SIC','SSD'}
                                        for i = 1:i_max-1   
                                            clear SINR
                                            [Z W U D] = LTE_common_get_precoding_matrix(4,nAtPort,i,rank_i,LTE_params);
                                            hk(:,:) = H*W;
                                            SINR(1) = abs(hk(:,1)'/(sigma_n2*eye(2)+hk(:,2)*hk(:,2)')*hk(:,1));
                                            SINR(2) = abs(hk(:,2)'/(sigma_n2*eye(2)+hk(:,1)*hk(:,1)')*hk(:,2));
                                            [C,I]=min(SINR);
                                            [C,I2]=max(SINR);
                                            SINR(3) = sum(abs(hk(:,I)).^2)/sigma_n2;
                                            SINR(I) = 1/(sigma_n2*sum(abs(pinv(hk(:,I))).^2)+abs((pinv(hk(:,I))*hk(:,I2)))^2*2);

                                            % 4QAM
                                            BEP_tmp4 = 2*qfunc(sqrt(SINR))-qfunc(sqrt(SINR)).^2;
                                            % 16QAM approximation
                                            BEP_tmp16 = 3/4*qfunc(sqrt(1/5*SINR));
                                            % 64QAM approximation
                                            BEP_tmp64 = 7/23*qfunc(sqrt(3/63*SINR));
                                            BEP4(i+1) = 1/2*(BEP_tmp4(I2)+(1-BEP_tmp4(I2))*BEP_tmp4(3)+BEP_tmp4(I2)*BEP_tmp4(I)); 
                                            BEP16(i+1) = 1/2*(BEP_tmp16(I2)+(1-BEP_tmp16(I2))*BEP_tmp16(3)+BEP_tmp16(I2)*BEP_tmp16(I)); 
                                            BEP64(i+1) = 1/2*(BEP_tmp64(I2)+(1-BEP_tmp64(I2))*BEP_tmp64(3)+BEP_tmp64(I2)*BEP_tmp64(I)); 
                                        end
                                end
                                switch alphabet     % NOTE: one should switch here between the QAM constellation needed depending on the BEP 
                                    case 2
                                        [C,PMI]=min(BEP4); 
                                    case 4
                                        [C,PMI]=min(BEP16); 
                                    case 6
                                        [C,PMI]=min(BEP64); 
                                end
                                PMI = PMI*ones(LTE_params.Nrb,2);
                        end
                end
            case 4
                i_max = 15;
                I = zeros(max_rank,i_max+1,interval_nr);
                interval_size_h =  ceil(freq_band_size/interval_nr);
                interval_size_l =  floor(freq_band_size/interval_nr);
                change = (freq_band_size-interval_nr*interval_size_l);
                if UE.RIandPMI_fb % wheter RI and PMI feedback is activated
                    rank_loop = 1:max_rank;
                else    % iterate only over the actual RI value if no RI feedback is activated         
                    rank_loop = LTE_params.scheduler.nLayers(uu);
                end
                for rr = rank_loop    % Iterate over all possible layer numbers, to find the one that delivers the maximum sum rate
                    for kk = 1:interval_nr
                        if kk <= change
                            H = reshape(mean(mean(H_est((kk-1)*interval_size_h+1:kk*interval_size_h,:,:,:),1),2),size(H_est,3),size(H_est,4));
                        else
                            temp_ind = change*interval_size_h;
                            H = reshape(mean(mean(H_est(temp_ind+(kk-change-1)*interval_size_l+1:temp_ind+(kk-change)*interval_size_l,:,:,:),1),2),size(H_est,3),size(H_est,4));
                        end
                        for i = 0:i_max
                          [Z W U D] = LTE_common_get_precoding_matrix(4,nAtPort,i,rr,LTE_params);
                          I(rr,i+1,kk) = I(rr,i+1,kk)+log2(det(eye(rr)+1/sigma_n2*W'*(H'*H)*W));
                        end
                    end
                end
                I_tot = sum(I(:,:,1:change),3)*interval_size_h+sum(I(:,:,change+1:end),3)*interval_size_l;
                I_g = I_g + max(I_tot,[],2);
                ind_tmp = (PMI_i-1)*UE.PMI_fb_gran+1:min(PMI_i*UE.PMI_fb_gran,LTE_params.Nrb);
                [C,I_tmp] = max(I_tot,[],2);
                PMI_temp(:,ind_tmp,slot_i) = repmat(I_tmp,1,numel(ind_tmp));
                if PMI_i == ceil(LTE_params.Nrb/UE.PMI_fb_gran) && slot_i == 2  % end of PMI and slot loop reached --> calculate feedback values
                    [C,rank_i] = max(I_g);
                    PMI = squeeze(PMI_temp(rank_i,:,:));
                end
        end
    end
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         