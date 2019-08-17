function [LLR_SD1,M1] = LTE_UL_detect_MUMIMO(MCS_and_scheduling,filtering,user_symbols,H_est_all,LTE_params,receiver,sigma_n2,UE_output, bb)
% MU MIMO detection
%
% (c) 2016 by ITC
% www.nt.tuwien.ac.at


nRX             = LTE_params.BS_config.nRX;

N_user = length(MCS_and_scheduling);
H_all_noref = cell(1);
% UEs_scheduled = cell(1);
% nLayers = cell(1);
% M = cell(1);
% bittable = cell(1);
% symbol_alphabet = cell(1);
% CH_mapping_UE_spec = cell(1);
% N_subc = cell(1);


for u_ = 1:N_user
    
    if max(max(MCS_and_scheduling{u_}.UE_mapping)) % if user is scheduled
        uG = utils.localToGlobalUser(LTE_params.connection_table, bb, u_);
        UEs_scheduled(u_)       = 1;
        nLayers(u_)             = MCS_and_scheduling{u_}.nLayers;
        M(u_)                   = MCS_and_scheduling{u_}.CQI_params.modulation_order;
        bittable{u_}            = LTE_params.bittable{M(u_)};
        symbol_alphabet{u_}     = LTE_params.SymbolAlphabet{M(u_)}.';
        CH_mapping_UE_spec{u_}  = UE_output(uG).UE_genie.CH_mapping(logical(kron(MCS_and_scheduling{u_}.UE_mapping,ones(LTE_params.Nsc,LTE_params.Ns))));
        CH_mapping_UE_spec{u_}  = repmat(reshape(CH_mapping_UE_spec{u_}, [], LTE_params.Nsub),[1 1 nRX]);
        N_subc(u_)              = sum(CH_mapping_UE_spec{u_}(:,1,1),1);
        H_est_user              = H_est_all(:,:,:,find(find(UEs_scheduled)==u_));
        H_all_noref{u_}         = reshape(H_est_user(CH_mapping_UE_spec{u_}), sum(CH_mapping_UE_spec{u_}(:,1,1),1), sum(CH_mapping_UE_spec{u_}(1,:,1),2), nRX) ;
    else
        M(u_) = 0;
    end
end


if (strcmp(filtering,'BlockFading')) % blockfading

    x = cell(N_user,1);
    for k = 1:LTE_params.Ntot
        % get number of scheduled users for current subcarrier k
        scheduled_UEs = [];
        for u_ = find(UEs_scheduled)
            uG = utils.localToGlobalUser(LTE_params.connection_table, bb, u_);
            scheduled_UEs(u_) = max(UE_output(uG).UE_genie.CH_mapping(k,:,1),[],2); %#ok
        end
        scheduled_UEs = find(scheduled_UEs);
        
        % build H matrix (for all scheduled UEs)
        H = zeros(nRX, length(scheduled_UEs));
        for u_ = 1:length(scheduled_UEs)
            H(:,u_) = squeeze(H_all_noref{scheduled_UEs(u_)}(k,1,:));
        end

        % calculate receiver
        switch receiver
            case 'ZF'
                inv_temp_k = inv(H'*H)*H';
            case 'MMSE'
                inv_temp_k = inv(H'*H + eye(size(H'*H))*sigma_n2)*H';
                Hg_tmp(k, :, :) = inv(H'*H + eye(size(H'*H))*sigma_n2)*H';
            otherwise
                error('unkown receiver type');
        end

        for n_ = 1:sum(max(CH_mapping_UE_spec{find(UEs_scheduled,1,'first')}(:,:,1),[],1))
            % apply receiver
            A  = inv_temp_k * squeeze(user_symbols(k,n_,:));
            
            inv_temp    = cell(N_user,1);
            
            us = 0;
            for u_  = scheduled_UEs
                us = us + 1;
                inv_temp{u_}(k) = mean(abs(inv_temp_k(us,:)));    % change this!
                x{u_}(k, n_) = A(us);
            end
        end
    end

    Hg = ones(N_user,1);

    
    for u_ = scheduled_UEs
        rx_layer_x = squeeze(x{u_});
        
        if LTE_params.DFT_spreading_off    % OFDMA
            rx_layer_x = rx_layer_x;
            noise_enhancement = sigma_n2*repmat(transpose(sum(abs(inv_temp{u_}).^2,2)),M(u_),1);
        else                               % SC-FDMA
            rx_layer_x = ifft(rx_layer_x,N_subc(u_))*sqrt(N_subc(u_));
            noise_enhancement = sigma_n2*mean(sum(abs(inv_temp{u_}).^2,2),1)*ones(M(u_),N_subc(u_)*sum(CH_mapping_UE_spec{u_}(1,:,1),2));
        end

        if LTE_params.DFT_spreading_off && ~strcmp(receiver,'ZF')  
                LLR_SD = zeros(size(noise_enhancement));
                s1 = size(rx_layer_x,1);
                s2 = size(rx_layer_x,2);
                s3 = size(rx_layer_x,3);
            for ii = 1:s1
                LLR_temp = LTE_demapper(reshape(rx_layer_x(ii,:,:),s2,s3).',symbol_alphabet{u_},bittable{u_},nLayers(u_),M(u_),Hg_full(:,ii),noise_enhancement(:,ii:s1:end),receiver);  
                LLR_SD(:,ii:s1:end) = LLR_temp;
            end
        else
                rx_layer_x = reshape(rx_layer_x,1,[]);
                LLR_SD = LTE_demapper(rx_layer_x,symbol_alphabet{u_},bittable{u_},nLayers(u_),M(u_),Hg(u_),noise_enhancement,receiver);
        end
        
        
        
        M1{u_} = M(u_);
        LLR_SD1{u_} = LLR_SD;
    end
    
else % FastFading
    x = cell(N_user,1);
    
    for k = 1:LTE_params.Ntot
        % get number of scheduled users for current subcarrier k
        scheduled_UEs = [];
        for u_ = find(UEs_scheduled)
            uG = utils.localToGlobalUser(LTE_params.connection_table, bb, u_);
            scheduled_UEs(u_) = max(UE_output(uG).UE_genie.CH_mapping(k,:,1),[],2); %#ok example 1 1 0 1 1
        end
        scheduled_UEs = find(scheduled_UEs); % example 1 2 4 5
        
        % build H matrix (for all scheduled UEs)
        H = zeros(nRX, length(scheduled_UEs));
        for n_ = 1:sum(max(CH_mapping_UE_spec{find(UEs_scheduled,1,'first')}(:,:,1),[],1))
            for u_ = 1:length(scheduled_UEs)
                H(:,u_) = squeeze(H_all_noref{scheduled_UEs(u_)}(k,n_,:));
            end

            % calculate receiver
            switch receiver
                case 'ZF'
                    inv_temp_k_n = inv(H'*H)*H';
                case 'MMSE'
                    inv_temp_k_n = inv(H'*H + eye(size(H'*H))*sigma_n2)*H';
                    Hg_tmp(k, :, :) = inv(H'*H + eye(size(H'*H))*sigma_n2);
                otherwise
                    error('unkown receiver type');
            end


            % apply receiver
            A  = inv_temp_k_n * squeeze(user_symbols(k,n_,:));
            
            inv_temp    = cell(N_user,1);
            for u_  = 1:length(scheduled_UEs)
                inv_temp{scheduled_UEs(u_)}(k) = mean(inv_temp_k_n(:,u_));   
                x{scheduled_UEs(u_)}(k, n_) = A(u_);
            end
        end
    end

    Hg = ones(N_user,1);

    
    for u_ = scheduled_UEs
        rx_layer_x = squeeze(x{u_});
        
        if LTE_params.DFT_spreading_off    % OFDMA
            noise_enhancement = sigma_n2*repmat(transpose(sum(abs(inv_temp{u_}).^2,2)),M(u_),1);
        else                               % SC-FDMA
            rx_layer_x = ifft(rx_layer_x,N_subc(u_))*sqrt(N_subc(u_));
            noise_enhancement = sigma_n2*mean(sum(abs(inv_temp{u_}).^2,2),1)*ones(M(u_),N_subc(u_)*sum(CH_mapping_UE_spec{u_}(1,:,1),2));
        end

        if LTE_params.DFT_spreading_off && ~strcmp(receiver,'ZF')  
                LLR_SD = zeros(size(noise_enhancement));
                s1 = size(rx_layer_x,1);
                s2 = size(rx_layer_x,2);
                s3 = size(rx_layer_x,3);
            for ii = 1:s1
                LLR_temp = LTE_demapper(reshape(rx_layer_x(ii,:,:),s2,s3).',symbol_alphabet{u_},bittable{u_},nLayers(u_),M(u_),Hg_full(:,ii),noise_enhancement(:,ii:s1:end),receiver);  
                LLR_SD(:,ii:s1:end) = LLR_temp;
            end
            else
                rx_layer_x = reshape(rx_layer_x,1,[]);
                LLR_SD = LTE_demapper(rx_layer_x,symbol_alphabet{u_},bittable{u_},nLayers(u_),M(u_),Hg(u_),noise_enhancement,receiver);
        end
        
        
        
        M1{u_} = M(u_);
        LLR_SD1{u_} = LLR_SD;
        
    end
end


end

