classdef lteScheduler < handle
    % Implements methods common to all the schedulers.
    % Stefan Pratschner, spratsch@nt.tuwien.ac.at
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at
    
    properties
        UEs                 % List of UEs to schedule
        nUEs                % Total number of UEs -> length(UEs)
        UE_scheduling       % output struct
        UE_config           % UE configuration
        BS_config           % eNodeB configuration
        ChanMod_config      % channel model configuration
        Ns_RB               % Number of symbols in one RB (needed to calculate the TB sizes)
        Ns                  % number of symbols per slot
        Ntot                % total number of uplink subcarriers
        Nfft                % FFT size
        Tg                  % guard time
        SamplingTime        % 
        ICI_power           % 
        MI_data             %
        
        max_HARQ_retx       % Max num of HARQ retransmissions, NOT including the original tx. 0, 1, 2 or 3
        RB_grid_size        % Size of the Resource Block grid
        
        CQI_params          % CQI parameters for all possible MCSs
        CQI_mapping         % CQI mapping parameters
        
        
        CQI_mapping_params  % parameters needed to perform the CQI mapping
        
        UE_specific         % direct reference to the HARQ processes from the eNodeB
        first               % indicator variable
        fairness            % desired fairness (for var. fair scheduler)
        av_const            % averaging constant for the exponential throughput averaging
        uplink              % for uplink is changed calculation of data bits, Prokopec
        fixed_assignment    % fixed RB assignment for fixed scheduler
        DFT_spreading_off   % to get OFDM performance
        SINR_averager       % SINR averyger type
    end
    
    methods
        % Superclass constructor which is called by all subclasses
        function obj = lteScheduler(params,UEs_to_be_scheduled,averager,mapping_data)

            obj.RB_grid_size        = params.Nrb;
            obj.Ns_RB               = params.Nsc*params.Ns;
            obj.Ns                  = params.Ns;
            obj.Ntot                = params.Ntot;
            obj.Nfft                = params.Nfft;
            obj.Tg                  = params.Tg;
            obj.SamplingTime        = params.SamplingTime;
            obj.ICI_power           = params.ICI_power;  
            obj.MI_data             = params.MI_data;
            
            obj.UEs                 = UEs_to_be_scheduled;
            obj.nUEs                = length(UEs_to_be_scheduled);
            obj.UE_config           = params.UE_config;
            obj.BS_config           = params.BS_config;
            obj.ChanMod_config      = params.ChanMod_config;
            
            obj.max_HARQ_retx       = params.scheduler.max_HARQ_retx;
            obj.CQI_params          = params.CQI_params;
            obj.CQI_mapping         = params.CQI_mapping;
            
            obj.CQI_mapping_params  = params.scheduler.CQI_mapping_params;
            obj.UE_specific         = params.scheduler.UE_specific;
%             obj.av_const            = params.scheduler.av_window;
%             obj.fairness            = params.scheduler.fairness;
%             obj.first               = true;
%             obj.uplink              = params.Simulation_type;
            obj.fixed_assignment    = params.scheduler.fixed_scheduler_assignment;
            
            obj.DFT_spreading_off   = params.DFT_spreading_off;
            obj.SINR_averager       = averager;
        end
        
        % Returns nUEs empty scheduler parameter objects
        function empty_UE_scheduling_params = get_new_UE_params(obj)
            empty_UE_scheduling_params = network_elements.UE_scheduling_params;
            for u_ = 2:obj.nUEs
                empty_UE_scheduling_params(u_) = network_elements.UE_scheduling_params;
            end
        end
        
        
        function calculate_allocated_bits(obj,UE_scheduling,subframe_corr,SyncUsedElements,ChUsedElements)
            
            % You can access the properties of the eNodeB like this
            UE_nTX = obj.UE_config.nTX;
            
            for u_=1:length(UE_scheduling)
                
                
                CHsyms = sum(sum(ChUsedElements(UE_scheduling(u_).UE_mapping)));
                % Calculation of the TB size
                if ~obj.uplink % PRokopec
                    if(subframe_corr == 1 || subframe_corr == 6)
                        %lenghts of primary and secondary synchronization channels (symbols)
                        sync_symbols = sum(sum(SyncUsedElements(UE_scheduling(u_).UE_mapping)));
                    else
                        sync_symbols = 0;
                    end
                else
                    sync_symbols = sum(sum(SyncUsedElements(UE_scheduling(u_).UE_mapping)));
                end % end Prokopec
                
                if(UE_scheduling(u_).assigned_RBs)
                    assigned_RBs = UE_scheduling(u_).assigned_RBs;
%                    switch UE_scheduling(u_).nCodewords
%                        case 1
%                            assigned_RBs = UE_scheduling(u_).assigned_RBs*2;
%                        case 2
%                            assigned_RBs = UE_scheduling(u_).assigned_RBs;
%                    end
                    UE_scheduling(u_).N_coded_bits = [UE_scheduling(u_).CQI_params.modulation_order] * (assigned_RBs * obj.Ns_RB - sync_symbols - CHsyms);
                    UE_scheduling(u_).N_data_bits  = 8*round(1/8 * UE_scheduling(u_).N_coded_bits .* [UE_scheduling(u_).CQI_params.coding_rate_x_1024] / 1024)-24; % calculate G based on TB_size and the target rate
                    for i_ = 1:length(UE_scheduling(u_).N_data_bits)
                        while 1
                            if(UE_scheduling(u_).N_data_bits(i_) < 16)
                                UE_scheduling(u_).cqi(i_)           = UE_scheduling(u_).cqi(i_)+1;
                                UE_scheduling(u_).CQI_params(i_)    = LTE_common_get_CQI_params(UE_scheduling(u_).cqi(i_),obj.CQI_params);
                                UE_scheduling(u_).N_coded_bits(i_)  = UE_scheduling(u_).CQI_params(i_).modulation_order * (assigned_RBs * obj.Ns_RB - sync_symbols - CHsyms);
                                UE_scheduling(u_).N_data_bits(i_)   = 8*round(1/8 * UE_scheduling(u_).N_coded_bits(i_) .* UE_scheduling(u_).CQI_params(i_).coding_rate_x_1024 / 1024)-24; % calculate G based on TB_size and the target rate
                            else
                                break
                            end
                        end
                    end
                    if sum(UE_scheduling(u_).N_data_bits) == 0
                        UE_scheduling(u_).assigned_RBs = 0;
                        UE_scheduling(u_).UE_mapping   = false(obj.RB_grid_size,2);
                    end
                else
                    UE_scheduling(u_).N_coded_bits = 0;
                    UE_scheduling(u_).N_data_bits = 0;
                end
                
                % Adjust number of bits for the SM case with 4 antennas
                if (obj.UE_config.mode ~= 1 && UE_nTX == 4)       % spratsch: why in 4 antenna case ?
                    switch UE_scheduling(u_).nLayers
                        case 2
                            if (UE_scheduling(u_).nCodewords == 1)
                                UE_scheduling(u_).N_coded_bits = UE_scheduling(u_).N_coded_bits*2;
                                UE_scheduling(u_).N_data_bits = UE_scheduling(u_).N_data_bits*2;
                            end
                        case 3
                            UE_scheduling(u_).N_coded_bits(2) = UE_scheduling(u_).N_coded_bits(2)*2;
                            UE_scheduling(u_).N_data_bits(2) = UE_scheduling(u_).N_data_bits(2)*2;
                        case 4
                            UE_scheduling(u_).N_coded_bits = UE_scheduling(u_).N_coded_bits*2;
                            UE_scheduling(u_).N_data_bits = UE_scheduling(u_).N_data_bits*2;
                    end
                end
                
                % Adjust information from the HARQ processes
                
                % %                if LTE_params.two_dim == 'off'
                %                    for cw_=1:UE_scheduling(u_).nCodewords
                %                        % Update HARQ process information
                %                        obj.UE_specific(u_).current_HARQ_process(cw_).assigned_RBs = UE_scheduling(u_).assigned_RBs;
                %                        obj.UE_specific(u_).current_HARQ_process(cw_).cqi          = UE_scheduling(u_).cqi;
                %                        % rv_idx information of the HARQ process is updated outside of the scheduler
                %                        % Update scheduler information only
                %                        UE_scheduling(u_).HARQ_process_id(cw_) = obj.UE_specific(u_).current_HARQ_process(cw_).id;
                %                        UE_scheduling(u_).rv_idx(cw_) = obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx;
                %                    end
                % %                end
                
                for cw_=UE_scheduling(u_).nCodewords+1:size(obj.UE_specific(u_).HARQ_processes,1)
                    obj.UE_specific(u_).current_HARQ_process(cw_).assigned_RBs = 0;
                    obj.UE_specific(u_).current_HARQ_process(cw_).cqi          = 0;
                    % rv_idx information of the HARQ process is updated outside of the scheduler
                    % Update scheduler information only
                    UE_scheduling(u_).rv_idx(cw_) = 0;
                end
                
            end
        end
        
      
        % link adaption
        function UE_scheduling = link_adaptation_SIxO(obj,UE_scheduling,BS,cqi_i,channel,sigma_n2,channel_MSE)
            N_UE    = obj.nUEs;
            nRX     = obj.BS_config.nRX;
            
            if obj.UE_config.CQI_fb
                % with link adaption
                for uu = 1:N_UE
                    if(UE_scheduling(uu).assigned_RBs)  % if UE was scheduled
                        for slot_i = 1
                            SNR = [];   % to avoid dimension mismatch
                            
                            % find UE mapping
                            ue_allocation = kron(UE_scheduling(uu).UE_mapping(:,1),ones(12,1));
                            
                            % get scheduled channel of UE
                            N           = sum(ue_allocation);
                            H_temp      = channel{uu};
%                             H_mean      = squeeze(mean(H_temp(:,(slot_i-1)*obj.Ns+1:slot_i*obj.Ns,:,:),2));
                            H_mean      = squeeze(mean(H_temp,2));
                            H_mean_user = reshape(H_mean(logical(repmat(ue_allocation,[1 nRX]))),N,[]);   
                            
                            % calculate ISI and ICI power
                            if obj.UE_config.ignore_ISI_ICI
                                % ISI
                                ISI_power       = zeros(N,1);
                                
                                % ICI
                                ICI_power_tmp   = zeros(1,N);
                            else
                                % ISI
                                H_cal           = LTE_UL_estimate_ISI_power(obj.ChanMod_config,obj.Tg,obj.SamplingTime,obj.Nfft,obj.Ntot,channel,uu);
                                ISI_power       = 2.*diag( squeeze(H_cal(:,:,1,1))'*squeeze(H_cal(:,:,1,1)) ) / obj.Nfft;
                                ISI_power       = ISI_power((obj.Nfft-obj.Ntot)/2+1:(obj.Nfft-obj.Ntot)/2+obj.Ntot);
                                ISI_power       = ISI_power(logical(ue_allocation));
                                
                                % ICI
                                ICI_power_tmp   = obj.ICI_power(logical(ue_allocation));
                            end
                            
                            % channel MSE
                            if obj.UE_config.ignore_channel_estimation
                                channel_MSE(uu) = 0;
                            end
                            
                            % calculate receiver                            
                            switch BS.receiver
                                case 'ZF'
%                                     F = pinv(diag(H_mean_user));
                                    F = conj(H_mean_user)./repmat(sum(H_mean_user.*conj(H_mean_user),2),1,nRX);
                                case 'MMSE'
%                                     F = diag(conj(H_mean_user)./sum(H_mean_user.*conj(H_mean_user)+sigma_n2(uu),2));
                                    F = (conj(H_mean_user).*H_mean_user +sigma_n2).^(-1).*conj(H_mean_user);
                                case 'IAMMSE'
                                    % for now just apply normal MMSE for
                                    % feedback. TODO: add IAMMSE here...
                                    F = (conj(H_mean_user).*H_mean_user +sigma_n2).^(-1).*conj(H_mean_user);
                                otherwise
                                    error('receiver type not supported');
                            end

                            % calculate effective SINR
                            FH_tmp = sum(F.*H_mean_user,2); % receiver times channel temp (also for RX diversity)
                            if obj.DFT_spreading_off    % OFDMA SINR
                                SNR(slot_i,:) = (FH_tmp.^2)./ ((ISI_power+ICI_power_tmp.'+sigma_n2).*sum(abs(F).^2,2));
                            else                        % SC-FDMA SNR
%                                 SNR(slot_i,:) = (1/N) * abs(sum(FH_tmp))^2 ./ ( norm(FH_tmp)^2 - (1/N)*abs(sum(FH_tmp))^2 + (F(:)'*F(:))*ISI_power + N*channel_MSE(uu) +(sigma_n2(uu)+ICI_power_tmp.')*norm(F(:))^2);
                                SNR(slot_i,:) = (1/N) * abs(sum(FH_tmp))^2 ./ ( norm(FH_tmp)^2 - (1/N)*abs(sum(FH_tmp))^2 + (F(:)'*F(:))*ISI_power + (channel_MSE(uu) + sigma_n2+ICI_power_tmp.')*norm(F(:))^2);
                            end
                        end

                        % appply SNR averaging for the case of ISI/ICI 
                        SNR_mean = obj.SINR_averager.average(SNR(:),1);

                        % get CQI from effective SINR
                        CQI_temp = LTE_common_CQI_mapping_table(obj.CQI_mapping_params,SNR_mean,(0:15)+1);
                        if ( CQI_temp == 20 || CQI_temp==0 || isnan(CQI_temp) )
                            UE_scheduling(uu).cqi       = 1;          % only take CQIs with throughput
                        else
                            UE_scheduling(uu).cqi       = CQI_temp;
                        end
                        UE_scheduling(uu).nCodewords    = 1;
                        UE_scheduling(uu).nLayers       = 1;
                        UE_scheduling(uu).CQI_params    = LTE_common_get_CQI_params(UE_scheduling(uu).cqi(1),obj.CQI_params);
                    end
                end     
            else
                % without link adaption
                for uu = 1:N_UE
                    if(UE_scheduling(uu).assigned_RBs)  % if UE was scheduled
                        UE_scheduling(uu).nCodewords    = 1;
                        UE_scheduling(uu).nLayers       = 1;
                        UE_scheduling(uu).cqi           = cqi_i;
                        UE_scheduling(uu).CQI_params    = LTE_common_get_CQI_params(UE_scheduling(uu).cqi(1),obj.CQI_params);
                    end
                end
            end
        end
        
        function UE_scheduling = link_adaptation_MU_MIMO(obj,UE_scheduling,BS,cqi_i,channel,sigma_n2,channel_MSE)
            % assumptions: no scheduling, only feedback, no isi, no ici, only nUE=nRX
            
            N_ue    = obj.nUEs;
            N_sub_total   = obj.Ntot; 
            N_rx    = obj.BS_config.nRX;
            
            for uu = 1:N_ue
                scheduled_UEs_total(uu) =  max(max(UE_scheduling(uu).UE_mapping));
            end
            scheduled_UEs_idx_total = find(scheduled_UEs_total);
            
            N_scheduled = length(scheduled_UEs_idx_total);
            
            if obj.UE_config.CQI_fb
                H_big = zeros(N_sub_total*N_rx, N_sub_total*N_scheduled);
                F_big = zeros(N_sub_total*N_scheduled, N_sub_total*N_rx);

                for uu = 1:N_ue
                    % assumption channel{uu} dimensions are: Nsub x Ntime x Nrx x Ntx
                    H_tmp = squeeze(mean(channel{uu}, 2));
                    H_basis{uu} = squeeze(H_tmp(:, :, 1));
                end

                % tmp variable
                E = eye(N_sub_total);

                for kk = 1:N_sub_total
                    
                    scheduled_UEs = [];
                    for u_ = 1:N_ue
                        scheduled_subcarriers = kron(max(UE_scheduling(u_).UE_mapping,[],2),ones(12,1)); % TODO replace 12 here
                        scheduled_UEs(u_) = scheduled_subcarriers(kk);
                    end
                    scheduled_UEs_idx = find(scheduled_UEs);
                    
                    % create sub matrix
                    H1 = zeros(N_rx, length(scheduled_UEs_idx));
                    
                    for uu = 1:length(scheduled_UEs_idx)
                        for rr = 1:N_rx
                            H1(rr,uu) = H_basis{scheduled_UEs_idx(uu)}(kk, rr);
                        end
                    end

                    switch BS.receiver
                        case 'ZF'
                             F1 = inv(H1'*H1)*H1';
                        case 'MMSE'
                             F1 = inv(H1'*H1 + eye(length(scheduled_UEs_idx))*sigma_n2)*H1';
                        otherwise
                            error('receiver type not supported');
                    end
                    
                    idx   = logical(kron(repmat(scheduled_UEs(scheduled_UEs~=0),[N_rx, 1]),diag(E(:,kk))));
                    H_big(idx) = H1;
                    F_big(idx') = F1;
                end

                
                
%                 s_l = zeros(N_ue*N_sub,1);
%                 for uu = 1:N_ue
%                     s_l(((uu-1)*N_sub)+1:uu*N_sub,1) = kron(max(UE_scheduling(uu).UE_mapping,[],2),ones(12,1));
%                 end
%                 s_x = diag(s_l);

                % -- now calculate SINR
                
                D = fft(eye(N_sub_total))/sqrt(N_sub_total);
                K=kron(eye(length(scheduled_UEs_idx_total)),D')*F_big*H_big*kron(eye(length(scheduled_UEs_idx_total)),D);
                for uu = 1:length(scheduled_UEs_idx_total)

                    % first create S_l the matrix that cuts out the relevant
                    % part
                    S_l = zeros( N_sub_total,N_scheduled*N_sub_total);
                    S_l(:,((uu-1)*N_sub_total)+1:uu*N_sub_total) = diag(kron(max(UE_scheduling(scheduled_UEs_idx_total(uu)).UE_mapping,[],2),ones(12,1)));%eye(N_sub);
                    
                    N_sub_user = sum(max(UE_scheduling(scheduled_UEs_idx_total(uu)).UE_mapping,[],2))*12;
                    
                    % calculate the parts of the formula
%                     P_signal = (1/N_sub_user)*abs(sum(S_l*diag(F_big*H_big)))^2;
%                     P_noise_enh = norm(S_l*F_big,'fro')^2 * sigma_n2; % this could be extended with channel_MSE % with or without S_l ? 
%                     P_interf= norm(S_l*F_big*H_big,'fro')^2-(1/N_sub_user)*abs(sum(S_l*diag(F_big*H_big)))^2;
                    
                    P_signal      = norm(S_l*diag(K),'fro')^2;
                    P_noise       = sigma_n2 * norm(D'*S_l*F_big,'fro')^2;
                    P_interf      = norm(S_l*K,'fro')^2 - norm(S_l*diag(K),'fro')^2;
                    if obj.UE_config.ignore_ISI_ICI
                        P_ici = 0;
                    else
                        P_ici         = sum(obj.ICI_power)*norm(D'*S_l*F_big,'fro')^2;
                    end
                    
                    SINR(scheduled_UEs_idx_total(uu)) = P_signal /(P_noise + P_interf + P_ici );
                end

                % -- end SINR

            
             
                for uu = 1:N_scheduled
                    SINR_AWGN = 10 * log10(SINR(scheduled_UEs_idx_total(uu)));
                    CQIs = (0:15);
                    SINReff = obj.SINR_averager.average(10.^(SINR_AWGN(:)/10),CQIs,[obj.CQI_params(20).modulation_order,obj.CQI_params(CQIs(2:end)).modulation_order]);
                
                    CQI_temp(scheduled_UEs_idx_total(uu)) = LTE_common_CQI_mapping_table(obj.CQI_mapping,SINReff,CQIs+1);
                end
                
            end



            
            for uu = 1:N_scheduled
                global_UE = scheduled_UEs_idx_total(uu);
                if(UE_scheduling(global_UE).assigned_RBs)  % if UE was scheduled
                    UE_scheduling(global_UE).nCodewords    = 1;
                    UE_scheduling(global_UE).nLayers       = 1;
                    if obj.UE_config.CQI_fb
                        if ( CQI_temp(global_UE) == 20 || CQI_temp(global_UE)==0 || isnan(CQI_temp(global_UE)) )
                            UE_scheduling(global_UE).cqi       = 1;          % only take CQIs with throughput
                        else
                            UE_scheduling(global_UE).cqi       = CQI_temp(global_UE);
                        end
                        UE_scheduling(global_UE).CQI_params    = LTE_common_get_CQI_params(UE_scheduling(global_UE).cqi(1), obj.CQI_params);
                    else
                        UE_scheduling(global_UE).cqi           = cqi_i;
                        UE_scheduling(global_UE).CQI_params    = LTE_common_get_CQI_params(UE_scheduling(global_UE).cqi(1),obj.CQI_params);
                    end
                end
            end

        end       
    
        function UE_scheduling = link_adaptation_CLSM(obj,UE_scheduling,channel,cqi_i,sigma_n2,channel_MSE)
            N_UE    = obj.nUEs;
            
           
            % call CLSM FB function
            for uu = 1:N_UE % we have to handle the inter user interference like inter layer interference in normal CLSM
                if(UE_scheduling(uu).assigned_RBs)  % if UE was scheduled
                    % find UE mapping
                    ue_allocation = kron(UE_scheduling(uu).UE_mapping(:,1),ones(12,1));
                    
                    % get scheduled UE channel
                    H_perf = channel{uu};
                    H_user = H_perf(logical( repmat(ue_allocation,[1 size(H_perf,2), size(H_perf,3), size(H_perf,4)]) ));
                    H_user = reshape(H_user, [], size(H_perf,2), size(H_perf,3), size(H_perf,4));
                    
                    % get scheduled ICI
                    ICI_power_user = obj.ICI_power;
                    ICI_power_user = ICI_power_user(logical(ue_allocation.'));
                    
                    % get ISI
                    if obj.UE_config.ignore_ISI_ICI            % ISI/ICI will be ignored if set so            
%                         H_cal = zeros(-1, LTE_params.Nfft, ChanMod.nRX, ChanMod.nTX);
                        H_cal = [];
                    else
                        H_cal = LTE_UL_estimate_ISI_power(obj.ChanMod_config,obj.Tg,obj.SamplingTime,obj.Nfft,obj.Ntot,channel,uu);
                    end
                    
                    % channel MSE
                    if obj.UE_config.ignore_channel_estimation
                        channel_MSE(uu) = 0;
                    end
                    
                    % calculate link adaption
                    UE_scheduling(uu).cqi = [];     % remove previous value
%                     if (obj.UE_config.RI_fb || obj.UE_config.PMI_fb || obj.UE_config.CQI_fb)
                        [RI,PMI,CQI] = LTE_UL_link_adaption_CLSM(obj.UE_config,obj.BS_config,obj.CQI_params,obj.CQI_mapping,obj.SINR_averager,obj.DFT_spreading_off,obj.MI_data,obj.Ntot,obj.Ns,obj.Nfft,obj.RB_grid_size,ICI_power_user,sigma_n2,channel_MSE(uu),H_user,H_cal,cqi_i,ue_allocation);
                        UE_scheduling(uu).nLayers       = RI;
                        UE_scheduling(uu).PMI           = PMI;
                        UE_scheduling(uu).nCodewords    = min(2,RI);
%                     else
%                         UE_scheduling(uu).nLayers       = obj.UE_config.RI;
%                         UE_scheduling(uu).PMI           = obj.UE_config.PMI;
%                         UE_scheduling(uu).nCodewords    = min(2,RI);
%                         CQI                             = cqi_i*ones(min(2,obj.UE_config.RI),obj.UE_config.RI);
%                     end
                    
                    for ii = 1:min(2,RI)  % only for used codewords
                        if( CQI(ii) == 20 || CQI(ii)==0 || isnan(CQI(ii)) )
                                UE_scheduling(uu).cqi(1,ii) = 1;          % only take CQIs with throughput
                            else
                                UE_scheduling(uu).cqi(1,ii) = CQI(ii);
                        end
                    end
                    
                    UE_scheduling(uu).CQI_params = LTE_common_get_CQI_params(UE_scheduling(uu).cqi,obj.CQI_params);
                end
            end
        end
        
        function [schedules] = generate_schedules(obj)
            % generate all possible UL schedules

            % initialization
            nRB = obj.RB_grid_size;
            nUE = obj.nUEs;
            schedules = [];
            temp_schedule = zeros(nRB,1);
            running = true;

            % generate schedules
            while(running)
                temp_schedule(1) = temp_schedule(1) + 1;                        % increment
                for RB_i=1:nRB-1
                    if(temp_schedule(RB_i)>nUE)                                 % check for carry
                        temp_schedule(RB_i) = 0;                                % mod nRB
                        temp_schedule(RB_i+1) = temp_schedule(RB_i+1) + 1;      % increment next bit
                    end
                end

                if(temp_schedule(nRB)>nUE)                                      % if last bit is full
                    running = false;                                            % break condition
                else
                    if(sum(temp_schedule>0)<=nRB)                               % if not too many RBs are already allocated
                        valid_schedule = 1;                                     % check if it is a valid schedule
                        for UE_i = 1:nUE
                            if( sum(abs(diff([0;temp_schedule;0]==UE_i)))>2)    % check that user is assigned consecutive RBs only
                                valid_schedule = 0;
                            end
                        end            

                        if valid_schedule                                       % when the generated schedule is a valid one
                            schedules = [schedules, temp_schedule];             % add scheulde
                        end
                    end
                end
            end
        end
    
    end

    
%     methods (Abstract)
%         % UE scheduling (to be implemented for each subclass
%         UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_specific_data)
%     end
end
