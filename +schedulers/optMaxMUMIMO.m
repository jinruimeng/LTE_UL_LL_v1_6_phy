classdef optMaxMUMIMO < network_elements.lteScheduler
    % Scheduler that maximizes the thorughput over many possible schedules by exhausitve search
    % Lukas Nagel, lnagel@nt.tuwien.ac.at
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at
    
    properties
        CQI_mapping_data
        alphabets
    end
    
    methods
        function obj = optMaxMUMIMO(params,UE,averager,mapping_data,varargin)
            % Fill in basic parameters (handled by the superclass constructor)
            obj = obj@network_elements.lteScheduler(params,UE,averager,mapping_data);
            obj.CQI_mapping_data = mapping_data;
            if ~isempty(varargin)
                obj.alphabets = varargin{1};
            else
                obj.alphabets = [];
            end

            % init UE_scheduling
            obj.UE_scheduling = obj.get_new_UE_params;
        end
        
        
        function UE_scheduling = scheduler_users(obj,subframe_corr,SyncUsedElements,BS,PBCHsyms,cqi_i,ChanMod_output,sigma_n2,channel_MSE)
            % assumptions: no scheduling, only feedback, no isi, no ici
            % additional assumption: the schedule stays constant for one
            % subframe ( otherwise the complexity would go through the sky! )
            % assumption N_ue > N_rx
            
            if ~obj.UE_config.CQI_fb
                error('MUMIMO exhaustiv search scheduler makes sense only with CQI_fb true');
            end

            N_rx    = obj.BS_config.nRX;             
            N_ue_total    = obj.nUEs;
            N_ue    = N_rx; % always use all receive antennas possible to schedule
            
            % generate all possible schedules
            N_combinations = nchoosek(N_ue_total, N_rx);
            possible_schedules = combnk(1:N_ue_total, N_ue);
            rates = zeros(1, N_combinations);
            
            for nn = 1:N_combinations
                current_schedule = possible_schedules(nn, :);

                for u_= 1:N_ue_total
                    if any(u_ == current_schedule)
                        obj.UE_scheduling(u_).UE_mapping    = true(obj.RB_grid_size,2);
                    else
                        obj.UE_scheduling(u_).UE_mapping    = false(obj.RB_grid_size,2);
                    end
                    obj.UE_scheduling(u_).assigned_RBs  = squeeze(sum(sum(obj.UE_scheduling(u_).UE_mapping,1),2));
                end

                UE_scheduling = obj.link_adaptation_MU_MIMO(obj.UE_scheduling,BS,cqi_i,ChanMod_output,sigma_n2,channel_MSE);

                for uu = 1:N_ue_total
                    if UE_scheduling(uu).assigned_RBs
                        rates(nn) = rates(nn) + UE_scheduling(uu).CQI_params.efficiency.*UE_scheduling(uu).assigned_RBs;
                    end
                end                
            end
            
            
            [~, index_opt] = max(rates); % find the best schedule
            optimal_schedule = possible_schedules(index_opt, :);
            % apply this schedule
            
            for u_= 1:N_ue_total
                if any(u_ == optimal_schedule)
                    obj.UE_scheduling(u_).UE_mapping    = true(obj.RB_grid_size,2);
                else
                    obj.UE_scheduling(u_).UE_mapping    = false(obj.RB_grid_size,2);
                end
                obj.UE_scheduling(u_).assigned_RBs  = squeeze(sum(sum(obj.UE_scheduling(u_).UE_mapping,1),2));
            end
            
            UE_scheduling = obj.link_adaptation_MU_MIMO(obj.UE_scheduling,BS,cqi_i,ChanMod_output,sigma_n2,channel_MSE);
            
            % calculate the number of bits
            obj.calculate_allocated_bits(UE_scheduling,subframe_corr,SyncUsedElements,PBCHsyms);
            
            for uu = 1:obj.nUEs
                %UE_scheduling(uu).total_no_refsym   = total_no_refsym;
                UE_scheduling(uu).SyncUsedElements  = SyncUsedElements;
                UE_scheduling(uu).CHusedElements    = PBCHsyms;
            end
                



        end


            
             %% end of: find the best SINR ---
            
        end
        
%         function UE_scheduling = link_adaptation_MU_MIMO(obj,UE_scheduling,BS,cqi_i,channel,sigma_n2,channel_MSE)
%             N_rx    = obj.BS_config.nRX;             
%             N_ue_total    = obj.nUEs;
%             N_ue    = N_rx; % always use all receive antennas possible to schedule
%             
%             
% 
%             % start of big loop!
% 
%             H_big = zeros(N_sub_total*N_ue, N_sub_total*N_rx);
%             F_big = zeros(N_sub_total*N_ue, N_sub_total*N_rx);
% 
%             for uu = current_schedule
%                 % assumption channel{uu} dimensions are: Nsub x Ntime x Nrx x Ntx
%                 H_tmp = squeeze(mean(channel{uu}, 2));
%                 H_basis{uu} = squeeze(H_tmp(:, :, 1));
%             end
% 
%             % tmp variable
%             E = eye(N_sub_total);
% 
%             for kk = 1:N_sub_total
% 
%                 scheduled_UEs = [];
%                 for u_ = current_schedule
%                     scheduled_subcarriers = kron(max(UE_scheduling(u_).UE_mapping,[],2),ones(12,1)); % TODO replace 12 here
%                     scheduled_UEs(u_) = scheduled_subcarriers(kk);
%                 end
%                 scheduled_UEs_idx = find(scheduled_UEs);
% 
%                 % create sub matrix
%                 H1 = zeros(N_rx, length(scheduled_UEs_idx));
% 
%                 for uu = 1:length(scheduled_UEs_idx)
%                     for rr = 1:N_rx
%                         H1(rr,uu) = H_basis{scheduled_UEs_idx(uu)}(kk, rr);
%                     end
%                 end
% 
%                 switch BS.receiver
%                     case 'ZF'
%                          F1 = inv(H1'*H1)*H1';
%                     case 'MMSE'
%                          F1 = inv(H1'*H1 + eye(length(scheduled_UEs_idx))*sigma_n2)*H1';
%                     otherwise
%                         error('receiver type not supported');
%                 end
% 
%                 %F1*H1
% 
%                 % DEBUG spratsch start
%                 K = F1*H1;
%                 for uu = 1:length(scheduled_UEs_idx)
%                     int_idx = scheduled_UEs(current_schedule);
%                     int_idx(uu) = 0;
%                     SINR_2(kk,uu) = abs(K(uu,uu)).^2 / ( sum(abs(K(uu,logical(int_idx))).^2) + sigma_n2*sum(abs(F1(uu,:)).^2) );
%                 end
%                 % DEBUG spratsch end
% 
%                 idx   = logical(kron(repmat(scheduled_UEs(current_schedule),[N_rx, 1]),diag(E(:,kk))));
%                 H_big(idx) = H1;
%                 F_big(idx') = F1;
% 
%             end
%             % DEBUG spratsch start
%             SINR_2 = harmmean(SINR_2);
%             % DEBUG spratsch end
% 
% 
%     %                 s_l = zeros(N_ue*N_sub,1);
%     %                 for uu = 1:N_ue
%     %                     s_l(((uu-1)*N_sub)+1:uu*N_sub,1) = kron(max(UE_scheduling(uu).UE_mapping,[],2),ones(12,1));
%     %                 end
%     %                 s_x = diag(s_l);
% 
%             % -- now calculate SINR
%             for uu = 1:N_ue
% 
%                 % first create S_l the matrix that cuts out the relevant
%                 % part
%                 S_l = zeros( N_sub_total,N_rx*N_sub_total);
%                 S_l(:,((uu-1)*N_sub_total)+1:uu*N_sub_total) = diag(kron(max(UE_scheduling(current_schedule(uu)).UE_mapping,[],2),ones(12,1)));%eye(N_sub);
% 
% 
% 
%                 N_sub_user = sum(max(UE_scheduling(current_schedule(uu)).UE_mapping,[],2))*12;
% 
%                 % calculate the parts of the formula
%                 P_signal = (1/N_sub_user)*abs(sum(S_l*diag(F_big*H_big)))^2;
%                 P_noise_enh = norm(S_l*F_big,'fro')^2 * sigma_n2; % this could be extended with channel_MSE % with or without S_l ? 
%                 P_interf= norm(S_l*F_big*H_big,'fro')^2-(1/N_sub_user)*abs(sum(S_l*diag(F_big*H_big)))^2;
% 
%                 SINR(uu) =  P_signal /(P_noise_enh + P_interf )   ;
% 
%             end
% 
%             % -- end SINR
% 
%             for uu = 1:N_ue
%                 SINR_AWGN = 10 * log10(SINR(uu));
%                 CQIs = (0:15);
%                 SINReff = obj.SINR_averager.average(10.^(SINR_AWGN(:)/10),CQIs,[obj.CQI_params(20).modulation_order,obj.CQI_params(CQIs(2:end)).modulation_order]);
% 
%                 CQI_temp(uu) = LTE_common_CQI_mapping_table(obj.CQI_mapping,SINReff,CQIs+1);
% 
%                 if(UE_scheduling(uu).assigned_RBs)  % if UE was scheduled
%                     UE_scheduling(uu).nCodewords    = 1;
%                     UE_scheduling(uu).nLayers       = 1;
%                     if obj.UE_config.CQI_fb
%                         if ( CQI_temp(uu) == 20 || CQI_temp(uu)==0 || isnan(CQI_temp(uu)) )
%                             UE_scheduling(uu).cqi       = 1;          % only take CQIs with throughput
%                         else
%                             UE_scheduling(uu).cqi       = CQI_temp(uu);
%                         end
%                         UE_scheduling(uu).CQI_params    = LTE_common_get_CQI_params(UE_scheduling(uu).cqi(1), obj.CQI_params);
%                     else
%                         UE_scheduling(uu).cqi           = cqi_i;
%                         UE_scheduling(uu).CQI_params    = LTE_common_get_CQI_params(UE_scheduling(uu).cqi(1),obj.CQI_params);
%                     end
%                 end
%             end
%         end    
     
end
