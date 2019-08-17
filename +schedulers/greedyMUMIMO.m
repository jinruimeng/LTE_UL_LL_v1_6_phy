classdef greedyMUMIMO < network_elements.lteScheduler
    % Scheduler that tries to improve the throughput by selecting a good user combination 
    % Lukas Nagel, lnagel@nt.tuwien.ac.at
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at
    
    properties
        CQI_mapping_data
        alphabets
    end
    
    methods
        function obj = greedyMUMIMO(params,UE,averager,mapping_data,varargin)
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
            % assumption N_ue >= N_rx
            

            N_rx    = obj.BS_config.nRX;             
            N_ue_total    = obj.nUEs;
            alpha = 0.4; % threshold
            
            % scheduling algorithm

            
            
             % assumption channel{uu} dimensions are: Nsub x Ntime x Nrx x Ntx
            H = zeros(N_rx, N_ue_total);
            
            for uu = 1:N_ue_total
                H(:, uu) = squeeze(mean(mean(ChanMod_output{uu},2),1)); % makes sense to do the mean over subcarrier this way?
            end
           
            % initialisation
            T = H;
            S = [];
            i = 1;
            G = [];
            
            
            while true
                
                B = eye(N_rx);
                for n = 1:i-1
                    %s = S(n);
                    % G instead of H
                    %B = B - H(:,s)*H(:,s)'/norm(H(:,s))^2;
                    B = B - G(:,n)*G(:,n)'/norm(G(:,n))^2;
                    
                end
                
                if i == 1
                    Gf = T;
                else
                    %G = B*T;
                    Gf = (eye(size(B))+B)*T;
                end
                
                Gn = zeros(1,size(Gf, 2));
                
                for n = 1:length(Gn)
                    Gn(n) = norm(Gf(:,n));
                end
                
                [~ , idx] = max(Gn);
                gnew = Gf(:,idx);
                G = [G gnew];
                S = [S, find(ismember(H.',T(:,idx).','rows'))];

                
                Tnew = [];
                
                for k = 1:size(T,2)
                    if k ~= idx
                        if T(:,k)'*gnew / ( norm(T(:,k)) * norm(gnew) ) < alpha
                            Tnew = [Tnew, T(:,k)];
                        end
                    end
                end
                
                i = i + 1;
                
                
                if length(S) == N_rx
                    %disp('done!')
                    break
                end
                
                if size(Tnew, 2) == 0%< N_rx - length(S)
                    %first pick the new ones
                    % for now pick don't pick something randomly
                    
                    %others = setdiff(1:N_ue_total, S);
                    %S = [S, others(1:(N_rx - length(S)))];

                    break;
                end
                
                T = Tnew;

                
            end
            
            
            
            % end scheduling algorithm
           
          
            % apply this schedule
            for u_= 1:N_ue_total
                if any(u_ == S)
                    obj.UE_scheduling(u_).UE_mapping    = true(obj.RB_grid_size,2);
                else
                    obj.UE_scheduling(u_).UE_mapping    = false(obj.RB_grid_size,2);
                end
                obj.UE_scheduling(u_).assigned_RBs  = squeeze(sum(sum(obj.UE_scheduling(u_).UE_mapping,1),2));
            end
            
            % apply link adaptation
            UE_scheduling = obj.link_adaptation_MU_MIMO(obj.UE_scheduling,BS,cqi_i,ChanMod_output,sigma_n2,channel_MSE);
            
            % calculate the number of bits
            obj.calculate_allocated_bits(UE_scheduling,subframe_corr,SyncUsedElements,PBCHsyms);
            
            for uu = 1:obj.nUEs
                %UE_scheduling(uu).total_no_refsym   = total_no_refsym;
                UE_scheduling(uu).SyncUsedElements  = SyncUsedElements;
                UE_scheduling(uu).CHusedElements    = PBCHsyms;
            end
        end
        

    end %end methods   

     
end
