classdef approxMaxThroughputScheduler < network_elements.lteScheduler
    % Scheduler that assigns a fixed and predefined amount of RBs to each user
    % Stefan Pratschner, spratsch@nt.tuwien.ac.at
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at
    
    properties
        CQI_mapping_data
        alphabets
    end
    
    methods
        function obj = approxMaxThroughputScheduler(params,UE,averager,mapping_data,varargin)
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
        
        
        function UE_scheduling = scheduler_users(obj,subframe_corr,SyncUsedElements,BS,PBCHsyms,cqi_i,channel,sigma_n2,channel_MSE)
            % heuristic resource allocation
            % similiar to recursive maximum expansion (RME)
            
            % initialization
            mean_SNR_RB         = NaN(obj.RB_grid_size,obj.nUEs);
            N_sub               = 12;                                   % make this dependent on LTE_params!
            nRX                 = size(channel{1},3);
            
            % calculate SNR for each RB and UE
            for uu = 1:obj.nUEs
                for RB_i = 1:obj.RB_grid_size
                    % get channel for each RB
                    H_mean = squeeze(mean(channel{uu},2));
                    H_mean_RB = H_mean((RB_i-1)*N_sub+1:RB_i*N_sub,:);            
                    N = N_sub;                                          % make this dependent on LTE_params!

                    % calculate effective SINR
                    switch BS.receiver
                        case 'ZF'
%                             F = pinv(diag(H_mean_RB));
                            F = conj(H_mean_RB)./repmat(sum(H_mean_RB.*conj(H_mean_RB),2),1,nRX);
                        case 'MMSE'
%                             F = diag(conj(H_mean_RB)./sum(H_mean_RB.*conj(H_mean_RB)+sigma_n2,2));
                            F = (conj(H_mean_RB).*H_mean_RB +sigma_n2).^(-1).*conj(H_mean_RB);
                        otherwise
                            error('receiver type not supported by approx. max throughput scheduler');
                    end

                    ISI_power = zeros(length(F),1);
                    ICI_power = zeros(length(F),1);
                    
                    FH_tmp = sum(F.*H_mean_RB,2);
                    if obj.DFT_spreading_off    % OFDMA SNR
%                         SNR = ((F*H_mean_RB).^2)./ ((ISI_power+ICI_power+sigma_n2).*abs(diag(F)).^2);
                        SNR = ((FH_tmp).^2)./ ((ISI_power+ICI_power+sigma_n2).*sum(abs(F).^2,2));
                    else                        % SC-FDMA SNR
%                         SNR = (1/N) * abs(sum(F*H_mean_RB))^2 ./ ( norm(F*H_mean_RB)^2 - (1/N)*abs(sum(F*H_mean_RB))^2 + (F'*F)*ISI_power + (sigma_n2+ICI_power)*norm(diag(F))^2 );
                        SNR = (1/N) * abs(sum(FH_tmp))^2 ./ ( norm(FH_tmp)^2 - (1/N)*abs(sum(FH_tmp))^2 + (F(:)'*F(:))*ISI_power + (sigma_n2+ICI_power)*norm(F(:))^2 );
                    end

                    % calculate SNR for each RB and UE
                    mean_SNR_RB(RB_i,uu) = 10.^(obj.SINR_averager.average(SNR,1)/10);
                end
            end
            
%             % heuristic resource allocation
%             sched_size          = size(mean_SNR_RB);
%             mask                = zeros(obj.RB_grid_size,1);
%             mean_SNR_RB_temp    = mean_SNR_RB;                          % this will be overwritten
%             
%             invalid_schedule    = 1;
%             while invalid_schedule
%                 [UE_SNR_max_values,UE_SNR_max] = max(mean_SNR_RB_temp,[],1);    % best SNR of each UE
%                 [~,sort_ind] = sort(UE_SNR_max_values);                         % order UEs by increasing max SNR
%                 invalid_schedule = 0; 
%                 for UE_i = sort_ind
%                 [~,schedule_RME]   = max(mean_SNR_RB_temp,[],2);            % give RB to UE with best SNR as intial guess
% 
%                     if( sum(abs(diff([0;schedule_RME;0]==UE_i)))>2 )        % if there is more than one cluster of RBs assigned to a UE
%                         invalid_schedule = 1;                               % the schedule is not a valid one yet
% 
%                         % calculate the assigned clusters
%                         clusters    = diff([0;schedule_RME;0]==UE_i);
%                         start       = find(clusters==1);
%                         stop        = find(clusters==-1)-1;
% 
%                         % find largest SNR within assigned clusters
%                         [~,max_SNR_ind]=max(mean_SNR_RB_temp(:,UE_i).*(schedule_RME==UE_i));
% 
%                         for ii=1:length(start)                                                      % look for the cluster where the UE has the best SNR
%                             if( max_SNR_ind>=start(ii,1) && max_SNR_ind<=stop(ii,1) )
%                                 mask(start(ii):stop(ii),1) = ones(stop(ii)-start(ii)+1,1);          % mask the best cluster of RBs for the UE
%                                 mean_SNR_RB_temp(:,UE_i) = mask .* mean_SNR_RB_temp(:,UE_i);
%                             end
%                         end        
%                     end
%                     if invalid_schedule
%                         break
%                     end
%                 end
%             end

            %% RME

            % take maxima of each UE
            mean_SNR_RB_temp                = mean_SNR_RB;
            [UE_SNR_max_values,UE_SNR_max]  = max(mean_SNR_RB_temp,[],1);   % best SNR of each UE
            [~,sort_ind]                    = sort(UE_SNR_max_values);      % order UEs by increasing max SNR
            [~,schedule_RME]                = max(mean_SNR_RB_temp,[],2);   % give RB to UE with best SNR as intial guess
            run = 1;
            while run
                for UE_i = sort_ind
                    if( sum(abs(diff([0;schedule_RME;0]==UE_i)))>2 )            % if there is more than one cluster of RBs assigned to a UE
                        % calculate the assigned clusters
                        clusters    = diff([0;schedule_RME;0]==UE_i);
                        start       = find(clusters==1);
                        stop        = find(clusters==-1)-1;

                        % find largest SNR within assigned clusters
                        [~,max_SNR_ind]=max(mean_SNR_RB_temp(:,UE_i).*(schedule_RME==UE_i));

                        mask = zeros(obj.RB_grid_size,1);
                        for ii=1:length(start)                                                      % look for the cluster where the UE has the best SNR
                            if( max_SNR_ind>=start(ii,1) && max_SNR_ind<=stop(ii,1) )
                                mask(start(ii):stop(ii),1) = ones(stop(ii)-start(ii)+1,1);          % mask the best cluster of RBs for the UE
                                UE_ind = schedule_RME==UE_i;
                                schedule_RME(UE_ind) = schedule_RME(UE_ind).*mask(UE_ind);
                            end
                        end        
                    end
                end

                % recursive expansion
                if sum(schedule_RME==0)
                    % search for highest SNR among unscheduled RBs
                    unassigned_mean_SNR_RB_temp = mean_SNR_RB_temp .* repmat(schedule_RME==0,[1 obj.nUEs]);
                    [max_unlocated,ind_unlocated] = max(unassigned_mean_SNR_RB_temp(:));
                    [RB_unloc,UE_unloc] = ind2sub(size(unassigned_mean_SNR_RB_temp),ind_unlocated);

                    % search upper neighbor
                    neighbor = false;
                    if(RB_unloc~=1)
                        for ii=0:obj.RB_grid_size
                            if(RB_unloc-ii==0)
                                break;
                            end
                            if(schedule_RME(RB_unloc-ii,1)~=0 && schedule_RME(RB_unloc-ii,1)==UE_unloc)
                                neighbor = true;
                                break;
                            elseif(schedule_RME(RB_unloc-ii,1)~=0)
                                break;                    
                            end
                        end
                    end

                    % search lower neighbor
                    if ~neighbor
                        if(RB_unloc~=obj.RB_grid_size)
                            for ii=0:obj.RB_grid_size
                                if(RB_unloc+ii==obj.RB_grid_size+1)
                                    break;
                                end
                                if(schedule_RME(RB_unloc+ii,1)~=0 && schedule_RME(RB_unloc+ii,1)==UE_unloc)
                                    neighbor = true;
                                    break;
                                elseif(schedule_RME(RB_unloc+ii,1)~=0)
                                    break; 
                                end
                            end
                        end
                    end

                    if ~neighbor
                        mean_SNR_RB_temp(RB_unloc,UE_unloc)=0;
                    else
                        clusters    = diff([1;schedule_RME;1]==0);
                        start       = find(clusters==1);
                        stop        = find(clusters==-1)-1;

                        mask = zeros(obj.RB_grid_size,1);
                        for ii=1:length(start)                                                      % look for the cluster where the UE has the best SNR
                            if( RB_unloc>=start(ii,1) && RB_unloc<=stop(ii,1) )
                                mask(start(ii):stop(ii),1) = ones(stop(ii)-start(ii)+1,1);          % mask the best cluster of RBs for the UE
                                unloc_ind = schedule_RME==0;
                                schedule_RME(unloc_ind) = schedule_RME(unloc_ind) + UE_unloc.*mask(unloc_ind);
                            end
                        end  
                    end
                else
                    run = 0;
                end
            end
            
            % assign the best schedule to the UEs
            for u_=1:obj.nUEs
                obj.UE_scheduling(u_).UE_mapping    = repmat(schedule_RME,[1 2])==u_;
                obj.UE_scheduling(u_).assigned_RBs  = squeeze(sum(sum(obj.UE_scheduling(u_).UE_mapping,1),2));
            end

            % calculate link apadtion
            if obj.UE_config.mode == 1    % assume all UEs have same transmission mode                
                UE_scheduling = obj.link_adaptation_SIxO(obj.UE_scheduling,BS,cqi_i,channel,sigma_n2,channel_MSE); 
            elseif obj.UE_config.mode == 4    
                UE_scheduling = obj.link_adaptation_CLSM(obj.UE_scheduling,channel,cqi_i,sigma_n2,channel_MSE);
            end

            % calculate the number of bits
            obj.calculate_allocated_bits(UE_scheduling,subframe_corr,SyncUsedElements,PBCHsyms);
            
            for uu = 1:obj.nUEs
%                 UE_scheduling(uu).total_no_refsym   = total_no_refsym;
                UE_scheduling(uu).SyncUsedElements  = SyncUsedElements;
                UE_scheduling(uu).CHusedElements    = PBCHsyms;
            end
        end
    end
end
