classdef randomMUMIMO < network_elements.lteScheduler
    % Scheduler that randomly picks subset of all the users and schedules them for one subframe 
    % Lukas Nagel, lnagel@nt.tuwien.ac.at
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at
    
    properties
        CQI_mapping_data
        alphabets
    end
    
    methods
        function obj = randomMUMIMO(params,UE,averager,mapping_data,varargin)
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
            % assumption N_ue >= N_rx otherwise it is not random because we
            % just schedule all users
            

            N_rx    = obj.BS_config.nRX;             
            N_ue_total    = obj.nUEs;
            
            if N_ue_total >= N_rx
                % generate all possible schedules
                possible_schedules = combnk(1:N_ue_total, N_rx);

                % randomly pick one of them
                random_schedule = possible_schedules(randi([1, size(possible_schedules,1)]),:);
            else
                % otherwise it is pretty pointless and not "random"
                random_schedule = 1:N_ue_total;
            end
          
            % apply this schedule
            for u_= 1:N_ue_total
                if any(u_ == random_schedule)
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
