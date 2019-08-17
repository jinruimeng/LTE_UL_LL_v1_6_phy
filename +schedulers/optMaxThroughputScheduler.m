classdef optMaxThroughputScheduler < network_elements.lteScheduler
    % Scheduler that maximizes the thorughput over all possible schedules by exhausitve search
    % Stefan Pratschner, spratsch@nt.tuwien.ac.at
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at
    
    properties
        CQI_mapping_data
        alphabets
        schedules           % all possible LTE uplink schedules
        N_schedules         % number of possible schedules
    end
    
    methods
        function obj = optMaxThroughputScheduler(params,UE,averager,mapping_data,varargin)
            % Fill in basic parameters (handled by the superclass constructor)
            obj = obj@network_elements.lteScheduler(params,UE,averager,mapping_data);
            obj.CQI_mapping_data = mapping_data;
            if ~isempty(varargin)
                obj.alphabets = varargin{1};
            else
                obj.alphabets = [];
            end
            
            % generate schedules
            obj.schedules = obj.generate_schedules;
            obj.N_schedules = size(obj.schedules,2);
            
            % init UE_scheduling
            obj.UE_scheduling = obj.get_new_UE_params;
        end
        
        
        function UE_scheduling = scheduler_users(obj,subframe_corr,SyncUsedElements,BS,PBCHsyms,cqi_i,ChanMod_output,sigma_n2,channel_MSE)
            % resource allocation for all UEs
            
            % initialization
            best_rate       = 0;
            best_schedule   = zeros(obj.RB_grid_size,1);

            % exhausitve search for all possible schedules
            for sched_i = 1:obj.N_schedules
                temp_schedule = obj.schedules(:,sched_i);
                
                % assign RBs to the users
                for u_=1:obj.nUEs
                    obj.UE_scheduling(u_).UE_mapping    = repmat(temp_schedule,[1 2])==u_;
                    obj.UE_scheduling(u_).assigned_RBs  = squeeze(sum(sum(obj.UE_scheduling(u_).UE_mapping,1),2));
                end

                % calculate link apadtion
                if obj.UE_config.mode == 1    % assume all UEs have same transmission mode                
                    UE_scheduling = obj.link_adaptation_SIxO(obj.UE_scheduling,BS,cqi_i,ChanMod_output,sigma_n2,channel_MSE); 
                elseif obj.UE_config.mode == 4    
                    UE_scheduling = obj.link_adaptation_CLSM(obj.UE_scheduling,ChanMod_output,cqi_i,sigma_n2,channel_MSE);
                end
                
                % calculate sum rate of all UEs
                temp_rate = 0;
                for uu = 1:obj.nUEs
                    if UE_scheduling(uu).assigned_RBs
                        temp_rate = temp_rate + UE_scheduling(uu).CQI_params.efficiency.*UE_scheduling(uu).assigned_RBs;
                    end
                end                

                % take best schedule
                if(temp_rate > best_rate)
                    best_rate = temp_rate;
                    best_schedule = temp_schedule;
                end
            end
            
            % assign the best schedule to the UEs
            for u_=1:obj.nUEs
                obj.UE_scheduling(u_).UE_mapping    = repmat(best_schedule,[1 2])==u_;
                obj.UE_scheduling(u_).assigned_RBs  = squeeze(sum(sum(obj.UE_scheduling(u_).UE_mapping,1),2));
            end

            % calculate link apadtion
            if obj.UE_config.mode == 1    % assume all UEs have same transmission mode                
                UE_scheduling = obj.link_adaptation_SIxO(obj.UE_scheduling,BS,cqi_i,ChanMod_output,sigma_n2,channel_MSE); 
            elseif obj.UE_config.mode == 4    
                UE_scheduling = obj.link_adaptation_CLSM(obj.UE_scheduling,ChanMod_output,cqi_i,sigma_n2,channel_MSE);
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
