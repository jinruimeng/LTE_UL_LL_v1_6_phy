classdef bestCQIScheduler < network_elements.lteScheduler
    % Scheduler that assignes all RBs in one subframe to the UE with the
    % best CQI value
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
        function obj = bestCQIScheduler(params,UE,averager,mapping_data,varargin)
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
            % resource allocation for all UEs
            
            % initialization
            best_CQI    = 0;
            best_ue     = 0;

            % exhausitve search for all possible schedules
            for sched_ue = 1:obj.nUEs
                
                temp_schedule = sched_ue*ones(obj.RB_grid_size,2);
                % assign RBs to the users
                for u_=1:obj.nUEs
                    obj.UE_scheduling(u_).UE_mapping    = temp_schedule==u_;
                    obj.UE_scheduling(u_).assigned_RBs  = squeeze(sum(sum(obj.UE_scheduling(u_).UE_mapping,1),2));
                end

                % calculate link apadtion
                if obj.UE_config.mode == 1    % assume all UEs have same transmission mode                
                    UE_scheduling = obj.link_adaptation_SIxO(obj.UE_scheduling,BS,cqi_i,ChanMod_output,sigma_n2,channel_MSE); 
                elseif obj.UE_config.mode == 4    
                    UE_scheduling = obj.link_adaptation_CLSM(obj.UE_scheduling,ChanMod_output,cqi_i,sigma_n2,channel_MSE);
                end             

                % take best schedule
                if(UE_scheduling(sched_ue).cqi > best_CQI)
                    best_CQI = UE_scheduling(sched_ue).cqi;
                    best_ue = sched_ue;
                end
            end
            
            % assign the best schedule to the UEs
            for u_=1:obj.nUEs
                if(u_ == best_ue)
                    obj.UE_scheduling(u_).UE_mapping    = true(obj.RB_grid_size,2);
                else
                    obj.UE_scheduling(u_).UE_mapping    = false(obj.RB_grid_size,2);
                end
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
