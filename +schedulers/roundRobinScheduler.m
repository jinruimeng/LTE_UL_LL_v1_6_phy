classdef roundRobinScheduler < network_elements.lteScheduler
    % Round Robin scheduler that assigns the same amount of RBs to each UE
    % Stefan Pratschner
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at
    
    properties
        CQI_mapping_data
        alphabets
    end
    
    methods
        function obj = roundRobinScheduler(params,UE,averager,mapping_data,varargin)
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
            
            % generate RB allocation
            UE_mapping_temp = zeros(obj.nUEs,1);
            RBs_to_schedule = obj.RB_grid_size;
            UEs_to_schedule = obj.nUEs;
            u_ = 1;            
            while RBs_to_schedule
                if ~mod(RBs_to_schedule,UEs_to_schedule)          % if the number of UEs divides the number of RBs
                    UE_mapping_temp(u_:obj.nUEs) = RBs_to_schedule/UEs_to_schedule .* ones(UEs_to_schedule,1);
                    RBs_to_schedule = 0;
                else
                    RBs_temp = round(obj.RB_grid_size/obj.nUEs);
                    UE_mapping_temp(u_,1) = RBs_temp;
                    RBs_to_schedule = RBs_to_schedule - RBs_temp;
                    UEs_to_schedule = UEs_to_schedule - 1;
                end
                u_ = u_ + 1;
            end
            
            % randomize the UE assignement
            UE_mapping_temp = UE_mapping_temp(randperm(size(UE_mapping_temp,1)),:);
            
            % change format of RB allocation
            UE_mapping_temp_big = [];
            for uu=1:obj.nUEs
                UE_mapping_temp_big = [UE_mapping_temp_big; uu.*ones(UE_mapping_temp(uu,1),2)];
            end
                        
            % assign the RBs to the users
            for u_=1:obj.nUEs
                obj.UE_scheduling(u_).UE_mapping = UE_mapping_temp_big==u_;
                obj.UE_scheduling(u_).assigned_RBs = squeeze(sum(sum(obj.UE_scheduling(u_).UE_mapping,1),2));
            end

            % calculate link adaption
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
