classdef fixed_MU_MIMO < network_elements.lteScheduler
    % Scheduler that assigns a fixed and predefined amount of RBs to each user
    % Stefan Pratschner, spratsch@nt.tuwien.ac.at
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at
    
    properties
        CQI_mapping_data
        alphabets
    end
    
    methods
        function obj = fixed_MU_MIMO(params,UE,averager,mapping_data,varargin)
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
            
            % check if the fixed assigned RBs make sense
%             if( obj.nUEs > obj.BS_config.nRX)
%                 error('there are more UEs than receive antennas');
%             end
                
            % assign the RBs to the users
            for u_=1:obj.nUEs
                obj.UE_scheduling(u_).UE_mapping    = logical(ones(obj.RB_grid_size,2));
                
                obj.UE_scheduling(u_).assigned_RBs  = squeeze(sum(sum(obj.UE_scheduling(u_).UE_mapping,1),2));
            end

            % calculate link adaption
            if obj.UE_config.mode == 1    % assume all UEs have same transmission mode                
                UE_scheduling = obj.link_adaptation_SIxO(obj.UE_scheduling,BS,cqi_i,ChanMod_output,sigma_n2,channel_MSE); 
            elseif obj.UE_config.mode == 4    
                UE_scheduling = obj.link_adaptation_CLSM(obj.UE_scheduling,ChanMod_output,cqi_i,sigma_n2,channel_MSE);
            elseif obj.UE_config.mode == 5
                UE_scheduling = obj.link_adaptation_MU_MIMO(obj.UE_scheduling,BS,cqi_i,ChanMod_output,sigma_n2,channel_MSE);
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
