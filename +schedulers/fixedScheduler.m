classdef fixedScheduler < network_elements.lteScheduler
    % Scheduler that assigns a fixed and predefined amount of RBs to each user
    % Stefan Pratschner, spratsch@nt.tuwien.ac.at
    % (c) 2016 by ITC
    % www.nt.tuwien.ac.at
    
    properties
        CQI_mapping_data
        alphabets
    end
    
    methods
        function obj = fixedScheduler(params,UE,averager,mapping_data,varargin)
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
            if( ~(size(obj.fixed_assignment,1)==1 && size(obj.fixed_assignment,2)==obj.nUEs) || ...
                    max(obj.fixed_assignment)>obj.RB_grid_size || ...
                    min(obj.fixed_assignment)<0 || ...
                    sum(obj.fixed_assignment)~=obj.RB_grid_size)
                error('the specified fixed RB assignement is not correct');
            end
            
            % generate RB allocation from fixed assignment
            UE_mapping_temp = [];
            for uu=1:obj.nUEs
                UE_mapping_temp = [UE_mapping_temp; uu.*ones(obj.fixed_assignment(uu),2)];
            end
                        
            % assign the fixed RBs to the users
            for u_=1:obj.nUEs
                obj.UE_scheduling(u_).UE_mapping = UE_mapping_temp==u_;
                obj.UE_scheduling(u_).assigned_RBs = squeeze(sum(sum(obj.UE_scheduling(u_).UE_mapping,1),2));
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
