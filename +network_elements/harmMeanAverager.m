classdef harmMeanAverager < network_elements.sinrAverager
    % Implements an harmonic mean averager.
    
    properties
        % BICM data
        MI_data
    end
    
    
    
    methods
        
        %class constructor
        function obj = harmMeanAverager(varargin)
            obj.MI_data = network_elements.MIdata_load;         %This is used in the lteScheduler.m
        end
        
        function effective_SINR = average(obj,SINR_vector,varargin)
            
            CQIs=varargin{1};
            effective_SINR = 10*log10(1./mean(1./SINR_vector))*ones(length(CQIs),1);
        end
        
    end
    
end


