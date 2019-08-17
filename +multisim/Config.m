classdef Config < handle
    %Config class contains non default parameters of a simulation
    
    properties
        label
        parameters
        cqi
    end
    
    methods
        function SC = Config(label)
            SC.label = label;
            SC.parameters = containers.Map;
            SC.cqi = NaN;
        end
        
        function add(SC, name, value)
                SC.parameters(name) = value;
        end 
        
        function setCQI(SC, cqi)
            SC.cqi = cqi;
        end
    end
    
end

