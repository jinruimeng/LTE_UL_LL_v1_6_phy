classdef fullbuffer < traffic_models.generic_tm
% This class is used for full buffer simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC
properties
end

methods
    function obj = fullbuffer
        obj = obj@traffic_models.generic_tm;
        obj.data_buffer = Inf;
        obj.type = 'fullbuffer';
    end
 
    function check_deadline(~,~)
    end
end

end