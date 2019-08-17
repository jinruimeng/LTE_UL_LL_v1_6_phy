classdef miscUtils
    % Implements miscellaneous functions.
    % Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
    % (c) 2011 by INTHFT
    % www.nt.tuwien.ac.at
    
    properties
    end
    
    methods(Static)

        function mod_angle = wrapTo359(the_angle)
            % Equivalent to wrapTo360, but it maps to the interval [0 359] (360 is set to zero).
            mod_angle = mod(the_angle, 360);
        end
    end
end

