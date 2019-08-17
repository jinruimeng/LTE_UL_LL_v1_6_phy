classdef antenna
    % This Abstract class represents an antenna
    % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016

    properties
        antenna_type
        max_antenna_gain
        pattern_is_3D = false;
    end

    methods (Abstract)
        % Print some info
        print(obj)
        % Returns antenna gain as a function of theta
        antenna_gain = gain(obj,theta)
        % Returns the maximum and minimum antenna gain [min max]
        minmaxgain = min_max_gain(obj)
    end
    
    methods (Static)
        function an_eNodeB = attach_antenna_to_eNodeB(an_eNodeB,LTE_config) 
            switch an_eNodeB.antenna_type 
                case 'TR36.873 3D antenna'
                    an_eNodeB.antenna = antennas.TR36873_3DAntenna(LTE_config);    
                    % Additional parameters 
                    an_eNodeB.tx_height              = LTE_config.tx_height;
                case 'TR36.873 3D antenna omnidirectional'
                    an_eNodeB.antenna = antennas.TR36873_3DAntenna_omnidir(LTE_config);    
                    % Additional parameters 
                    an_eNodeB.tx_height              = LTE_config.tx_height;
                otherwise
                    error('This antenna is not supported');
            end
        end
    end
end
