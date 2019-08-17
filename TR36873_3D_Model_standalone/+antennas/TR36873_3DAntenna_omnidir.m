classdef TR36873_3DAntenna_omnidir < antennas.TR36873_3DAntenna
    methods
        function obj = TR36873_3DAntenna_omnidir(LTE_config)
            obj@antennas.TR36873_3DAntenna(LTE_config);
            obj.antenna_type = 'TR36.873 3D antenna omnidirectional';
        end
        
    end
    methods(Static)
        % Returns antenna gain in [dB] => 0
        function antenna_gain = gain(~, ~, ~)
            antenna_gain = 0; % omnidirectional
        end
    end
end