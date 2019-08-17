classdef TR36873PathlossModel < macroscopic_pathloss_models.generalPathlossModel
    % Propagation conditions as proposed by TR 36.873 V12.0.0
    % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016
  
    properties
        frequency               % Frequency in HERTZs (for consistency, although this 
                                % model makes calculations with frequencies set in GHz
        environment             % Environment that this instance represents
                       
        W_street                % Street width in m - The applicability range according to TR 36.873 is 20m < W < 50m
        
        B_height                % Avg. building height - The applicability range according to TR 36.873 is 20m < B < 50m
        
        h_effective             % Effective environment height
        
        light_speed             % Speed of light 299792458 m/s
        pathloss_scnenario      % UMa or UMi scenario
    end

    methods
        % Class constructor
        function obj = TR36873PathlossModel(frequency,environment)
            obj.frequency   = frequency;
            obj.environment = environment;
            obj.name        = 'TR 36.873';
            obj.light_speed = 299792458;
            switch environment
                case '3D_UMa'
                    obj.name = [obj.name ' 3D_UMa'];
                    obj.W_street    = 20;
                    obj.B_height    = 20;
                    obj.pathloss_scnenario = 1;
                case '3D_UMi'
                    obj.name = [obj.name ' 3D_UMi'];  
                    obj.pathloss_scnenario = 2;
                otherwise
                    error(['"' environment '"" environment not valid']);
            end
        end
               
        
        
        % Returns the macroscopic pathloss in dB. Note: distance in METERS
        % Calculates pathloss in dB combining LOS and NLOS pathlosses using
        % the condition specified with function LOS_positions  
        %Changed from 'TR36873PathlossModel' lines 46 to 58
         function [pathloss_in_dB] = pathloss(obj,distance,config,networkMacroscopicPathlossMap, pathloss_LOS_positions)
             % pathloss   -  Returns the macroscopic pathloss in dB
             % Calculates the pathloss in dB by combining LOS (Line of sight)
             % and NLOS (Non-Line of sight) pathloss using the propagation
             % condition specified with the the parameter
             % pathloss_LOS_positions.
             % If the UE is indoors, the pathloss is generated based on
             % O-to-I scenario 
             %
             % input:    distance               ... actual 2D distance in [m]
             %           pathloss_LOS_positions ... actual propagation condition [boolean type]
             % output:   pathloss_in_dB         ... macroscopic pathloss in dB
              % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016

             wall_pathloss = 20;             % Loss thorugh wall           
             indoor_pathloss = 0.5.*networkMacroscopicPathlossMap.UE_indoor_map;   % Penetration loss, where distance_indoor is assumed uniformly distributed betweeen 0 and 25 [m]
             distance_3D = sqrt(distance.^2 + (config.tx_height - networkMacroscopicPathlossMap.UE_height_map).^2);
             pathloss_in_dB_LOS = obj.pathloss_LOS(distance,config,networkMacroscopicPathlossMap, distance_3D);
             pathloss_in_dB_NLOS = obj.pathloss_NLOS(distance,config,networkMacroscopicPathlossMap);
             pathloss_in_dB_NLOS = max(pathloss_in_dB_NLOS, pathloss_in_dB_LOS);
           % Pathloss for both LOS and NLOS parts of the link
             pathloss_LOS_positions(isnan(pathloss_LOS_positions)) = 0;
             pathloss_in_dB = networkMacroscopicPathlossMap.UE_indoor_map.*((pathloss_LOS_positions.*pathloss_in_dB_LOS + ~pathloss_LOS_positions.*pathloss_in_dB_NLOS) + wall_pathloss + indoor_pathloss)  +...
                            ~networkMacroscopicPathlossMap.UE_indoor_map.*(pathloss_LOS_positions.*pathloss_in_dB_LOS + ~pathloss_LOS_positions.*pathloss_in_dB_NLOS);
         end
                
        % Generate distance in 3D: 
        function distance_3D = calculate_3D_distance(obj,distance,LTE_config,networkMacroscopicPathlossMap)
            % Generates the 3D distance between transmitter and receiver
            % Reffers to the d_3D in 3GPP TR 36.873 
            UE_height_map = networkMacroscopicPathlossMap.UE_height_map;
            distance_3D = sqrt(distance.^2 + (LTE_config.tx_height - UE_height_map).^2);
        end
        
        % Calculate the variable C as specified in TR 36.873 that is used
        % to determine LOS probability and environment probability
        %Changed from 'TR36873PathlossModel' lines 65 to 77
        function C_variable = calculate_probability(obj,distance,networkMacroscopicPathlossMap)
            % Calculates the variable C used to determine the probability
            % of being LOS. 
            % Refference: Table 7.2-2, TR 36.873
            UE_height_map = networkMacroscopicPathlossMap.UE_height_map;
            compare_distance = distance > 18;
            distance_compared = distance.*compare_distance;
            g = (1.25e-6).*distance_compared.^3.*exp(-distance_compared./150);
            
            UE_is_high     = (UE_height_map > 13) & (UE_height_map < 23);
            UE_is_high(isnan(UE_is_high))=0;
            C_variable_all = ((UE_height_map - 13)./10).^(1.5).*g;
            C_variable     = zeros(size(UE_height_map)) + C_variable_all.*UE_is_high;
               
        end
        
        %Changed from 'TR36873PathlossModel' line 80
        function environment_probability = effective_probability(obj,distance,networkMacroscopicPathlossMap)
            % Determines the environment probability used to generate the
            % effective environmet height h_E
            % environment_probability = 1/(1+C(d_2D,UE_height))
            % The variable C is derived from calculate_probability
            environment_probability = 1./(1 + obj.calculate_probability(max(distance-networkMacroscopicPathlossMap.UE_indoor_map,0),networkMacroscopicPathlossMap)); % consider distance outdoor in case the UE is indoor
        end       

        % Probability of the link to be LOS as specified in TR 36.873
        function p_LOS = LOS_probability(obj,distance,networkMacroscopicPathlossMap)
            % Generates the probability of being LOS for 3D_UMa and 3D_UMi
            % scenarios
            % Reference: Table 7.2-2, TR 36.873 
            switch obj.pathloss_scnenario
                case 1
                    p_LOS = (min(18./distance,1).*(1-exp(-distance./63))+exp(-distance./63)).*(1+obj.calculate_probability(distance,networkMacroscopicPathlossMap));
                case 2
                    p_LOS = min(18./distance,1).*(1-exp(-distance./36))+exp(-distance./36);
                otherwise
                    error(['"' obj.environment '"" environment not valid']);
            end
        end
 
        function LOS_positions = generate_LOS_positions(obj,distance,networkMacroscopicPathlossMap) 
            % Creates LOS/NLOS positions pixel-wise based on the 
            % probability of being LOS 
            %
            % output:      LOS_positions ... boolean type 
            LOS_positions = binornd(1,obj.LOS_probability(distance,networkMacroscopicPathlossMap));
            LOS_positions (isnan(LOS_positions))=0; 
        end
        
        function pl  = pathloss_LOS(obj,distance,LTE_config,networkMacroscopicPathlossMap, distance_3D) 
            % Urban macro cell LOS (Line of sight) scenario
            % with high UE density according to TR 36.873
            % where BS is above sorrounding buildings.
            % input:    distance ... actual 2D distance in [m]
            %           distance_3D ... actual 3D distance in [m]
            %
            % output:   pl       ... LOS pathloss in dB
            switch obj.pathloss_scnenario
                case 1
                    rand_var = binornd(1,obj.effective_probability(distance,networkMacroscopicPathlossMap));
                    if rand_var == 1
                        obj.h_effective = 1;
                    else
                        obj.h_effective = randsrc(1,1,[12,15,18,21]);  % Uniform distribution uniform(12,15,18,21), heights in [m]
                    end
                case 2
                    obj.h_effective = 1;
                otherwise
                    error(['"' obj.environment '"" environment not valid']);
            end
            
            tx_height_effective = LTE_config.tx_height - obj.h_effective;
            rx_height_effective = networkMacroscopicPathlossMap.UE_height_map - obj.h_effective;
            d_break_point       = 4*tx_height_effective*rx_height_effective*obj.frequency/obj.light_speed;  
            frequency           = obj.frequency/1000000000;    % Calculations are done in freq in MHz 

            distance_smaller_than_breakpoint = (distance > 10) & (distance <= d_break_point);
            pl_lower_distance = 22.0.*log10(distance_3D) + 28.0 + 20.*log10(frequency);
            pl_higher_distance = 40.*log10(distance_3D)...
                                        + 28 + 20.*log10(frequency) - 9.*log10((d_break_point).^2 + (LTE_config.tx_height - networkMacroscopicPathlossMap.UE_height_map).^2);            
            pl = pl_lower_distance.*distance_smaller_than_breakpoint +  pl_higher_distance.*~distance_smaller_than_breakpoint;                 
        end

        % NLOS pathloss
        %Changed from 'TR36873PathlossModel' lines 102 to 132 
        %changing the distance_3D
        function pl_NLOS = pathloss_NLOS(obj,distance,config,networkMacroscopicPathlossMap)
            % Urban macro cell NLOS (Non-Line of sight) scenario
            % with high UE density according to TR 36.873
            % input:    distance    ... actual 2D distance in [m]
            %           distance_3D ... actual 3D distance in [m]
            %
            % output:   pl_NLOS     ... NLOS pathloss in dB

            frequency   = obj.frequency/1000000000;                                       % Calculations are done in freq in MHz
            switch obj.pathloss_scnenario
                case 1
                    pl_NLOS = 161.04 - 7.1*log10(obj.W_street) + 7.5*log10(obj.B_height) - (24.37 - 3.7*(obj.B_height./config.tx_height)^2)...
                        *log10(config.tx_height)+(43.42 - 3.1*log10(config.tx_height))*(log10(obj.calculate_3D_distance(distance,config,networkMacroscopicPathlossMap)) - 3) + 20*log10(frequency)...
                        -(3.2*(log10(17.625))^2 - 4.97) - 0.6*(networkMacroscopicPathlossMap.UE_height_map - 1.5);
                case 2
                    pl_NLOS = 36.7*log10(distance) +22.7 + 26*log10(frequency) - 0.3*(networkMacroscopicPathlossMap.UE_height_map - 1.5);
                otherwise
                    error(['"' obj.environment '"" environment not valid']);
            end
            
           % pl          = max(pl_NLOS,obj.pathloss_LOS(distance,LTE_config,networkMacroscopicPathlossMap));
        end
        %Not in the simulator
        function pl_O_to_I = pathloss_O_to_I(obj,distance,LTE_config,networkMacroscopicPathlossMap)
            % Urban macro cell O-to-I (Outdoor-to-Indoor) scenario
            % with high UE density according to TR 36.873
            % input:    distance   ... actual 2D distance in [m]
            %
            % output:   pl_O_to_I  ... NLOS pathloss in dB
            pl_wall = 20;             % Loss thorugh wall
            distance_indoor = unifrnd(0,25);    % Indoor distance assumed uniformly distributed between 0 and 25 m
            pl_indoor = 0.5.*distance_indoor; % LOss inside , where distance_indoor is assumed uniformly distributed betweeen 0 and 25 [m]
            pl_basic = obj.pathloss(distance,LTE_config,networkMacroscopicPathlossMap).*(obj.calculate_3D_distance(distance,LTE_config,networkMacroscopicPathlossMap)); % Basic pathloss as a multiplication of pathloss outdoor and distnces_3D (indoor and outdoor)
            pl_O_to_I = pl_basic + pl_wall + pl_indoor;
        end        

            
    end
end
        


        
       
        
        


       
      