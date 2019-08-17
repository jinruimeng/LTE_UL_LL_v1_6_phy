
    %% Generate the LOS/NLOS probability
  
    function is_LOS = generate_LOS_probability(distance, UE_is_high, scenario)
    % Generates the probability of being LOS (Line of sight)
    % Used for link level simulation when no pathloss-and shadowign maps
    % are generated
    % input:  distance        ... actual 2D distance between transmitter and
    %         receiver in [m]
    %         UE_is_heigh     ... actual UE height in [m]
    %         scenario        ... actual simulated scenario: 3D_UMa or 3D_UMi
    
    % output: is_LOS          ... boolean type: 1->LOS; 0->NLOS
    % Refference: Table 7.2-2 TR 36.873
    % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016
     
    if distance > 18
        g = (1.25e-6).*distance.^3.*exp(-distance./150);
    else
        g = 0;
    end
    
    if UE_is_high < 13
        C_variable = 0;
    else
        C_variable = ((UE_is_high - 13)/10)^(1.5)*g;
    end
    
    %
    switch scenario
        case '3D_UMi_fading'
            probability_LOS = (min(18./distance,1).*(1-exp(-distance./63))+exp(-distance./63)).*(1+C_variable);
        case '3D_UMa_fading'
            probability_LOS = min(18./distance,1).*(1-exp(-distance./36))+exp(-distance./36);
        otherwise
            error('Environment scenario not valid')
    end
    %
    is_LOS = binornd(1,probability_LOS);
    
    end
