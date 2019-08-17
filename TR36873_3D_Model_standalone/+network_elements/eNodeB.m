classdef eNodeB < handle
    % This object represents the eNodeB sectors with properties:
    %  id                          ... sector id: 1, 2 or 3
    %  eNodeB_id                   ... eNodeB sector id depending on the parent eNodeB id: form 1 to 21
    %  parent_eNodeB               ... parent eNodeB site
    %  neighbors_eNodeB            ... a list of neighboring eNodebs
    %  pos                         ... eNodeB position
    %  max_power                   ... maximum transmit power
    %  nTX                         ... number of antenna ports
    %  tx_height                   ... height in [m]
    %  boresight                   ... boresight angle of antenna pointing in [°]
    %  attached_UEs                ... number of attached UEs
    %  attached_UEs_vector         ... a list of attached UEs
    %  clock                       ... the world-wide clock of the simulation
    %  macroscopic_pathloss_model  ... a list of macroscopic path loss
    %  antenna_type                ... antenna type 
    %  antenna                     ... alist of antenna type properties
    %  electrical_downtilt         ... electrical downtilt steering angle of antenna element
    %  mechanical_downtilt         ... mechanichal downtilt angle of antenna element
    %
    % NOTE: The correlation of Large Scale Parameters as give in 3GPP TR 36.873 
    % is generated for each eNodeB.
    % Assumption from 3GPP TR 36.873: correlation applies to the UEs served by the
    % same eNodeB. No correlation exists between UEs served by different
    % eNodeBs
    
    % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016
    
    properties
        id
        eNodeB_id
        parent_eNodeB
        neighbors_eNodeB
        pos
        max_power
        nTX     
        tx_height
        boresight               
        attached_UEs = 0;            % Number of attached UEs
        attached_UEs_vector          % A list of attached UEs
        clock
        macroscopic_pathloss_model
        antenna_type
        antenna
        electrical_downtilt
        mechanical_downtilt
    end
    
    methods
        
        % Attachs a user to this eNodeB, first checking that the node is
        % not already in the list. It will update the UE's
        % 'attached_site' variable, effectively binding the UE to this
        % eNodeB. 
        function attachUser(obj,user)
            if isempty(obj.attached_UEs_vector)
                % If the user list is empty
                obj.attached_UEs_vector  = user;
                obj.attached_UEs         = obj.attached_UEs + 1;
            else
                % If there are already some users, check if the UE is new
                % to the list. If yes, add him.
                current_UEs = [obj.attached_UEs_vector.id];
                if ~sum(current_UEs==user.id)
                    obj.attached_UEs_vector = [obj.attached_UEs_vector user];
                    obj.attached_UEs        = obj.attached_UEs + 1;
                end
            end
            
            % Fill eNodeB-attachment info from the UE
            user.attached_sector_idx = obj.id;
            user.attached_site       = obj.parent_eNodeB;
            user.attached_eNodeB     = obj;  
        end
        

        
                % Calculate correlation of large scale parameters for LOS case
        function correlated_large_scale_param_LOS = calculate_correlated_large_scale_param_LOS(obj,LTE_config, large_scale_params_LOS)
            % Cross-correlation generation %
            % Extract cross correlation parameters from input
            a = LTE_config.asD_ds_LOS;
            b = LTE_config.asA_ds_LOS;
            c = LTE_config.asA_sf_LOS;
            d = LTE_config.asD_sf_LOS;
            e = LTE_config.ds_sf_LOS;
            f = LTE_config.asD_asA_LOS;
            g = LTE_config.asD_kf_LOS;
            h = LTE_config.asA_kf_LOS;
            k = LTE_config.ds_kf_LOS;
            l = LTE_config.sf_kf_LOS;
            m = LTE_config.zsD_sf_LOS;
            n = LTE_config.zsA_sf_LOS;
            o = LTE_config.zsD_kf_LOS;
            p = LTE_config.zsA_kf_LOS;
            q = LTE_config.zsD_ds_LOS;
            r = LTE_config.zsA_ds_LOS;
            s = LTE_config.zsD_asD_LOS;
            t = LTE_config.zsA_asD_LOS;
            u = LTE_config.zsD_asA_LOS;
            v = LTE_config.zsA_asA_LOS;
            w = LTE_config.zsD_zsA_LOS;
            
              R_LOS = [1 l e d c m n;...
                   l 1 k g h o p;...
                   e k 1 a b q r;...
                   d g a 1 f s t;...
                   c h b f 1 u v;...
                   m o q s u 1 w;...
                   n p r t v w 1];

%                R_LOS_sqrt = chol(R_LOS);
               R_LOS_sqrt = sqrtm(R_LOS);
               correlated_large_scale_param_LOS = R_LOS_sqrt*large_scale_params_LOS;
        end
        
        

        
        
        %Calculate correlation of large scale parameters for NLOS case
        function correlated_large_scale_param_NLOS = calculate_correlated_large_scale_param_NLOS(obj,LTE_config, large_scale_params_NLOS)
            % Cross-correlation generation %
            % Extract cross correlation parameters from input
            a = LTE_config.asD_ds_NLOS;
            b = LTE_config.asA_ds_NLOS;
            c = LTE_config.asA_sf_NLOS;
            d = LTE_config.asD_sf_NLOS;
            e = LTE_config.ds_sf_NLOS;
            f = LTE_config.asD_asA_NLOS;
%             g = 0;                       % LTE_config.asD_kf_NLOS N/A;
%             h = 0;                       % LTE_config.asA_kf_NLOS N/A;
%             k = 0;                       % LTE_config.ds_kf_NLOS N/A;
%             l = 0;                       % LTE_config.sf_kf_NLOS N/A;
            m = LTE_config.zsD_sf_NLOS;
            n = LTE_config.zsA_sf_NLOS;
%             o = 0;                       % LTE_config.zsD_kf_NLOS N/A;
%             p = 0;                       % LTE_config.zsA_kf_NLOS N/A;
            q = LTE_config.zsD_ds_NLOS;
            r = LTE_config.zsA_ds_NLOS;
            s = LTE_config.zsD_asD_NLOS;
            t = LTE_config.zsA_asD_NLOS;
            u = LTE_config.zsD_asA_NLOS;
            v = LTE_config.zsA_asA_NLOS;
            w = LTE_config.zsD_zsA_NLOS;
            
              R_NLOS = [1 e d c m n;...
                   e 1 a b q r;...
                   d a 1 f s t;...
                   c b f 1 u v;...
                   m q s u 1 w;...
                   n r t v w 1];
               
%                R_NLOS_sqrt = chol(R_NLOS);
                R_NLOS_sqrt = sqrtm(R_NLOS);
               correlated_large_scale_param_NLOS = R_NLOS_sqrt*large_scale_params_NLOS;     
        end

        %Added from 'eNodeB_sector.m' lines 648 to 683
        %Calculate correlation of large scale parameters for OTOI case
        function correlated_large_scale_param_OTOI = calculate_correlated_large_scale_param_OTOI(obj,LTE_config, large_scale_params_OTOI)
            % Extract cross correlation parameters from input
            a = LTE_config.asD_ds_OTOI;
            b = LTE_config.asA_ds_OTOI;
            c = LTE_config.asA_sf_OTOI;
            d = LTE_config.asD_sf_OTOI;
            e = LTE_config.ds_sf_OTOI;
            f = LTE_config.asD_asA_OTOI;
%             g = 0;                       % LTE_config.asD_kf_NLOS N/A;
%             h = 0;                       % LTE_config.asA_kf_NLOS N/A;
%             k = 0;                       % LTE_config.ds_kf_NLOS N/A;
%             l = 0;                       % LTE_config.sf_kf_NLOS N/A;
            m = LTE_config.zsD_sf_OTOI;
            n = LTE_config.zsA_sf_OTOI;
%             o = 0;                       % LTE_config.zsD_kf_NLOS N/A;
%             p = 0;                       % LTE_config.zsA_kf_NLOS N/A;
            q = LTE_config.zsD_ds_OTOI;
            r = LTE_config.zsA_ds_OTOI;
            s = LTE_config.zsD_asD_OTOI;
            t = LTE_config.zsA_asD_OTOI;
            u = LTE_config.zsD_asA_OTOI;
            v = LTE_config.zsA_asA_OTOI;
            w = LTE_config.zsD_zsA_OTOI;
            
              R_OTOI = [1 e d c m n;...
                   e 1 a b q r;...
                   d a 1 f s t;...
                   c b f 1 u v;...
                   m q s u 1 w;...
                   n r t v w 1];
               
                R_OTOI_sqrt = sqrtm(R_OTOI);
               correlated_large_scale_param_OTOI = R_OTOI_sqrt*large_scale_params_OTOI;     
   
        end
        
        %Transfromation of large scale parameters from normal Gaussian to log-normal
        function [sigmas_SF_LOS, sigmas_KF_LOS, sigmas_DS_LOS, sigmas_ASD_LOS, sigmas_ASA_LOS, sigmas_ZSD_LOS, sigmas_ZSA_LOS] = sigmas_LOS(obj,LTE_config, correlated_large_scale_params,ZOD_parameters)
            DS_mu_LOS             = LTE_config.DS_mu_LOS;
            DS_sigma_LOS          = LTE_config.DS_sigma_LOS;
            AS_D_mu_LOS           = LTE_config.AS_D_mu_LOS;
            AS_D_sigma_LOS        = LTE_config.AS_D_sigma_LOS;
            AS_A_mu_LOS           = LTE_config.AS_A_mu_LOS;
            AS_A_sigma_LOS        = LTE_config.AS_A_sigma_LOS;
            SF_sigma_LOS          = LTE_config.SF_sigma_LOS;
            KF_mu_LOS             = LTE_config.KF_mu_LOS;
            KF_sigma_LOS          = LTE_config.KF_sigma_LOS;
            %
            if LTE_config.generate_uplink
                ZS_A_mu_LOS = ZOD_parameters(1);
                ZS_A_sigma_LOS = ZOD_parameters(2);
                ZS_D_mu_LOS = LTE_config.ZS_A_mu_LOS;
                ZS_D_sigma_LOS = LTE_config.ZS_A_sigma_LOS;
            else
                ZS_A_mu_LOS = LTE_config.ZS_A_mu_LOS;
                ZS_A_sigma_LOS = LTE_config.ZS_A_sigma_LOS;
                ZS_D_mu_LOS = ZOD_parameters(1);
                ZS_D_sigma_LOS = ZOD_parameters(2);
            end
%             [ZS_D_mu_LOS,ZS_D_sigma_LOS,~] = obj.generate_ZSD_ZoD_offset_parameters_LOS(LTE_config, rx_height, distance);
            %Generate sigmas for LOS case
            sigmas_SF_LOS       = correlated_large_scale_params(1).*SF_sigma_LOS;                         %[dB]
            sigmas_KF_LOS       = KF_mu_LOS + correlated_large_scale_params(2).*KF_sigma_LOS;             %[dB]
            sigmas_DS_LOS       = 10.^(DS_mu_LOS + correlated_large_scale_params(3).*DS_sigma_LOS);       %[s]
            sigmas_ASD_LOS      = 10.^(AS_D_mu_LOS + correlated_large_scale_params(4).*AS_D_sigma_LOS);  %[°]
            sigmas_ASA_LOS      = 10.^(AS_A_mu_LOS + correlated_large_scale_params(5).*AS_A_sigma_LOS);  %[°]
            sigmas_ZSA_LOS      = 10.^(ZS_A_mu_LOS + correlated_large_scale_params(7).*ZS_A_sigma_LOS);  %[°]
            %Generate sigma_ZSD using offset parameters in Table 7.3-7 for 3D-UMa
            sigmas_ZSD_LOS      = 10.^(ZS_D_mu_LOS + correlated_large_scale_params(6).*ZS_D_sigma_LOS);  %[°]
            %Limit random rms azimuth and elevation arrival/departure spread values
            sigmas_ASD_LOS      = min(sigmas_ASD_LOS,104);             %[°]
            sigmas_ASA_LOS      = min(sigmas_ASA_LOS,104);             %[°]
            sigmas_ZSD_LOS      = min(sigmas_ZSD_LOS,52);              %[°]
            sigmas_ZSA_LOS      = min(sigmas_ZSA_LOS,52);              %[°]
        end
        
        %ADDED non_value_param
        function [sigmas_SF_NLOS,non_value_param, sigmas_DS_NLOS, sigmas_ASD_NLOS, sigmas_ASA_NLOS, sigmas_ZSD_NLOS, sigmas_ZSA_NLOS] = sigmas_NLOS(obj,LTE_config, correlated_large_scale_params, ZOD_parameters)
            DS_mu_NLOS            = LTE_config.DS_mu_NLOS;
            DS_sigma_NLOS         = LTE_config.DS_sigma_NLOS;
            AS_D_mu_NLOS          = LTE_config.AS_D_mu_NLOS;
            AS_D_sigma_NLOS       = LTE_config.AS_D_sigma_NLOS;
            AS_A_mu_NLOS          = LTE_config.AS_A_mu_NLOS;
            AS_A_sigma_NLOS       = LTE_config.AS_A_sigma_NLOS;
            SF_sigma_NLOS         = LTE_config.SF_sigma_NLOS;
            %
            
            if LTE_config.generate_uplink
                ZS_A_mu_NLOS = ZOD_parameters(1);
                ZS_A_sigma_NLOS = ZOD_parameters(2);
                ZS_D_mu_NLOS = LTE_config.ZS_A_mu_NLOS;
                ZS_D_sigma_NLOS = LTE_config.ZS_A_sigma_NLOS;
            else
                ZS_A_mu_NLOS = LTE_config.ZS_A_mu_NLOS;
                ZS_A_sigma_NLOS = LTE_config.ZS_A_sigma_NLOS;
                ZS_D_mu_NLOS = ZOD_parameters(1);
                ZS_D_sigma_NLOS = ZOD_parameters(2);
            end

%             [ZS_D_mu_NLOS,ZS_D_sigma_NLOS,~] = obj.generate_ZSD_ZoD_offset_parameters_NLOS(LTE_config,rx_height, distance);
            %Generate sigmas for NLOS case
            sigmas_SF_NLOS        = correlated_large_scale_params(1).*SF_sigma_NLOS;                          %[dB]
            sigmas_DS_NLOS        = 10^(DS_mu_NLOS + correlated_large_scale_params(2).*DS_sigma_NLOS);       %[s]
            sigmas_ASD_NLOS       = 10.^(AS_D_mu_NLOS + correlated_large_scale_params(3).*AS_D_sigma_NLOS);  %[°]
            sigmas_ASA_NLOS       = 10.^(AS_A_mu_NLOS + correlated_large_scale_params(4).*AS_A_sigma_NLOS);  %[°]
            sigmas_ZSA_NLOS       = 10.^(ZS_A_mu_NLOS + correlated_large_scale_params(6).*ZS_A_sigma_NLOS);  %[°]
            %Generate sigma_ZSD using offset parameters in Table 7.3-7 for 3D-UMa
            sigmas_ZSD_NLOS       = 10.^(ZS_D_mu_NLOS + correlated_large_scale_params(5).*ZS_D_sigma_NLOS);  %[°]
            %Limit random rms azimuth and elevation arrival/departure spread values
            sigmas_ASD_NLOS       = min(sigmas_ASD_NLOS,104);             %[°]
            sigmas_ASA_NLOS       = min(sigmas_ASA_NLOS,104);             %[°]
            sigmas_ZSD_NLOS       = min(sigmas_ZSD_NLOS,52);              %[°]
            sigmas_ZSA_NLOS       = min(sigmas_ZSA_NLOS,52);              %[°]
            non_value_param       =0; %ADDED
        end
        
        %Added from 'eNodeB_sector.m' lines 751 to 778
        function [sigmas_SF_OTOI, non_value_param, sigmas_DS_OTOI, sigmas_ASD_OTOI, sigmas_ASA_OTOI, sigmas_ZSD_OTOI, sigmas_ZSA_OTOI] = sigmas_OTOI(obj,LTE_config, correlated_large_scale_params,ZOD_parameters)
            DS_mu_OTOI            = LTE_config.DS_mu_OTOI;
            DS_sigma_OTOI         = LTE_config.DS_sigma_OTOI;
            AS_D_mu_OTOI         = LTE_config.AS_D_mu_OTOI;
            AS_D_sigma_OTOI       = LTE_config.AS_D_sigma_OTOI;
            AS_A_mu_OTOI          = LTE_config.AS_A_mu_OTOI;
            AS_A_sigma_OTOI       = LTE_config.AS_A_sigma_OTOI;

            SF_sigma_OTOI         = LTE_config.SF_sigma_OTOI;
            %changed from the simulator
            
            if LTE_config.generate_uplink
                ZS_A_mu_OTOI = ZOD_parameters(1);
                ZS_A_sigma_OTOI = ZOD_parameters(2);
                ZS_D_mu_OTOI = LTE_config.ZS_A_mu_OTOI;
                ZS_D_sigma_OTOI = LTE_config.ZS_A_sigma_OTOI;
            else
                ZS_A_mu_OTOI = LTE_config.ZS_A_mu_OTOI;
                ZS_A_sigma_OTOI = LTE_config.ZS_A_sigma_OTOI;
                ZS_D_mu_OTOI = ZOD_parameters(1);
                ZS_D_sigma_OTOI = ZOD_parameters(2);
            end

%             [ZS_D_mu_OTOI,ZS_D_sigma_OTOI,~] = obj.generate_ZSD_ZoD_offset_parameters_NLOS(LTE_config,rx_height, distance);
%             
            sigmas_SF_OTOI        = correlated_large_scale_params(1).*SF_sigma_OTOI;                          %[dB]
            sigmas_DS_OTOI        = 10^(DS_mu_OTOI + correlated_large_scale_params(2).*DS_sigma_OTOI);       %[s]
            sigmas_ASD_OTOI       = 10.^(AS_D_mu_OTOI + correlated_large_scale_params(3).*AS_D_sigma_OTOI);  %[°]
            sigmas_ASA_OTOI       = 10.^(AS_A_mu_OTOI+ correlated_large_scale_params(4).*AS_A_sigma_OTOI);  %[°]
            sigmas_ZSA_OTOI       = 10.^(ZS_A_mu_OTOI + correlated_large_scale_params(6).*ZS_A_sigma_OTOI);  %[°]
            %Generate sigma_ZSD using offset parameters in Table 7.3-7 for 3D-UMa
            sigmas_ZSD_OTOI       = 10.^(ZS_D_mu_OTOI + correlated_large_scale_params(5).*ZS_D_sigma_OTOI);  %[°]
            %Limit random rms azimuth and elevation arrival/departure spread values
            sigmas_ASD_OTOI       = min(sigmas_ASD_OTOI,104);             %[°]
            sigmas_ASA_OTOI       = min(sigmas_ASA_OTOI,104);             %[°]
            sigmas_ZSD_OTOI       = min(sigmas_ZSD_OTOI,52);              %[°]
            sigmas_ZSA_OTOI       = min(sigmas_ZSA_OTOI,52);              %[°]
            non_value_param       = 0;
        end
        
        function [ZS_D_mu_NLOS,ZS_D_sigma_NLOS,ZOD_mu_offset_NLOS] = generate_ZSD_ZoD_offset_parameters_NLOS(obj, LTE_config, rx_height, distance)
            %3D-UMa
            if strcmp(LTE_config.channel_model.type,'3D_UMa_fading')
                ZS_D_mu_NLOS          = max(-0.5,(-2.1.*(distance./1000)-0.01.*(rx_height-1.5)+0.9));        % zenith arrival angle spread, mean [log10(deg)]
                ZS_D_sigma_NLOS       = 0.49;                                                                % zenith arrival angle spread, std [log10(deg)]
                ZOD_mu_offset_NLOS    = -10.^(-0.62*log10(max(10,distance))+1.93-0.07*(rx_height-1.5));
            end
            %3D_UMi
            if strcmp(LTE_config.channel_model.type,'3D_UMi_fading')
                ZS_D_mu_NLOS          = max(-0.5,(-2.1.*(distance./1000)+0.01*max(rx_height-LTE_config.tx_height,0)+0.9));     % zenith arrival angle spread, mean [log10(deg)]
                ZS_D_sigma_NLOS       = 0.6;                                                                 % zenith arrival angle spread, std [log10(deg)]
                ZOD_mu_offset_NLOS    = -10.^(-0.55*log10(max(10,distance))+1.6);
            end
        end
        
        %Added from 'eNodeB_sector.m' lines 797 to 810
        function [ZS_D_mu_OTOI,ZS_D_sigma_OTOI,ZOD_mu_offset_OTOI] = generate_ZSD_ZoD_offset_parameters_OTOI_NLOS(obj,LTE_config, rx_height, distance)
            %3D-UMa
            if strcmp(LTE_config.channel_model.type,'3D_UMa_fading')
                ZS_D_mu_OTOI          = max(-0.5,(-2.1.*(distance./1000)-0.01.*(rx_height-1.5)+0.9));        % zenith arrival angle spread, mean [log10(deg)]
                ZS_D_sigma_OTOI       = 0.49;                                                                % zenith arrival angle spread, std [log10(deg)]
                ZOD_mu_offset_OTOI    = -10.^(-0.62*log10(max(10,distance))+1.93-0.07*(rx_height-1.5));
            end
            %3D_UMi
            if strcmp(LTE_config.channel_model.type,'3D_UMi_fading')
                ZS_D_mu_OTOI          = max(-0.5, -2.1.*(distance./1000)+0.01*max(rx_height-LTE_config.tx_height,0)+0.9);     % zenith arrival angle spread, mean [log10(deg)]
                ZS_D_sigma_OTOI       = 0.6;                                                                 % zenith arrival angle spread, std [log10(deg)]
                ZOD_mu_offset_OTOI    = -10.^(-0.55*log10(max(10,distance))+1.6);
            end
        end
        
        %Added from 'eNodeB_sector.m' lines 812 to 825
        function [ZS_D_mu_OTOI,ZS_D_sigma_OTOI,ZOD_mu_offset_OTOI] = generate_ZSD_ZoD_offset_parameters_OTOI_LOS(obj,LTE_config, rx_height, distance)
            %3D-UMa
            if strcmp(LTE_config.channel_model.type,'3D_UMa_fading')
                ZS_D_mu_OTOI          = max(-0.5,(-2.1.*(distance./1000)-0.01.*(rx_height-1.5)+0.75));        % zenith arrival angle spread, mean [log10(deg)]
                ZS_D_sigma_OTOI       = 0.4;                                                                % zenith arrival angle spread, std [log10(deg)]
                ZOD_mu_offset_OTOI    = 0;
            end
            %3D_UMi
            if strcmp(LTE_config.channel_model.type,'3D_UMi_fading')
                ZS_D_mu_OTOI          = max(-0.5, -2.1.*(distance./1000)+0.01*abs(rx_height-LTE_config.tx_height)+0.75);     % zenith arrival angle spread, mean [log10(deg)]
                ZS_D_sigma_OTOI       = 0.4;                                                                 % zenith arrival angle spread, std [log10(deg)]
                ZOD_mu_offset_OTOI    = 0;
            end
        end
        
        function [ZS_D_mu_LOS,ZS_D_sigma_LOS,ZOD_mu_offset_LOS] = generate_ZSD_ZoD_offset_parameters_LOS(obj, LTE_config, rx_height, distance)
            %3D_UMa
            if strcmp(LTE_config.channel_model.type,'3D_UMa_fading')
                ZS_D_mu_LOS         = max((-0.5),(-2.1.*(distance./1000)-0.01.*(rx_height-1.5)+0.75));    % zenith arrival angle spread, mean [log10(deg)]
                ZS_D_sigma_LOS      = 0.40;                                                               % zenith arrival angle spread, std [log10(deg)]
                ZOD_mu_offset_LOS       = 0;
            end
            %3D_UMi
            if strcmp(LTE_config.channel_model.type,'3D_UMi_fading')
                ZS_D_mu_LOS         = max((-0.5),(-2.1.*(distance./1000)+0.01*abs(rx_height-LTE_config.tx_height)+0.75));   % zenith arrival angle spread, mean [log10(deg)]
                ZS_D_sigma_LOS      = 0.4;                                                                % zenith arrival angle spread, std [log10(deg)]
                ZOD_mu_offset_LOS       = 0;
            end
        end

    end
    
end

