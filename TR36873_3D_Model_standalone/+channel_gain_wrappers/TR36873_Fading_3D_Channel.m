classdef TR36873_Fading_3D_Channel < handle
    % Class that represents fast fading calculation and channel matrix
    % generation,as specified in TR 36.873.
    % Includes Step 5 - Step 11 of the channel coefficient generation procedure
    %
    % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016

    properties
        % Common parameters
        channel_type
        generate_uplink
        % Scenario specific parameters
        C_constant
        C_constant_elevation
        PerClusterRays
        NumClusters_LOS
        r_DS_LOS
        KF_mu_LOS
        PerClusterAS_D_LOS
        PerClusterAS_A_LOS
        PerClusterZS_A_LOS
        LNS_ksi_LOS
        xpr_mu_LOS
        xpr_sigma_LOS
        cluster_ASD_LOS
        cluster_ASA_LOS
        cluster_ZSA_LOS
        %%%%%%%%%%%%%%%%%%
        NumClusters_NLOS
        r_DS_NLOS
        PerClusterAS_D_NLOS
        PerClusterAS_A_NLOS
        PerClusterZS_A_NLOS
        LNS_ksi_NLOS
        xpr_mu_NLOS
        xpr_sigma_NLOS
        cluster_ASD_NLOS
        cluster_ASA_NLOS
        cluster_ZSA_NLOS
        offset_vec
        %%%%%%%%%%%%%%%%%%
        NumClusters_OTOI
        r_DS_OTOI
        PerClusterAS_D_OTOI
        PerClusterAS_A_OTOI
        PerClusterZS_A_OTOI
        LNS_ksi_OTOI
        xpr_mu_OTOI
        xpr_sigma_OTOI
        cluster_ASD_OTOI
        cluster_ASA_OTOI
        cluster_ZSA_OTOI
        % Save this parameters from the fast fading
        delay_LOS = [];
        delay_NLOS = [];
        delay_OTOI = [];
        % Outdoor-to-Indoor scenario: OTOI
        cl_power_OTOI                                = [];
        cl_power_OTOI_general                        = [];
        cl_power_OTOI_per_ray                        = [];
        xpr_OTOI                                     = [];
        RX_spherical_unit_vector_OTOI                = [];
        TX_spherical_unit_vector_OTOI                = [];
        Doppler_speed_OTOI                           = [];
        azimuth_angle_of_departure_OTOI_once         = [];
        azimuth_angle_of_arrival_OTOI_once           = [];
        zenith_angle_of_departure_OTOI_once          = [];
        zenith_angle_of_arrival_OTOI_once            = [];
        UE_theta_AntennaField_global_OTOI            = [];
        UE_phi_AntennaField_global_OTOI              = [];
        theta_AntennaField_tx_global_OTOI            = [];
        phi_AntennaField_tx_global_OTOI              = [];
        field_and_phases_over_rays_and_clusters_OTOI = [];
        exponential_part_without_Doppler_OTOI        = [];
        Doppler_over_time_OTOI                       = [];
        channel_coefficients_over_clusters_OTOI      = [];
        int_sampled_channel_OTOI                     = [];
        delay_tap_size_OTOI                          = [];
        chan_norm_OTOI                               = [];
        
        % Save NLOS parameters
        cl_power_NLOS                                = [];
        cl_power_NLOS_general                        = [];
        cl_power_NLOS_per_ray                        = [];
        xpr_NLOS                                     = [];
        RX_spherical_unit_vector_NLOS                = [];
        TX_spherical_unit_vector_NLOS                = [];
        Doppler_speed_NLOS                           = [];
        azimuth_angle_of_departure_NLOS_once         = [];
        azimuth_angle_of_arrival_NLOS_once           = [];
        zenith_angle_of_departure_NLOS_once          = [];
        zenith_angle_of_arrival_NLOS_once            = [];
        UE_theta_AntennaField_global_NLOS            = [];
        UE_phi_AntennaField_global_NLOS              = [];
        theta_AntennaField_tx_global_NLOS            = [];
        phi_AntennaField_tx_global_NLOS              = [];
        field_and_phases_over_rays_and_clusters_NLOS = [];
        exponential_part_without_Doppler_NLOS        = [];
        Doppler_over_time_NLOS                       = [];
        channel_coefficients_over_clusters_NLOS      = [];
        int_sampled_channel_NLOS                     = [];
        delay_tap_size_NLOS                          = [];
        chan_norm_NLOS                               = [];
        
        % Save LOS parameters
        xpr_LOS                                      = [];
        cl_power_LOS                                 = [];
        cl_power_LOS_general                         = [];
        cl_power_LOS_per_ray                         = [];
        ray_power_without_LOS_final                  = [];
        RX_spherical_unit_vector_LOS                 = [];
        TX_spherical_unit_vector_LOS                 = [];
        RX_spherical_unit_vector_LOS_direct_ray      = [];
        TX_spherical_unit_vector_LOS_direct_ray      = [];
        Doppler_speed_LOS                            = [];
        Doppler_speed_LOS_direct_ray                 = [];
        azimuth_angle_of_departure_LOS_once          = [];
        azimuth_angle_of_arrival_LOS_once            = [];
        zenith_angle_of_departure_LOS_once           = [];
        zenith_angle_of_arrival_LOS_once             = [];
        UE_theta_AntennaField_global_LOS             = [];
        UE_phi_AntennaField_global_LOS               = [];
        theta_AntennaField_tx_global_LOS             = [];
        phi_AntennaField_tx_global_LOS               = [];
        UE_theta_AntennaField_global_LOS_direct_ray  = [];
        UE_phi_AntennaField_global_LOS_direct_ray    = [];
        theta_AntennaField_tx_global_LOS_direct_ray  = [];
        phi_AntennaField_tx_global_LOS_direct_ray    = [];
        field_and_phases_single_LOS_ray              = [];
        field_and_phases_over_rays_and_clusters_LOS  = [];
        Doppler_over_time_LOS                        = [];
        Doopler_LOS_direct_ray                       = [];
        exp_part_without_Doppler_LOS                 = [];
        exp_part_without_Doppler_direct_ray          = [];
        channel_coefficients_over_clusters_LOS       = [];
        channel_coefficients_over_clusters_LOS_ray   = [];
        int_sampled_channel_LOS                      = [];
        delay_tap_size_LOS                           = [];  
        chan_norm_LOS                                = [];
        chan_norm_LOS_ray                            = [];
        
        % For uplink
        ZS_A_mu_LOS
        ZS_A_mu_NLOS
        ZS_A_mu_OTOI
        
        %Save direct path angles between Tx and Rx
        theta_arrival              = [];
        theta_departure            = [];
        phi_arrival                = [];
        phi_departure              = [];
        
        BS_boresight;
        
        % Properies to generate PDP
        bandwidth
        subcarrierSpacing = 15e3;  % By default
        resourceBlock     = 180e3; % Fixed badwidth of resource block in Hz, page 33
        Nsc                        % number of subcarriers in one resource block, fixed length of resource block in Hz, page 33
        Nrb                        % number of resource blocks, transmission BW is 90% of the total BW unless for 1.4 MHz
        Ntot                       % Total number of subcarriers not NULL
        Nfft                       % number of FFT points
        Tb                         % useful Symbol Time
        fs                         % Sampling frequency
        FFT_sampling_interval = 6; % Sampling interval in the frequency domain when getting RB FFT traces
        
    end
    
    % Associated eNodeB methods
    methods
        % Class constructor
        function obj = TR36873_Fading_3D_Channel(LTE_config)
            obj.channel_type          = LTE_config.channel_model.type;
            obj.generate_uplink       = LTE_config.generate_uplink;
            % LOS
            obj.C_constant            = LTE_config.C_constant;
            obj.C_constant_elevation  = LTE_config.C_constant_elevation;
            obj.PerClusterRays        = LTE_config.PerClusterRays;
            obj.NumClusters_LOS       = LTE_config.NumClusters_LOS;
            obj.KF_mu_LOS             = LTE_config.KF_mu_LOS;
            obj.r_DS_LOS              = LTE_config.r_DS_LOS;
            obj.PerClusterAS_D_LOS    = LTE_config.PerClusterAS_D_LOS;
            obj.PerClusterAS_A_LOS    = LTE_config.PerClusterAS_A_LOS;
            obj.PerClusterZS_A_LOS    = LTE_config.PerClusterZS_A_LOS;
            obj.LNS_ksi_LOS           = LTE_config.LNS_ksi_LOS;
            obj.xpr_mu_LOS            = LTE_config.xpr_mu_LOS;
            obj.xpr_sigma_LOS         = LTE_config.xpr_sigma_LOS;
            obj.cluster_ASD_LOS       = LTE_config.cluster_ASD_LOS;
            obj.cluster_ASA_LOS       = LTE_config.cluster_ASA_LOS;
            obj.cluster_ZSA_LOS       = LTE_config.cluster_ZSA_LOS;
            % NLOS
            obj.NumClusters_NLOS      = LTE_config.NumClusters_NLOS;
            obj.r_DS_NLOS             = LTE_config.r_DS_NLOS;
            obj.PerClusterAS_D_NLOS   = LTE_config.PerClusterAS_D_NLOS;
            obj.PerClusterAS_A_NLOS   = LTE_config.PerClusterAS_A_NLOS;
            obj.PerClusterZS_A_NLOS   = LTE_config.PerClusterZS_A_NLOS;
            obj.LNS_ksi_NLOS          = LTE_config.LNS_ksi_NLOS;
            obj.xpr_mu_NLOS           = LTE_config.xpr_mu_NLOS;
            obj.xpr_sigma_NLOS        = LTE_config.xpr_sigma_NLOS;
            obj.cluster_ASD_NLOS      = LTE_config.cluster_ASD_NLOS;
            obj.cluster_ASA_NLOS      = LTE_config.cluster_ASA_NLOS;
            obj.cluster_ZSA_NLOS      = LTE_config.cluster_ZSA_NLOS;   
            obj.offset_vec            = [0.0447,0.0447,0.1413,0.1413,0.2492,0.2492,0.3715,0.3715,0.5129,0.5129,0.6797,0.6797,...
                                                     0.8844,0.8844,1.1481,1.1481,1.5195,1.5195,2.1551,2.1551].';
            % OTOI                                     
            obj.NumClusters_OTOI      = LTE_config.NumClusters_OTOI;
            obj.r_DS_OTOI             = LTE_config.r_DS_OTOI;
            obj.PerClusterAS_D_OTOI   = LTE_config.PerClusterAS_D_OTOI;
            obj.PerClusterAS_A_OTOI   = LTE_config.PerClusterAS_A_OTOI;
            obj.PerClusterZS_A_OTOI   = LTE_config.PerClusterZS_A_OTOI;
            obj.LNS_ksi_OTOI          = LTE_config.LNS_ksi_OTOI;
            obj.xpr_mu_OTOI           = LTE_config.xpr_mu_OTOI;
            obj.xpr_sigma_OTOI        = LTE_config.xpr_sigma_OTOI;
            obj.cluster_ASD_OTOI      = LTE_config.cluster_ASD_OTOI;
            obj.cluster_ASA_OTOI      = LTE_config.cluster_ASA_OTOI;
            obj.cluster_ZSA_OTOI      = LTE_config.cluster_ZSA_OTOI;
            
            % For uplink
            obj.ZS_A_mu_LOS = LTE_config.ZS_A_mu_LOS;
            obj.ZS_A_mu_NLOS = LTE_config.ZS_A_mu_NLOS;
            obj.ZS_A_mu_OTOI = LTE_config.ZS_A_mu_OTOI;
            
            obj.bandwidth = LTE_config.bandwidth;
            obj.Nsc       = obj.resourceBlock/obj.subcarrierSpacing;
            %%%%
            if(obj.bandwidth == 1.4e6)
                obj.Nrb = 6;
            else
                obj.Nrb = (obj.bandwidth*0.9) / obj.resourceBlock;
            end
            obj.Ntot = obj.Nsc*obj.Nrb;
            if(obj.bandwidth == 15e6 && obj.subcarrierSpacing == 15e3)
                obj.Nfft = 1536;
            elseif(obj.bandwidth == 15e6 && obj.subcarrierSpacing == 7.5e3)
                obj.Nfft = 1536*2;
            else
                obj.Nfft =  2^ceil(log2(obj.Ntot));
            end
            obj.Tb = 1/obj.subcarrierSpacing;
            obj.fs = obj.subcarrierSpacing*obj.Nfft;                                     
        end
                 
           %% Generate LOS elevation and azimuth angles of arrival/departure 
           % downlink: angles of departure defined at eNodeB, angles of arrival defined at UE
           function [theta_arrival,theta_departure,phi_arrival,phi_departure] = eNodeB_UE_LOS_direction_angles_signal(obj,LTE_config,rx_height,attached_site_pos,UE_pos)   
               % Generates angles of departure- and arrival of the direct
               % link between eNodeB and an outdoor UE 
               %
               % input:     tx_height       ... eNodeB height in [m]
               %            rx_height       ... UE height in [m]
               %            eNodeB_pos      ... actual eNodeB position [x y]
               %            UE_pos          ... actual UE position [x y]
               %
               % output:    theta_arrival   ... angle of arrival in elevation
               %            theta_departure ... angle of departure in elevation
               %            phi_arrival     ... angle of arrival in azimuth
               %            phi_departure   ... angle of departure in azimuth
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               tx_height             = LTE_config.tx_height;
               UE_distance_to_site   = sqrt(sum((attached_site_pos - UE_pos).^2));
               distance_3D           = sqrt(UE_distance_to_site^2 +(rx_height - tx_height).^2);
               theta_arrival         = acosd((tx_height-rx_height)./distance_3D);                                           %[°]
               theta_departure       = 180 - theta_arrival;                                                                 %[°]
               phi_departure         = atan2d(UE_pos(:,2)-attached_site_pos(:,2),UE_pos(:,1)-attached_site_pos(:,1));       %[°]
               phi_arrival           = atan2d(attached_site_pos(:,2)-UE_pos(:,2),attached_site_pos(:,1)-UE_pos(:,1));       %[°]
           end
           
           function [theta_arrival,theta_departure,phi_arrival,phi_departure] = eNodeB_UE_LOS_direction_angles_signal_OTOI(obj,LTE_config,rx_height,attached_site_pos,UE_pos, dist_indoor)    
               % Generates angles of departure- and arrival of the direct
               % link between eNodeB and an indoor UE 
               % NOTE: Used for O-to-I (Outdoor-to-Indoor) scenario 
               %
               % input:     tx_height       ... eNodeB height in [m]
               %            rx_height       ... UE height in [m]
               %            eNodeB_pos      ... actual eNodeB position [x y]
               %            UE_pos          ... actual UE position [x y]
               %
               % output:    theta_arrival   ... angle of arrival in elevation
               %            theta_departure ... angle of departure in elevation
               %            phi_arrival     ... angle of arrival in azimuth
               %            phi_departure   ... angle of departure in azimuth
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               
               
               tx_height             = LTE_config.tx_height;
               UE_distance_to_site   = abs(sqrt(sum((attached_site_pos - UE_pos).^2))- dist_indoor);
               distance_3D           = sqrt(UE_distance_to_site^2 +(rx_height - tx_height).^2);
               theta_arrival         = acosd((tx_height-rx_height)./distance_3D);                                           %[°]
               theta_departure       = 180 - theta_arrival;                                                                 %[°]
               phi_departure         = atan2d(UE_pos(:,2)-attached_site_pos(:,2),UE_pos(:,1)-attached_site_pos(:,1));       %[°]
               phi_arrival           = atan2d(attached_site_pos(:,2)-UE_pos(:,2),attached_site_pos(:,1)-UE_pos(:,1));       %[°]
           end
           

           %% Step 5: Generate delays
           function  delay_LOS = generate_delay_LOS(obj,sigmas_LOS)
               % Generates the delay for LOS (Line of Sight) scenario
               % See Step 5 in TR 36.873
               % input:    sigmas_LOS        ... correlated LSP (Large Scale Parameter) vector:
               %                                 [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %
               % output:   delay_LOS         ... scaled delay for LOS eq(7.5) TR 36.873
               %     
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sigmas_KF_LOS         = obj.KF_mu_LOS;   
               sigmas_DS_LOS         = sigmas_LOS(3);
               param_D               = 0.7705 - 0.0433.*(sigmas_KF_LOS) + 0.0002.*(sigmas_KF_LOS).^2+...
                                         0.000017.*(sigmas_KF_LOS).^3;                                            
               delay_LOS             = -obj.r_DS_LOS.*sigmas_DS_LOS.*log(rand(1,obj.NumClusters_LOS));             % Delay for each cluster Eq(7.2)
               delay_LOS             = sort((delay_LOS - min(delay_LOS)),'ascend');                                % Sort delays, Eq(7.3)
               delay_LOS             = delay_LOS./param_D;                                                         % Scaling LOS delay with paramter D, Eq(7.5)
               obj.delay_LOS = delay_LOS;                                                                  
           end
           
           function delay_NLOS = generate_delay_NLOS(obj,sigmas_NLOS)
               % Generates the delay for NLOS (Non-Line of Sight) scenario
               % See Step 5 in TR 36.873
               % input:    sigmas_NLOS        ... correlated LSP (Large Scale Parameter) vector:
               %                                 [SF, 0, DS, ASD, ASA, ZSD, ZSA]
               %
               % output:   delay_NLOS         ... scaled delay for NLOS eq(7.5) TR 36.873
               %     
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sigmas_DS_NLOS        = sigmas_NLOS(3); 
               delay_NLOS            = -obj.r_DS_NLOS.*sigmas_DS_NLOS.*log(rand(1,obj.NumClusters_NLOS));         % Delay for each cluster, Eq(7.2)
               delay_NLOS            = sort((delay_NLOS - min(delay_NLOS)),'ascend');                             % Sort delays, Eq(7.3)
               obj.delay_NLOS = delay_NLOS; 
           end
           
           function delay_OTOI = generate_delay_OTOI(obj,sigmas_OTOI)
               % Generates the delay for O-to-I (Outdoor-to-Indoor) scenario
               % See Step 5 in TR 36.873
               % input:    sigmas_OTOI        ... correlated LSP (Large Scale Parameter) vector:
               %                                 [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %
               % output:   delay_OTOI         ... scaled delay for O-to-I eq(7.5) TR 36.873
               %     
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               
               sigmas_DS_OTOI        = sigmas_OTOI(3);                                                            %after changing the sorting of LSP vector
               delay_OTOI            = -obj.r_DS_OTOI.*sigmas_DS_OTOI.*log(rand(1,obj.NumClusters_OTOI));         % Delay for each cluster, Eq(7.2)
               delay_OTOI            = sort((delay_OTOI - min(delay_OTOI)),'ascend');                             % Sort delays, Eq(7.3)  
               obj.delay_OTOI = delay_OTOI; 
           end
              
           %% Step 6: Generate cluster powers over number of clusers - LOS case
           function [cl_power_LOS,cl_power_LOS_general,cl_power_LOS_per_ray, ray_power_without_LOS_final] = generate_cluster_power_LOS(obj,sigmas_LOS)
               % Generates the cluster powers for LOS (Line of Sight) scenario
               % See Step 6 in TR 36.873
               % input:    sigmas_LOS        ... correlated LSP (Large Scale Parameter) vector:
               %                                 [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %
               % output:   cl_power_LOS                ... cluster powers for LOS: for each cluster and each ray within the cluster
               %                                           a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                        nCl- number of clusters
               %           cl_power_LOS_general        ... power per cluster as given in eq(7.6) TR 36.873
               %           cl_power_LOS_per_ray        ... cluster powers for each ray without applying any threshold
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               
               sigmas_KF_LOS             = obj.KF_mu_LOS;         
               sigmas_DS_LOS             = sigmas_LOS(3);
               delay_LOS_             = obj.generate_delay_LOS(sigmas_LOS);   
               cl_power_LOS_general      = exp(-delay_LOS_*((obj.r_DS_LOS - 1)/(obj.r_DS_LOS*sigmas_DS_LOS))).*10.^(-(normrnd(0,3,1,obj.NumClusters_LOS))/10); % Eq(7.6)
               single_ray_power          = (10.^(sigmas_KF_LOS/10))./(10.^(sigmas_KF_LOS/10) + 1);          % Power of the single LOS ray
               dirac_vector              = zeros(1,obj.NumClusters_LOS);
               dirac_vector(:,1)         = 1;
               
               cl_power_without_LOS      = (1./(10.^(sigmas_KF_LOS/10) + 1)).*(cl_power_LOS_general./sum(cl_power_LOS_general));

               ray_power_without_LOS     = cl_power_without_LOS./obj.PerClusterRays;
               ray_power_without_LOS_final = repmat(ray_power_without_LOS,obj.PerClusterRays,1); 
               cl_power_LOS              = (1./(10.^(sigmas_KF_LOS/10) + 1)).*(cl_power_LOS_general./sum(cl_power_LOS_general))+...
                                              dirac_vector.*single_ray_power;                        % Averaged cluster power for LOS 
               ray_power                 = cl_power_LOS_general./obj.PerClusterRays;
               cl_power_LOS_per_ray      = repmat(ray_power,obj.PerClusterRays,1);                   % Power of each ray within a cluster                        
               max_cluster_power_LOS     = max(cl_power_LOS);
               threshold_LOS             = max_cluster_power_LOS -max_cluster_power_LOS;             % No threshold is used, can be adapted                                          
               cl_power_LOS              = cl_power_LOS(cl_power_LOS > threshold_LOS);  
           end
            
           %Generate cluster powers for nr of clusers - NLOS case
           function [cl_power_NLOS,cl_power_NLOS_general,cl_power_NLOS_per_ray] = generate_cluster_power_NLOS(obj,sigmas_NLOS)
               % Generates the cluster powers for NLOS (Non-Line of Sight) scenario
               % See Step 6 in TR 36.873
               % input:    sigmas_NLOS        ... correlated LSP (Large Scale Parameter) vector:
               %                                 [SF, 0, DS, ASD, ASA, ZSD, ZSA]
               %
               % output:   cl_power_NLOS                ... cluster powers for NLOS: for each cluster and each ray within the cluster
               %                                           a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                        nCl- number of clusters
               %           cl_power_NLOS_general        ... power per cluster as given in eq(7.6) TR 36.873
               %           cl_power_NLOS_per_ray        ... cluster powers for each ray without applying any threshold
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sigmas_DS_NLOS            = sigmas_NLOS(3);
               cl_power_NLOS_general     = exp(-obj.generate_delay_NLOS(sigmas_NLOS)*((obj.r_DS_NLOS-1)/(obj.r_DS_NLOS*sigmas_DS_NLOS))).*...
                                            10.^(-(normrnd(0,3,1,obj.NumClusters_NLOS))/10);                    % Cluster power determined with exponential delay distribution
               cl_power_NLOS_general     = cl_power_NLOS_general./sum(cl_power_NLOS_general);                   % Averaged cluster power for NLOS
               cl_power_NLOS_general_per_ray     = cl_power_NLOS_general./obj.PerClusterRays; 
               cl_power_NLOS_per_ray     = repmat(cl_power_NLOS_general_per_ray,obj.PerClusterRays,1);          % Power assigned to each ray within a cluster
               max_cluster_power_NLOS    = max(cl_power_NLOS_general);
               threshold_NLOS            = max_cluster_power_NLOS - max_cluster_power_NLOS;                     % No threshold is used  -(Remove clusters with less than -25 dB power compared to max cl. power)
               cl_power_NLOS             = cl_power_NLOS_general(cl_power_NLOS_general > threshold_NLOS);        
           end
           
           function [cl_power_OTOI,cl_power_OTOI_general,cl_power_OTOI_per_ray] = generate_cluster_power_OTOI(obj,sigmas_OTOI)
               % Generates the cluster powers for O-to-I (Outdoor-to-Indoor) scenario
               % See Step 6 in TR 36.873
               % input:    sigmas_NLOS        ... correlated LSP (Large Scale Parameter) vector:
               %                                 [SF, 0, DS, ASD, ASA, ZSD, ZSA]
               %
               % output:   cl_power_OTOI                ... cluster powers for O-to-I: for each cluster and each ray within the cluster
               %                                           a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                        nCl- number of clusters
               %           cl_power_OTOI_general        ... power per cluster as given in eq(7.6) TR 36.873
               %           cl_power_OTOI_per_ray        ... cluster powers for each ray without applying any threshold
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               
               sigmas_DS_OTOI            = sigmas_OTOI(3);
               cl_power_OTOI_general     = exp(-obj.generate_delay_OTOI(sigmas_OTOI)*((obj.r_DS_OTOI-1)/(obj.r_DS_OTOI*sigmas_DS_OTOI))).*...
                                            10.^(-(normrnd(0,3,1,obj.NumClusters_OTOI))/10);                    % Cluster power determined with exponential delay distribution
               cl_power_OTOI_general     = cl_power_OTOI_general./sum(cl_power_OTOI_general);                   % Averaged cluster power for NLOS
               cl_power_OTOI_general_per_ray     = cl_power_OTOI_general./obj.PerClusterRays; 
               cl_power_OTOI_per_ray     = repmat(cl_power_OTOI_general_per_ray,obj.PerClusterRays,1);          % Power assigned to each ray within a cluster
               max_cluster_power_OTOI    = max(cl_power_OTOI_general);
               threshold_OTOI            = max_cluster_power_OTOI - max_cluster_power_OTOI;                     % No threshold is used -(Remove clusters with less than -25 dB power compared to max cl. power)
               cl_power_OTOI             = cl_power_OTOI_general(cl_power_OTOI_general > threshold_OTOI);        
           end
           
           %% Step 7: Generate arrival and departure angles for both azimuth and elevation in [°] 
           % Generate Azimuth-of-Arrival (AOA), LOS
           function azimuth_angle_of_arrival_LOS = azimuth_angle_of_arrival_LOS(obj,sigmas_LOS,cl_power_LOS,pp)
               % Generates LOS (Line of Sight) angles of arrival in azimuth
               % See Step 7 in TR 36.873
               % input:    sigmas_LOS                    ... correlated LSP (Large Scale Parameter) vector:
               %                                             [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %           cl_power_LOS                  ... cluster powers 
               %           pp                            ... a counter used when interfering links channel impulse response is generated
               %
               % output:   azimuth_angle_of_arrival_LOS  ... azimuth angle of arrival in LOS: for each cluster and each ray within the cluster
               %                                              a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                           nCl- number of clusters
               %
               % NOTE: After that some clusters are filtered (-25 dB less power
               % than max power in cluster array),it is expected that the number 
               % of clusters is reduced. In Table 7.3-2 are given scaling factors
               % depending on nr of clusters, but not all cluster numbers are
               % considered. For simplicity, the scaling factor will be based 
               % on nr of clusters before filtering as sugested in TR.36.873
               % C_const_LOS finds out this scaling factor
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               
               sigmas_KF_LOS                         = obj.KF_mu_LOS;    
               sigmas_ASA_LOS                        = sigmas_LOS(5);   
               idx                                   = (obj.C_constant(1,:) == obj.NumClusters_LOS);
               C_const_LOS                           = obj.C_constant(2,idx);
               C_LOS                                 = C_const_LOS.*(1.1035-0.028.*sigmas_KF_LOS -0.002.*sigmas_KF_LOS.^2 + 0.0001.*sigmas_KF_LOS.^3);     % Scaling factor for LOS
               azimuth_angle_of_arrival_LOS          = (2.*(sigmas_ASA_LOS./1.4).*sqrt(-log(cl_power_LOS./max(cl_power_LOS))))./C_LOS;         
               rand_factor_add                       = (2.*round(rand(1,length(azimuth_angle_of_arrival_LOS)))-1).*azimuth_angle_of_arrival_LOS ...
                                                             + normrnd(0,sigmas_ASA_LOS/7,1,length(azimuth_angle_of_arrival_LOS));             % Assign a positive or negative sign to the angles and add 
                                                                                                                                               % component or random normal distributiom normrnd
               azimuth_angle_of_arrival_LOS          = rand_factor_add - (rand_factor_add(1) - obj.phi_arrival(pp));                           % Eq (7.11)
               azimuth_angle_of_arrival_LOS_per_ray  = repmat(azimuth_angle_of_arrival_LOS,obj.PerClusterRays,1);                              % Assign AOA per each ray in each cluster
               sign_offset                           = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))  = -1;
               cl_wise_alpha_LOS                     = obj.cluster_ASA_LOS.*obj.offset_vec.*sign_offset;
               azimuth_angle_of_arrival_LOS          = azimuth_angle_of_arrival_LOS_per_ray + repmat(cl_wise_alpha_LOS,1,length(azimuth_angle_of_arrival_LOS));   % Add offset angles per ray in each cluster
               azimuth_angle_of_arrival_LOS          = wrapTo180(azimuth_angle_of_arrival_LOS);
           end
           
           
           % Generate Azimuth-of-Arrival (AOA), NLOS
           function azimuth_angle_of_arrival_NLOS = azimuth_angle_of_arrival_NLOS(obj,sigmas_NLOS,cl_power_NLOS,pp)
               % Generates NLOS (Non-Line of Sight) angles of arrival in azimuth
               % See Step 7 in TR 36.873
               % input:    sigmas_NLOS                    ... correlated LSP (Large Scale Parameter) vector:
               %                                             [SF, 0, DS, ASD, ASA, ZSD, ZSA]
               %           cl_power_NLOS                  ... cluster powers 
               %           pp                             ... a counter used when interfering links channel impulse response is generated
               %
               % output:   azimuth_angle_of_arrival_NLOS  ... azimuth angle of arrival in NLOS: for each cluster and each ray within the cluster
               %                                              a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                           nCl- number of clusters
               %
               % NOTE: After that some clusters are filtered (-25 dB less power
               % than max power in cluster array),it is expected that the number 
               % of clusters is reduced. In Table 7.3-2 are given scaling factors
               % depending on nr of clusters, but not all cluster numbers are
               % considered. For simplicity, the scaling factor will be based 
               % on nr of clusters before filtering as sugested in TR.36.873
               % C_const finds out this scaling factor
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
                sigmas_ASA_NLOS                       = sigmas_NLOS(5);
               idx                                   = (obj.C_constant(1,:)== obj.NumClusters_NLOS);
               C_const                               = obj.C_constant(2,idx);
               azimuth_angle_of_arrival_NLOS         = (2.*(sigmas_ASA_NLOS./1.4).*sqrt(-log(cl_power_NLOS./max(cl_power_NLOS))))./C_const;        
               azimuth_angle_of_arrival_NLOS         = (2.*round(rand(1,length(azimuth_angle_of_arrival_NLOS)))-1).*azimuth_angle_of_arrival_NLOS...
                                                         + normrnd(0,sigmas_ASA_NLOS/7,1,length(azimuth_angle_of_arrival_NLOS))+ obj.phi_arrival(pp); 
               %Assign AOA per each ray in each cluster
               azimuth_angle_of_arrival_NLOS_per_ray = repmat(azimuth_angle_of_arrival_NLOS,obj.PerClusterRays,1);                                                                                  
               sign_offset                           = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))  = -1;
               cl_wise_alpha_NLOS                    = obj.cluster_ASA_NLOS.*obj.offset_vec.*sign_offset;
               azimuth_angle_of_arrival_NLOS         = azimuth_angle_of_arrival_NLOS_per_ray + repmat(cl_wise_alpha_NLOS,1,length(azimuth_angle_of_arrival_NLOS)); % Add offset angles per ray in each cluster
               azimuth_angle_of_arrival_NLOS         = wrapTo180(azimuth_angle_of_arrival_NLOS);
           end
           
           function azimuth_angle_of_arrival_OTOI = azimuth_angle_of_arrival_OTOI(obj,sigmas_OTOI,cl_power_OTOI,pp)
               % Generates angles of arrival in azimuth for O-to-I (Outdoor-to-Indoor) scenario
               % See Step 7 in TR 36.873
               % input:    sigmas_OTOI                    ... correlated LSP (Large Scale Parameter) vector:
               %                                             [SF, 0, DS, ASD, ASA, ZSD, ZSA]
               %           cl_power_OTOI                  ... cluster powers 
               %           pp                             ... a counter used when interfering links channel impulse response is generated
               %
               % output:   azimuth_angle_of_arrival_OTOI  ... azimuth angle of arrival in O-to-I: for each cluster and each ray within the cluster
               %                                              a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                           nCl- number of clusters
               %
               % NOTE: After that some clusters are filtered (-25 dB less power
               % than max power in cluster array),it is expected that the number 
               % of clusters is reduced. In Table 7.3-2 are given scaling factors
               % depending on nr of clusters, but not all cluster numbers are
               % considered. For simplicity, the scaling factor will be based 
               % on nr of clusters before filtering as sugested in TR.36.873
               % C_const finds out this scaling factor
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sigmas_ASA_OTOI                       = sigmas_OTOI(5);
               idx                                   = (obj.C_constant(1,:)== obj.NumClusters_OTOI);
               C_const                               = obj.C_constant(2,idx);
               azimuth_angle_of_arrival_OTOI         = (2.*(sigmas_ASA_OTOI./1.4).*sqrt(-log(cl_power_OTOI./max(cl_power_OTOI))))./C_const;        
               azimuth_angle_of_arrival_OTOI         = (2.*round(rand(1,length(azimuth_angle_of_arrival_OTOI)))-1).*azimuth_angle_of_arrival_OTOI...
                                                         + normrnd(0,sigmas_ASA_OTOI/7,1,length(azimuth_angle_of_arrival_OTOI))+ obj.phi_arrival(pp); 
               %Assign AOA per each ray in each cluster
               azimuth_angle_of_arrival_OTOI_per_ray = repmat(azimuth_angle_of_arrival_OTOI,obj.PerClusterRays,1);                                                                                  
               sign_offset                           = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))  = -1;
               cl_wise_alpha_OTOI                    = obj.cluster_ASA_OTOI.*obj.offset_vec.*sign_offset;
               azimuth_angle_of_arrival_OTOI         = azimuth_angle_of_arrival_OTOI_per_ray + repmat(cl_wise_alpha_OTOI,1,length(azimuth_angle_of_arrival_OTOI)); % Add offset angles per ray in each cluster
               azimuth_angle_of_arrival_OTOI         = wrapTo180(azimuth_angle_of_arrival_OTOI);
           end
           
           function azimuth_angle_of_departure_LOS = azimuth_angle_of_departure_LOS(obj,sigmas_LOS,cl_power_LOS,pp)
               % Generates LOS (Line of Sight) angle of departure in azimuth
               % See Step 7 in TR 36.873
               % input:    sigmas_LOS                       ... correlated LSP (Large Scale Parameter) vector:
               %                                             [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %           cl_power_LOS                     ... cluster powers 
               %           pp                               ... a counter used when interfering links channel impulse response is generated
               %
               % output:   azimuth_angle_of_departure_LOS   ... azimuth angle of departure in LOS: for each cluster and each ray within the cluster
               %                                               a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                           nCl- number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               
               sigmas_KF_LOS                          = obj.KF_mu_LOS; 
               sigmas_ASD_LOS                         = sigmas_LOS(4);
               idx                                    = (obj.C_constant(1,:) == obj.NumClusters_LOS);
               C_const_LOS                            = obj.C_constant(2,idx);
               C_LOS                                  = C_const_LOS.*(1.1035-0.028.*sigmas_KF_LOS -0.002.*sigmas_KF_LOS.^2 + 0.0001.*sigmas_KF_LOS.^3);
               azimuth_angle_of_departure_LOS         = (2.*(sigmas_ASD_LOS./1.4).*sqrt(-log(cl_power_LOS./max(cl_power_LOS))))./C_LOS;  
               rand_factor_add_AOD                    =  (2.*round(rand(1,length(azimuth_angle_of_departure_LOS)))-1).*azimuth_angle_of_departure_LOS ...
                                                       + normrnd(0,sigmas_ASD_LOS/7,1,length(azimuth_angle_of_departure_LOS));         % Assign a positive or negative sign to the angles and add 
                                                                                                                                       % component or random normal distributiom 
               azimuth_angle_of_departure_LOS         = rand_factor_add_AOD - rand_factor_add_AOD(1) + obj.phi_departure(pp);                                            
               %Assign AOD per each ray in each cluster
               azimuth_angle_of_departure_LOS_per_ray = repmat(azimuth_angle_of_departure_LOS,obj.PerClusterRays,1);            
               sign_offset                            = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))   = -1;
               cl_wise_alpha_LOS_AOD                  = obj.cluster_ASD_LOS.*obj.offset_vec.*sign_offset;
               azimuth_angle_of_departure_LOS         = azimuth_angle_of_departure_LOS_per_ray + repmat(cl_wise_alpha_LOS_AOD,1,length(azimuth_angle_of_departure_LOS)); 
               azimuth_angle_of_departure_LOS         = wrapTo180(azimuth_angle_of_departure_LOS);
           end
           
           
           %Generate Azimuth-of-Departure (AOD), NLOS
           function azimuth_angle_of_departure_NLOS = azimuth_angle_of_departure_NLOS(obj,sigmas_NLOS,cl_power_NLOS,pp)
               % Generates NLOS (Non-Line of Sight) angle of departure in azimuth
               % See Step 7 in TR 36.873
               % input:    sigmas_NLOS                      ... correlated LSP (Large Scale Parameter) vector:
               %                                               [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %           cl_power_NLOS                    ... cluster powers 
               %           pp                               ... a counter used when interfering links channel impulse response is generated
               %
               % output:   azimuth_angle_of_departure_NLOS  ... azimuth angle of departure in NLOS: for each cluster and each ray within the cluster
               %                                               a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                           nCl- number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               
               sigmas_ASD_NLOS                         = sigmas_NLOS(4);
               idx                                     = (obj.C_constant(1,:) == obj.NumClusters_NLOS);
               C_const                                 = obj.C_constant(2,idx);
               azimuth_angle_of_departure_NLOS         = (2.*(sigmas_ASD_NLOS./1.4).*sqrt(-log(cl_power_NLOS./max(cl_power_NLOS))))./C_const; 
               azimuth_angle_of_departure_NLOS         = (2.*round(rand(1,length(azimuth_angle_of_departure_NLOS)))-1).*azimuth_angle_of_departure_NLOS...
                                                             + normrnd(0,sigmas_ASD_NLOS/7,1,length(azimuth_angle_of_departure_NLOS))+ obj.phi_departure(pp);  % Assign a positive or negative sign to the angles 
               %Assign the value of AOD for each ray within a cluster
               azimuth_angle_of_departure_NLOS_per_ray = repmat(azimuth_angle_of_departure_NLOS,obj.PerClusterRays,1);
               sign_offset                             = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))    = -1;
               cl_wise_alpha_NLOS_AOD                  = obj.cluster_ASD_NLOS.*obj.offset_vec.*sign_offset;
               azimuth_angle_of_departure_NLOS         = azimuth_angle_of_departure_NLOS_per_ray + repmat(cl_wise_alpha_NLOS_AOD,1,length(azimuth_angle_of_departure_NLOS)); % Add offset angles per ray in each cluster        
               azimuth_angle_of_departure_NLOS         = wrapTo180(azimuth_angle_of_departure_NLOS);
           end
           
            function azimuth_angle_of_departure_OTOI = azimuth_angle_of_departure_OTOI(obj,sigmas_OTOI,cl_power_OTOI,pp)
               % Generates O-to-I (Outdoor-to-Indoor) angle of departure in azimuth
               % See Step 7 in TR 36.873
               % input:    sigmas_OTOI                      ... correlated LSP (Large Scale Parameter) vector:
               %                                               [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %           cl_power_OTOI                    ... cluster powers 
               %           pp                               ... a counter used when interfering links channel impulse response is generated
               %
               % output:   azimuth_angle_of_departure_OTOI  ... azimuth angle of departure in O-to-I scenario: for each cluster and each ray within the cluster
               %                                                a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                             nCl- number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               
               sigmas_ASD_OTOI                         = sigmas_OTOI(4);
               idx                                     = (obj.C_constant(1,:) == obj.NumClusters_OTOI);
               C_const                                 = obj.C_constant(2,idx);
               azimuth_angle_of_departure_OTOI         = (2.*(sigmas_ASD_OTOI./1.4).*sqrt(-log(cl_power_OTOI./max(cl_power_OTOI))))./C_const; 
               azimuth_angle_of_departure_OTOI         = (2.*round(rand(1,length(azimuth_angle_of_departure_OTOI)))-1).*azimuth_angle_of_departure_OTOI...
                                                             + normrnd(0,sigmas_ASD_OTOI/7,1,length(azimuth_angle_of_departure_OTOI))+ obj.phi_departure(pp);  % Assign a positive or negative sign to the angles 
               %Assign the value of AOD for each ray within a cluster
               azimuth_angle_of_departure_OTOI_per_ray = repmat(azimuth_angle_of_departure_OTOI,obj.PerClusterRays,1);
               sign_offset                             = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))    = -1;
               cl_wise_alpha_OTOI_AOD                  = obj.cluster_ASD_OTOI.*obj.offset_vec.*sign_offset;
               azimuth_angle_of_departure_OTOI         = azimuth_angle_of_departure_OTOI_per_ray + repmat(cl_wise_alpha_OTOI_AOD,1,length(azimuth_angle_of_departure_OTOI)); % Add offset angles per ray in each cluster        
               azimuth_angle_of_departure_OTOI         = wrapTo180(azimuth_angle_of_departure_OTOI);
           end
           
           % Generate Zenith-of-Arrival and Zenith-of-Departure angles in [°]
           % Generate Zenith-of-Arrival (ZOA), LOS
           function zenith_angle_of_arrival_LOS = zenith_angle_of_arrival_LOS(obj,sigmas_LOS,cl_power_LOS,pp)
               % Generates LOS (Line of Sight) zenith angle of arrival in [°] 
               % See Step 7 in TR 36.873
               % input:    sigmas_LOS                       ... correlated LSP (Large Scale Parameter) vector:
               %                                             [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %           cl_power_LOS                     ... cluster powers 
               %           pp                               ... a counter used when interfering links channel impulse response is generated
               %
               % output:   zenith_angle_of_arrival_LOS   ... zenith angle of arrival in LOS: for each cluster and each ray within the cluster
               %                                             a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                          nCl- number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sigmas_KF_LOS                        = obj.KF_mu_LOS; 
               sigmas_ZSA_LOS                      = sigmas_LOS(7);
               idx                                 = (obj.C_constant_elevation(1,:) == obj.NumClusters_LOS);
               C_const_el_LOS                      = obj.C_constant_elevation(2,idx); 
               C_LOS_el                            = C_const_el_LOS.*(1.3086+0.0339.*sigmas_KF_LOS -0.0077.*sigmas_KF_LOS.^2 + 0.0002.*sigmas_KF_LOS.^3);               
               zenith_angle_of_arrival_LOS         = (-sigmas_ZSA_LOS.*log(cl_power_LOS./max(cl_power_LOS)))./C_LOS_el;                                    
               rand_factor_add_el                  = (2.*round(rand(1,length(zenith_angle_of_arrival_LOS)))-1).*zenith_angle_of_arrival_LOS+ ...
                                                      normrnd(0,sigmas_ZSA_LOS/7,1,length(zenith_angle_of_arrival_LOS));     
 
               zenith_angle_of_arrival_LOS         =  rand_factor_add_el - (rand_factor_add_el(1) - obj.theta_arrival(pp));  
               zenith_angle_of_arrival_LOS_per_ray = repmat(zenith_angle_of_arrival_LOS,obj.PerClusterRays,1);        % Assign AOA per each ray in each cluster               
               sign_offset                         = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset)) = -1;
               if obj.generate_uplink
                   cl_wise_alpha_LOS_el = (3/8).*10.^(obj.ZS_A_mu_LOS).*obj.offset_vec.*sign_offset;  %%check ZS_D_mu_LOS
               else
                   cl_wise_alpha_LOS_el = obj.cluster_ZSA_LOS.*obj.offset_vec.*sign_offset;
               end
               zenith_angle_of_arrival_LOS         = zenith_angle_of_arrival_LOS_per_ray + repmat(cl_wise_alpha_LOS_el,1,length(zenith_angle_of_arrival_LOS)); 
               zenith_angle_of_arrival_LOS         = wrapTo360(zenith_angle_of_arrival_LOS);
               % Check if ZOA angle values are within interval [180°,360°] and perform theta_ZOA = 360° - theta_ZOA
               for i=1:length(zenith_angle_of_arrival_LOS)
                   for j=1:obj.NumClusters_LOS
                       if (zenith_angle_of_arrival_LOS(i,j)>=180) && (zenith_angle_of_arrival_LOS(i,j)<=360)
                           zenith_angle_of_arrival_LOS(i,j) = 360 - zenith_angle_of_arrival_LOS(i,j);
                       end
                   end
               end
           end
           
           % Generate Zenith-of-Arrival (ZOA), NLOS
           function zenith_angle_of_arrival_NLOS = zenith_angle_of_arrival_NLOS(obj,sigmas_NLOS,cl_power_NLOS,pp)
               % Generates NLOS (Non-Line of Sight) zenith angle of arrival in [°]
               % See Step 7 in TR 36.873
               % input:    sigmas_NLOS                      ... correlated LSP (Large Scale Parameter) vector:
               %                                               [SF, 0, DS, ASD, ASA, ZSD, ZSA]
               %           cl_power_NLOS                    ... cluster powers 
               %           pp                               ... a counter used when interfering links channel impulse response is generated
               %
               % output:   zenith_angle_of_arrival_NLOS     ... zenith angle of arrival in NLOS: for each cluster and each ray within the cluster
               %                                               a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                           nCl- number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sigmas_ZSA_NLOS                       = sigmas_NLOS(7);
               idx                                   = (obj.C_constant_elevation(1,:) == obj.NumClusters_NLOS);
               C_const_el                            = obj.C_constant_elevation(2,idx);
               zenith_angle_of_arrival_NLOS          = (-sigmas_ZSA_NLOS.*log(cl_power_NLOS./max(cl_power_NLOS)))./C_const_el;  

               zenith_angle_of_arrival_NLOS          = (2.*round(rand(1,length(zenith_angle_of_arrival_NLOS)))-1).*zenith_angle_of_arrival_NLOS+...
                                                         normrnd(0,sigmas_ZSA_NLOS/7,1,length(zenith_angle_of_arrival_NLOS)) + obj.theta_arrival(pp);
               zenith_angle_of_arrival_NLOS_per_ray   = repmat(zenith_angle_of_arrival_NLOS,obj.PerClusterRays,1);                                                                                                         % and add component or random normal distributiom N(0,sigma_ZSA^2/7^2)
               sign_offset = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))  = -1;
               if obj.generate_uplink
                   cl_wise_alpha_NLOS_el = (3/8).*10.^(obj.ZS_A_mu_NLOS)*obj.offset_vec.*sign_offset;
               else
                   cl_wise_alpha_NLOS_el = obj.cluster_ZSA_NLOS.*obj.offset_vec.*sign_offset;
               end
               zenith_angle_of_arrival_NLOS          = zenith_angle_of_arrival_NLOS_per_ray + repmat(cl_wise_alpha_NLOS_el,1,length(zenith_angle_of_arrival_NLOS)); % Add offset angles per ray in each cluster
               zenith_angle_of_arrival_NLOS          = wrapTo360(zenith_angle_of_arrival_NLOS);
               % Check if ZOA angle values are within interval [180°,360°] and perform theta_ZOA = 360° - theta_ZOA
               for i=1:length(zenith_angle_of_arrival_NLOS)
                   for j=1:obj.NumClusters_NLOS
                       if (zenith_angle_of_arrival_NLOS(i,j)>=180) && (zenith_angle_of_arrival_NLOS(i,j)<=360)
                           zenith_angle_of_arrival_NLOS(i,j) = 360 - zenith_angle_of_arrival_NLOS(i,j);
                       end
                   end
               end
           end
           
           function zenith_angle_of_arrival_OTOI = zenith_angle_of_arrival_OTOI(obj,sigmas_OTOI,cl_power_OTOI,UE_is_indoor,pp)
               % Generates O-to-I (Outdoor-to-Indoor) zenith angle of arrival in [°]
               % See Step 7 in TR 36.873
               % input:    sigmas_OTOI                      ... correlated LSP (Large Scale Parameter) vector:
               %                                               [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %           cl_power_OTOI                    ... cluster powers 
               %           UE_is_indoor                     ... actual user condition indoor/outdoor (boolean type)  
               %           pp                               ... a counter used when interfering links channel impulse response is generated
               %
               % output:   zenith_angle_of_arrival_OTOI     ... zenith angle of arrival in O-to-I scenario: for each cluster and each ray within the cluster
               %                                                a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                             nCl- number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sigmas_ZSA_OTOI                       = sigmas_OTOI(7);
               idx                                   = (obj.C_constant_elevation(1,:) == obj.NumClusters_OTOI);
               C_const_el                            = obj.C_constant_elevation(2,idx);
               zenith_angle_of_arrival_OTOI          = (-sigmas_ZSA_OTOI.*log(cl_power_OTOI./max(cl_power_OTOI)))./C_const_el;  
               
               if UE_is_indoor == true
                   theta_ZOA_OTOI = 90;
               else
                   theta_ZOA_OTOI = obj.theta_arrival(pp);
               end
               
               zenith_angle_of_arrival_OTOI          = (2.*round(rand(1,length(zenith_angle_of_arrival_OTOI)))-1).*zenith_angle_of_arrival_OTOI+...
                                                         normrnd(0,sigmas_ZSA_OTOI/7,1,length(zenith_angle_of_arrival_OTOI)) + theta_ZOA_OTOI;
               zenith_angle_of_arrival_OTOI_per_ray   = repmat(zenith_angle_of_arrival_OTOI,obj.PerClusterRays,1);                                                                                                         % and add component or random normal distributiom N(0,sigma_ZSA^2/7^2)
               sign_offset = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))  = -1;
               if obj.generate_uplink
                   cl_wise_alpha_OTOI_el = (3/8).*10.^(obj.ZS_A_mu_OTOI)*obj.offset_vec.*sign_offset;
               else
                   cl_wise_alpha_OTOI_el = obj.cluster_ZSA_OTOI.*obj.offset_vec.*sign_offset;
               end
               zenith_angle_of_arrival_OTOI          = zenith_angle_of_arrival_OTOI_per_ray + repmat(cl_wise_alpha_OTOI_el,1,length(zenith_angle_of_arrival_OTOI)); % Add offset angles per ray in each cluster
               zenith_angle_of_arrival_OTOI          = wrapTo360(zenith_angle_of_arrival_OTOI);
               % Check if ZOA angle values are within interval [180°,360°] and perform theta_ZOA = 360° - theta_ZOA
               for i=1:length(zenith_angle_of_arrival_OTOI)
                   for j=1:obj.NumClusters_OTOI
                       if (zenith_angle_of_arrival_OTOI(i,j)>=180) && (zenith_angle_of_arrival_OTOI(i,j)<=360)
                           zenith_angle_of_arrival_OTOI(i,j) = 360 - zenith_angle_of_arrival_OTOI(i,j);
                       end
                   end
               end
           end
           
           
           % Generate Zenith-of-Departure (ZOD), LOS
           function zenith_angle_of_departure_LOS = zenith_angle_of_departure_LOS(obj,sigmas_LOS,ZOD_parameters,cl_power_LOS,pp)
               % Generates LOS (Line of Sight) zenith angle of departure in [°] 
               % See Step 7 in TR 36.873
               % input:    sigmas_LOS                       ... correlated LSP (Large Scale Parameter) vector:
               %                                                [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %           ZOD_parameters                   ... zenith spread offset parameters generated from Table 7.3-7/8 as
               %                                                a vector: [mean,std,offset_value] of the actual zenith spread
               %           cl_power_LOS                     ... cluster powers 
               %           pp                               ... a counter used when interfering links channel impulse response is generated
               %
               % output:   zenith_angle_of_departure_LOS    ... zenith angle of departure in LOS: for each cluster and each ray within the cluster
               %                                                a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                             nCl- number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sigma_KF_LOS                            = obj.KF_mu_LOS; 
               sigmas_ZSD_LOS                          = sigmas_LOS(6);
               ZS_D_mu_LOS                             = ZOD_parameters(1);
               ZS_D_mu_offset_LOS                      = ZOD_parameters(3); 
               idx                                     = (obj.C_constant_elevation(1,:) == obj.NumClusters_LOS);
               C_const_el_LOS                          = obj.C_constant_elevation(2,idx); 
               C_LOS_el                                = C_const_el_LOS.*(1.3086-0.0339.*sigma_KF_LOS -0.0077.*sigma_KF_LOS.^2 + 0.0002.*sigma_KF_LOS.^3);           
               zenith_angle_of_departure_LOS           = (-sigmas_ZSD_LOS.*log(cl_power_LOS./max(cl_power_LOS)))./C_LOS_el;                                    
               rand_factor_add_el                      = (2.*round(rand(1,length(zenith_angle_of_departure_LOS)))-1).*zenith_angle_of_departure_LOS+ ...
                                                         normrnd(0,sigmas_ZSD_LOS/7,1,length(zenith_angle_of_departure_LOS));                                 
                                                                                                                             
               zenith_angle_of_departure_LOS           = rand_factor_add_el - rand_factor_add_el(1) + obj.theta_departure(pp) + ZS_D_mu_offset_LOS; %deifne_angle changes for indoor 90°  %LOS case equation 7.11 
               zenith_angle_of_departure_LOS_per_ray   = repmat(zenith_angle_of_departure_LOS,obj.PerClusterRays,1);                
               sign_offset                             = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))    = -1;
               if obj.generate_uplink
                   cl_wise_alpha_LOS_el = obj.cluster_ZSA_LOS.*obj.offset_vec.*sign_offset;
               else
                   cl_wise_alpha_LOS_el = (3/8).*10.^(ZS_D_mu_LOS).*obj.offset_vec.*sign_offset;  %%check ZS_D_mu_LOS
               end
               
               zenith_angle_of_departure_LOS           = zenith_angle_of_departure_LOS_per_ray + repmat(cl_wise_alpha_LOS_el,1,length(zenith_angle_of_departure_LOS));   
               zenith_angle_of_departure_LOS           = wrapTo360(zenith_angle_of_departure_LOS);
               % Check if ZOA angle values are within interval [180°,360°] and perform theta_ZOA = 360° - theta_ZOA
               for i=1:length(zenith_angle_of_departure_LOS)
                   for j=1:obj.NumClusters_LOS
                       if (zenith_angle_of_departure_LOS(i,j)>=180) && (zenith_angle_of_departure_LOS(i,j)<=360)
                           zenith_angle_of_departure_LOS(i,j) = 360 - zenith_angle_of_departure_LOS(i,j);
                       end
                   end
               end 
           end
 
           
           % Generate Zenith-of-Departure (ZOD), NLOS
           function zenith_angle_of_departure_NLOS = zenith_angle_of_departure_NLOS(obj,sigmas_NLOS,ZOD_parameters,cl_power_NLOS,pp) 
               % Generates NLOS (Non-Line of Sight) zenith angle of departure in [°] 
               % See Step 7 in TR 36.873
               % input:    sigmas_NLOS                      ... correlated LSP (Large Scale Parameter) vector:
               %                                                [SF, 0, DS, ASD, ASA, ZSD, ZSA]
               %           ZOD_parameters                   ... zenith spread offset parameters generated from Table 7.3-7/8 as
               %                                                a vector: [mean,std,offset_value] of the actual zenith spread
               %           cl_power_NLOS                    ... cluster powers 
               %           pp                               ... a counter used when interfering links channel impulse response is generated
               %
               % output:   zenith_angle_of_departure_NLOS    ... zenith angle of departure in NLOS: for each cluster and each ray within the cluster
               %                                                 a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                             nCl - number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sigmas_ZSD_NLOS                        = sigmas_NLOS(6); 
               ZS_D_mu_NLOS                           = ZOD_parameters(1);
               ZOD_mu_offset_NLOS                     = ZOD_parameters(3);
               idx                                    = (obj.C_constant_elevation(1,:) == obj.NumClusters_NLOS);
               C_const_el                             = obj.C_constant_elevation(2,idx);
               zenith_angle_of_departure_NLOS         = (-sigmas_ZSD_NLOS.*log(cl_power_NLOS./max(cl_power_NLOS)))./C_const_el;                               
               zenith_angle_of_departure_NLOS         = (2.*round(rand(1,length(zenith_angle_of_departure_NLOS)))-1).*zenith_angle_of_departure_NLOS...
                                                          + normrnd(0,sigmas_ZSD_NLOS/7,1,length(zenith_angle_of_departure_NLOS))...
                                                           + obj.theta_departure(pp) + ZOD_mu_offset_NLOS;                         % Assign a positive or negative sign to the angles eq.7.18
               zenith_angle_of_departure_NLOS_per_ray = repmat(zenith_angle_of_departure_NLOS,obj.PerClusterRays,1);                                                                                                                                  
               sign_offset                            = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))   = -1;
               if obj.generate_uplink
                   cl_wise_alpha_NLOS_el = obj.cluster_ZSA_NLOS.*obj.offset_vec.*sign_offset;
               else
                   cl_wise_alpha_NLOS_el = (3/8).*10.^(ZS_D_mu_NLOS)*obj.offset_vec.*sign_offset;
               end
               zenith_angle_of_departure_NLOS         = zenith_angle_of_departure_NLOS_per_ray + repmat(cl_wise_alpha_NLOS_el,1,length(zenith_angle_of_departure_NLOS));       
               zenith_angle_of_departure_NLOS         = wrapTo360(zenith_angle_of_departure_NLOS);
               % Check if ZOD angle values are within interval [180°,360°] and perform theta_ZOD = 360° - theta_ZOD
               for i=1:length(zenith_angle_of_departure_NLOS)
                   for j=1:obj.NumClusters_NLOS
                       if (zenith_angle_of_departure_NLOS(i,j)>=180) && (zenith_angle_of_departure_NLOS(i,j)<=360)
                           zenith_angle_of_departure_NLOS(i,j) = 360 - zenith_angle_of_departure_NLOS(i,j);
                       end
                   end
               end
           end
                
            function zenith_angle_of_departure_OTOI = zenith_angle_of_departure_OTOI(obj,sigmas_OTOI,ZOD_parameters,cl_power_OTOI,pp) 
               % Generates O-to-I (Outdoor-to-Indoor) zenith angle of departure in [°] 
               % See Step 7 in TR 36.873
               % input:    sigmas_OTOI                      ... correlated LSP (Large Scale Parameter) vector:
               %                                                [SF, 0, DS, ASD, ASA, ZSD, ZSA]
               %           ZOD_parameters                   ... zenith spread offset parameters generated from Table 7.3-7/8 as
               %                                                a vector: [mean,std,offset_value] of the actual zenith spread
               %           cl_power_OTOI                    ... cluster powers 
               %           pp                               ... a counter used when interfering links channel impulse response is generated
               %
               % output:   zenith_angle_of_departure_OTOI   ... zenith angle of departure in O-to-I: for each cluster and each ray within the cluster
               %                                                 a matrix of size [nR x nCl]: nR - number of rays per cluster
               %                                                                             nCl - number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sigmas_ZSD_OTOI                        = sigmas_OTOI(6); 
               ZS_D_mu_OTOI                           = ZOD_parameters(1);
               ZOD_mu_offset_OTOI                     = ZOD_parameters(3);
               idx                                    = (obj.C_constant_elevation(1,:) == obj.NumClusters_OTOI);
               C_const_el                             = obj.C_constant_elevation(2,idx);
               zenith_angle_of_departure_OTOI         = (-sigmas_ZSD_OTOI.*log(cl_power_OTOI./max(cl_power_OTOI)))./C_const_el;                               
               zenith_angle_of_departure_OTOI         = (2.*round(rand(1,length(zenith_angle_of_departure_OTOI)))-1).*zenith_angle_of_departure_OTOI...
                                                          + normrnd(0,sigmas_ZSD_OTOI/7,1,length(zenith_angle_of_departure_OTOI))...
                                                           + obj.theta_departure(pp) + ZOD_mu_offset_OTOI;                         % Assign a positive or negative sign to the angles eq.7.18
               zenith_angle_of_departure_OTOI_per_ray = repmat(zenith_angle_of_departure_OTOI,obj.PerClusterRays,1);                                                                                                                                  
               sign_offset                            = ones(obj.PerClusterRays,1);
               sign_offset(2:2:length(sign_offset))   = -1;
               if obj.generate_uplink
                   cl_wise_alpha_OTOI_el = obj.cluster_ZSA_OTOI.*obj.offset_vec.*sign_offset;
               else
                   cl_wise_alpha_OTOI_el = (3/8).*10.^(ZS_D_mu_OTOI)*obj.offset_vec.*sign_offset;
               end
               zenith_angle_of_departure_OTOI         = zenith_angle_of_departure_OTOI_per_ray + repmat(cl_wise_alpha_OTOI_el,1,length(zenith_angle_of_departure_OTOI));       
               zenith_angle_of_departure_OTOI         = wrapTo360(zenith_angle_of_departure_OTOI);
               % Check if ZOD angle values are within interval [180°,360°] and perform theta_ZOD = 360° - theta_ZOD
               for i=1:length(zenith_angle_of_departure_OTOI)
                   for j=1:obj.NumClusters_OTOI
                       if (zenith_angle_of_departure_OTOI(i,j)>=180) && (zenith_angle_of_departure_OTOI(i,j)<=360)
                           zenith_angle_of_departure_OTOI(i,j) = 360 - zenith_angle_of_departure_OTOI(i,j);
                       end
                   end
               end
           end
           %% Step 8: Coupling rays within a cluster for both azimuth and
           % elevation
            function [azimuth_angle_of_departure_LOS,...
                     azimuth_angle_of_arrival_LOS,...
                     zenith_angle_of_departure_LOS,...
                     zenith_angle_of_arrival_LOS]      = coupling_rays_LOS(obj,sigmas_LOS,ZOD_parameters,cl_power_LOS,pp)
               % Couplning of rays within a cluster for both azimuth- and zenith angles of arrival- and departure in LOS (Line of Sight)
               % Couples randlomly:     - angles of departure with angles of arrival in azimuth within each cluster 
               %                        - angles of departure with angles of arrival in zenith within each cluster 
               % See Step 8 in TR 36.873
               % input:    sigmas_LOS                       ... correlated LSP (Large Scale Parameter) vector:
               %                                                [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %           ZOD_parameters                   ... zenith spread offset parameters generated from Table 7.3-7/8 as
               %                                                a vector: [mean,std,offset_value] of the actual zenith spread
               %           cl_power_LOS                     ... cluster powers 
               %           pp                               ... a counter used when interfering links channel impulse response is generated
               %
               % output:   azimuth_angle_of_departure_LOS   ... azimuth angle of departure 
               %           azimuth_angle_of_arrival_LOS     ... azimuth angle of arrival
               %           zenith_angle_of_departure_LOS    ... zenith angle of departure
               %           zenith_angle_of_arrival_LOS      ... zenith angle of arrival
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016  
                 
               [dummy, h]                              = sort(rand(obj.PerClusterRays,obj.NumClusters_LOS),1);                                         
               inds                                    = h + repmat([1:obj.PerClusterRays:obj.NumClusters_LOS*obj.PerClusterRays],obj.PerClusterRays,1)-1;
               % Couple randomly AOD to AOA angles within a cluster
               azimuth_angle_of_departure_LOS          = obj.azimuth_angle_of_departure_LOS(sigmas_LOS,cl_power_LOS,pp);
               azimuth_angle_of_departure_LOS          = azimuth_angle_of_departure_LOS(inds);
               azimuth_angle_of_arrival_LOS            = obj.azimuth_angle_of_arrival_LOS(sigmas_LOS,cl_power_LOS,pp);
               azimuth_angle_of_arrival_LOS            = azimuth_angle_of_arrival_LOS(inds);
                % Couple randomly ZOD to ZOA angles within a cluster
               zenith_angle_of_departure_LOS           = obj.zenith_angle_of_departure_LOS(sigmas_LOS,ZOD_parameters,cl_power_LOS,pp);
               zenith_angle_of_departure_LOS           = zenith_angle_of_departure_LOS(inds);
               zenith_angle_of_arrival_LOS             = obj.zenith_angle_of_arrival_LOS(sigmas_LOS,cl_power_LOS,pp);
               zenith_angle_of_arrival_LOS             = zenith_angle_of_arrival_LOS(inds);
            end
           
     
           function [azimuth_angle_of_departure_NLOS,...
                     azimuth_angle_of_arrival_NLOS,...
                     zenith_angle_of_departure_NLOS,...
                     zenith_angle_of_arrival_NLOS]     = coupling_rays_NLOS(obj,sigmas_NLOS,ZOD_parameters,cl_power_NLOS,pp)
               % Couplning of rays within a cluster for both azimuth- and zenith angles of arrival- and departure in NLOS (Non-Line of Sight)
               % Couples randlomly:     - angles of departure with angles of arrival in azimuth within each cluster 
               %                        - angles of departure with angles of arrival in zenith within each cluster 
               % See Step 8 in TR 36.873
               % input:    sigmas_NLOS                       ... correlated LSP (Large Scale Parameter) vector:
               %                                                [SF, 0, DS, ASD, ASA, ZSD, ZSA]
               %           ZOD_parameters                    ... zenith spread offset parameters generated from Table 7.3-7/8 as
               %                                                a vector: [mean,std,offset_value] of the actual zenith spread
               %           cl_power_NLOS                     ... cluster powers 
               %           pp                                ... a counter used when interfering links channel impulse response is generated
               % 
               % output:   azimuth_angle_of_departure_NLOS   ... azimuth angle of departure 
               %           azimuth_angle_of_arrival_NLOS     ... azimuth angle of arrival
               %           zenith_angle_of_departure_NLOS    ... zenith angle of departure
               %           zenith_angle_of_arrival_NLOS      ... zenith angle of arrival
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016  

               [k, h]                                  = sort(rand(obj.PerClusterRays,obj.NumClusters_NLOS),1);    % Create Nr_of_clusters random permutations of integers [1:NR-of_rays]
               inds                                    = h + repmat([1:obj.PerClusterRays:obj.NumClusters_NLOS*obj.PerClusterRays],obj.PerClusterRays,1)-1;
               % Couple randomly AOD to AOA angles within a cluster
               azimuth_angle_of_departure_NLOS        = obj.azimuth_angle_of_departure_NLOS(sigmas_NLOS,cl_power_NLOS,pp);
               azimuth_angle_of_departure_NLOS        = azimuth_angle_of_departure_NLOS(inds);
               azimuth_angle_of_arrival_NLOS          = obj.azimuth_angle_of_arrival_NLOS(sigmas_NLOS,cl_power_NLOS,pp);
               azimuth_angle_of_arrival_NLOS          = azimuth_angle_of_arrival_NLOS(inds);
                % Couple randomly ZOD to ZOA angles within a cluster
               zenith_angle_of_departure_NLOS         = obj.zenith_angle_of_departure_NLOS(sigmas_NLOS,ZOD_parameters,cl_power_NLOS,pp);
               zenith_angle_of_departure_NLOS         = zenith_angle_of_departure_NLOS(inds);
               zenith_angle_of_arrival_NLOS           = obj.zenith_angle_of_arrival_NLOS(sigmas_NLOS,cl_power_NLOS,pp);
               zenith_angle_of_arrival_NLOS           = zenith_angle_of_arrival_NLOS(inds);
           end
          
           function [azimuth_angle_of_departure_OTOI,...
                     azimuth_angle_of_arrival_OTOI,...
                     zenith_angle_of_departure_OTOI,...
                     zenith_angle_of_arrival_OTOI]     = coupling_rays_OTOI(obj,sigmas_OTOI,ZOD_parameters,cl_power_OTOI,UE_is_indoor,pp)
               % Couplning of rays within a cluster for both azimuth- and zenith angles of arrival- and departure in O-to-I (Outdoor-to-Indoor) scenario
               % Couples randlomly:     - angles of departure with angles of arrival in azimuth within each cluster 
               %                        - angles of departure with angles of arrival in zenith within each cluster 
               % See Step 8 in TR 36.873
               % input:    sigmas_OTOI                       ... correlated LSP (Large Scale Parameter) vector:
               %                                                [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
               %           ZOD_parameters                    ... zenith spread offset parameters generated from Table 7.3-7/8 as
               %                                                a vector: [mean,std,offset_value] of the actual zenith spread
               %           cl_power_OTOI                     ... cluster powers 
               %           pp                                ... a counter used when interfering links channel impulse response is generated
               %
               % output:   azimuth_angle_of_departure_OTOI   ... azimuth angle of departure 
               %           azimuth_angle_of_arrival_OTOI     ... azimuth angle of arrival
               %           zenith_angle_of_departure_OTOI    ... zenith angle of departure
               %           zenith_angle_of_arrival_OTOI      ... zenith angle of arrival
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016  
                 
               [k, h]                                  = sort(rand(obj.PerClusterRays,obj.NumClusters_OTOI),1);    % Create Nr_of_clusters random permutations of integers [1:NR-of_rays]
               inds                                    = h + repmat([1:obj.PerClusterRays:obj.NumClusters_OTOI*obj.PerClusterRays],obj.PerClusterRays,1)-1;
               % Couple randomly AOD to AOA angles within a cluster
               azimuth_angle_of_departure_OTOI        = obj.azimuth_angle_of_departure_OTOI(sigmas_OTOI,cl_power_OTOI,pp);
               azimuth_angle_of_departure_OTOI        = azimuth_angle_of_departure_OTOI(inds);
               azimuth_angle_of_arrival_OTOI          = obj.azimuth_angle_of_arrival_OTOI(sigmas_OTOI,cl_power_OTOI,pp);
               azimuth_angle_of_arrival_OTOI          = azimuth_angle_of_arrival_OTOI(inds);
                % Couple randomly ZOD to ZOA angles within a cluster
               zenith_angle_of_departure_OTOI         = obj.zenith_angle_of_departure_OTOI(sigmas_OTOI,ZOD_parameters,cl_power_OTOI,pp);
               zenith_angle_of_departure_OTOI         = zenith_angle_of_departure_OTOI(inds);
               zenith_angle_of_arrival_OTOI           = obj.zenith_angle_of_arrival_OTOI(sigmas_OTOI,cl_power_OTOI,UE_is_indoor,pp);
               zenith_angle_of_arrival_OTOI           = zenith_angle_of_arrival_OTOI(inds);
           end
           
           %% Step 9: Generate Cross Polarization Power Rations (XPR)
           % LOS - XPR
           function xpr_LOS = generate_XPR_LOS(obj)
               % Generate Cross Polarization Power Rations for LOS 
               % See Step 9 in TR 36.873
               %
               % output:    xpr_LOS ... XPR value generated by log-Normal distribution
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016  
               X_LOS    = normrnd(obj.xpr_mu_LOS,obj.xpr_sigma_LOS,obj.PerClusterRays,obj.NumClusters_LOS);
               xpr_LOS  = 10.^(X_LOS./10);  % Log-Normal distribution of XPR
           end
           
           % NLOS - XPR
           function xpr_NLOS = generate_XPR_NLOS(obj)
               % Generate Cross Polarization Power Rations for NLOS 
               % See Step 9 in TR 36.873
               %
               % output:    xpr_NLOS ... XPR value generated by log-Normal distribution
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016  
               X_NLOS   = normrnd(obj.xpr_mu_NLOS,obj.xpr_sigma_NLOS,obj.PerClusterRays,obj.NumClusters_NLOS);
               xpr_NLOS = 10.^(X_NLOS./10); 
           end
           
           function xpr_OTOI = generate_XPR_OTOI(obj)
               % Generate Cross Polarization Power Rations for O-to-I 
               % See Step 9 in TR 36.873
               %
               % output:    xpr_OTOI ... XPR value generated by log-Normal distribution
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016  
               X_OTOI   = normrnd(obj.xpr_mu_OTOI,obj.xpr_sigma_OTOI,obj.PerClusterRays,obj.NumClusters_OTOI);
               xpr_OTOI = 10.^(X_OTOI./10); 
           end

           %% Draw initial phases for each ray and each cluster for four
           % different polarization combinations
           function initial_phases_LOS = generate_initial_phases_LOS(obj)
               % Generates initial phases for each ray of each cluster and for four polarization combinations
               % The distribution of intial phases is uniform in [-pi,pi] 
               % See Step 10 in TR 36.873
               %
               % output        initial_phases_LOS ... a matrix of size [nR x nCl x 4] 
               %                                       nR - number of rays per cluster
               %                                       nCl - number of clusters
               %                                       4 - number of polarization combinations 
               %
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016  
               initial_phases_LOS = 360*rand(obj.PerClusterRays,obj.NumClusters_LOS,4);               
               initial_phases_LOS = wrapTo180(initial_phases_LOS);    % Wrap angles from interval [0°,360°] to [-180°,180°]
           end
           
           function initial_phases_NLOS  = generate_initial_phases_NLOS(obj)
               % Generates initial phases for each ray of each cluster and for four polarization combinations
               % The distribution of intial phases is uniform in [-pi,pi] 
               % See Step 10 in TR 36.873
               %
               % output        initial_phases_NLOS ... a matrix of size [nR x nCl x 4] 
               %                                       nR - number of rays per cluster
               %                                       nCl - number of clusters
               %                                       4 - number of polarization combinations 
               %
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               initial_phases_NLOS = 360*rand(obj.PerClusterRays,obj.NumClusters_NLOS,4);            
               initial_phases_NLOS = wrapTo180(initial_phases_NLOS);                                    
           end
           
            function initial_phases_OTOI  = generate_initial_phases_OTOI(obj)
               % Generates initial phases for each ray of each cluster and for four polarization combinations
               % The distribution of intial phases is uniform in [-pi,pi] 
               % See Step 10 in TR 36.873
               %
               % output        initial_phases_OTOI ... a matrix of size [nR x nCl x 4] 
               %                                       nR - number of rays per cluster
               %                                       nCl - number of clusters
               %                                       4 - number of polarization combinations 
               %
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               initial_phases_OTOI = 360*rand(obj.PerClusterRays,obj.NumClusters_OTOI,4);            
               initial_phases_OTOI = wrapTo180(initial_phases_OTOI);                                    
           end
           
           % In LOS case draw initial phases just for two polarization
           % combinations
           function initial_phases_LOS_direct_ray = generate_initial_phases_LOS_direct_ray(obj)
               % Generates initial phases for the direct LOS path and for two polarization combinations
               % The distribution of intial phases is uniform in [-pi,pi] 
               % See Step 10 in TR 36.873
               %
               % output        initial_phases_LOS_direct_ray ... two output phases [theta', phi']: 
               %                                                 theta' - for the polarization in zenith
               %                                                 phi' - for the polarization in azimuth
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               initial_phases_LOS_direct_ray = 360*rand(1,2);                                        
               initial_phases_LOS_direct_ray = wrapTo180(initial_phases_LOS_direct_ray);  % Wrap angles from interval [0°,360°] to [-180°,180°]
           end
                    
           %% Generate spherical unit vectors with azimuth/elevation arrival
           % angle AOA,ZOA for receiver and AOD,ZOD for transmiter
           % -->RX LOS
           
           function RX_spherical_unit_vector_LOS = receiver_spherical_unit_vector_LOS(obj,pp)
               % Generates the spherical unit vector r at the receiver from
               % angle of arrival in azimuth and zenith
               % See eq. 7.22 in TR 36.873
               %
               % input:     zenith_angle_of_arrival_LOS
               %            azimuth_angle_of_arrival_LOS
               %
               % output:    RX_spherical_unit_vector_LOS     % a matrix of three dimensions:
               %                                               - 1st dimension is x,y,z coordinates
               %                                               - 2nd dimension is number of rays
               %                                               - 3rd dimension is number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               
               sin_cos_temp_LOS = shiftdim(sind(obj.zenith_angle_of_arrival_LOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp)),-1);
               sin_sin_temp_LOS = shiftdim(sind(obj.zenith_angle_of_arrival_LOS_once(:,:,pp)).*sind(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp)),-1);
               cos_temp_LOS = shiftdim(cosd(obj.zenith_angle_of_arrival_LOS_once(:,:,pp)),-1);
               RX_spherical_unit_vector_LOS = [sin_cos_temp_LOS; sin_sin_temp_LOS; cos_temp_LOS];  % 1st dimension is x,y,z coordinates, 2nd dimension is rays, 3rd dimension is clusters
           end
           
           function RX_spherical_unit_vector_NLOS = receiver_spherical_unit_vector_NLOS(obj,pp)
               % Generates the spherical unit vector r at the receiver from
               % angle of arrival in azimuth and zenith
               % See eq. 7.22 in TR 36.873
               %
               % input:     zenith_angle_of_arrival_NLOS
               %            azimuth_angle_of_arrival_NLOS
               %
               % output:    RX_spherical_unit_vector_NLOS     % a matrix of three dimensions:
               %                                               - 1st dimension is x,y,z coordinates
               %                                               - 2nd dimension is number of rays
               %                                               - 3rd dimension is number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sin_cos_temp_NLOS = shiftdim(sind(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp)),-1);
               sin_sin_temp_NLOS = shiftdim(sind(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp)).*sind(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp)),-1);
               cos_temp_NLOS = shiftdim(cosd(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp)),-1);
               RX_spherical_unit_vector_NLOS = [sin_cos_temp_NLOS;sin_sin_temp_NLOS; cos_temp_NLOS];
           end
           
           function RX_spherical_unit_vector_OTOI = receiver_spherical_unit_vector_OTOI(obj,pp)
               % Generates the spherical unit vector r at the receiver from
               % angle of arrival in azimuth and zenith
               % See eq. 7.22 in TR 36.873
               %
               % input:     zenith_angle_of_arrival_OTOI
               %            azimuth_angle_of_arrival_OTOI
               %
               % output:    RX_spherical_unit_vector_OTOI     % a matrix of three dimensions:
               %                                               - 1st dimension is x,y,z coordinates
               %                                               - 2nd dimension is number of rays
               %                                               - 3rd dimension is number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sin_cos_temp_OTOI = shiftdim(sind(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp)).*cosd(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp)),-1);
               sin_sin_temp_OTOI = shiftdim(sind(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp)).*sind(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp)),-1);
               cos_temp_OTOI = shiftdim(cosd(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp)),-1);
               RX_spherical_unit_vector_OTOI = [sin_cos_temp_OTOI;sin_sin_temp_OTOI; cos_temp_OTOI];
           end
           
           % TX LOS
           function TX_spherical_unit_vector_LOS = transmitter_spherical_unit_vector_LOS(obj,pp)
               % Generates the spherical unit vector r at the transmitter from
               % angle of departure in azimuth and zenith
               % See eq. 7.23 in TR 36.873
               %
               % input:     zenith_angle_of_departure_LOS
               %            azimuth_angle_of_departure_LOS
               %
               % output:    TX_spherical_unit_vector_LOS     % a matrix of three dimensions:
               %                                               - 1st dimension is x,y,z coordinates
               %                                               - 2nd dimension is number of rays
               %                                               - 3rd dimension is number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sin_cos_temp_LOS = shiftdim(sind(obj.zenith_angle_of_departure_LOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_departure_LOS_once(:,:,pp)),-1);
               sin_sin_temp_LOS = shiftdim(sind(obj.zenith_angle_of_departure_LOS_once(:,:,pp)).*sind(obj.azimuth_angle_of_departure_LOS_once(:,:,pp)),-1);
               cos_temp_LOS = shiftdim(cosd(obj.zenith_angle_of_departure_LOS_once(:,:,pp)),-1);
               TX_spherical_unit_vector_LOS = [sin_cos_temp_LOS;sin_sin_temp_LOS;cos_temp_LOS];
           end
           
           function TX_spherical_unit_vector_NLOS = transmitter_spherical_unit_vector_NLOS(obj,pp)
               % Generates the spherical unit vector r at the transmitter from
               % angle of departure in azimuth and zenith
               % See eq. 7.23 in TR 36.873
               %
               % input:     zenith_angle_of_departure_NLOS
               %            azimuth_angle_of_departure_NLOS
               %
               % output:    TX_spherical_unit_vector_NLOS     % a matrix of three dimensions:
               %                                               - 1st dimension is x,y,z coordinates
               %                                               - 2nd dimension is number of rays
               %                                               - 3rd dimension is number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sin_cos_temp_NLOS = shiftdim(sind(obj.zenith_angle_of_departure_NLOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp)),-1);
               sin_sin_temp_NLOS = shiftdim(sind(obj.zenith_angle_of_departure_NLOS_once(:,:,pp)).*sind(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp)),-1);
               cos_temp_NLOS = shiftdim(cosd(obj.zenith_angle_of_departure_NLOS_once(:,:,pp)),-1);
               TX_spherical_unit_vector_NLOS = [sin_cos_temp_NLOS;sin_sin_temp_NLOS;cos_temp_NLOS];
           end
           
            function TX_spherical_unit_vector_OTOI = transmitter_spherical_unit_vector_OTOI(obj,pp)
               % Generates the spherical unit vector r at the transmitter from
               % angle of departure in azimuth and zenith
               % See eq. 7.23 in TR 36.873
               %
               % input:     zenith_angle_of_departure_OTOI
               %            azimuth_angle_of_departure_OTOI
               %
               % output:    TX_spherical_unit_vector_OTOI     % a matrix of three dimensions:
               %                                               - 1st dimension is x,y,z coordinates
               %                                               - 2nd dimension is number of rays
               %                                               - 3rd dimension is number of clusters
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sin_cos_temp_OTOI = shiftdim(sind(obj.zenith_angle_of_departure_OTOI_once(:,:,pp)).*cosd(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp)),-1);
               sin_sin_temp_OTOI = shiftdim(sind(obj.zenith_angle_of_departure_OTOI_once(:,:,pp)).*sind(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp)),-1);
               cos_temp_OTOI = shiftdim(cosd(obj.zenith_angle_of_departure_OTOI_once(:,:,pp)),-1);
               TX_spherical_unit_vector_OTOI = [sin_cos_temp_OTOI;sin_sin_temp_OTOI;cos_temp_OTOI];
           end
           
           % Generate spherical unit vectors with azimuth/elevation for
           % single LOS ray using the direct (LOS) angles between eNodeB and UE
           function RX_spherical_unit_vector_LOS_direct_ray = receiver_spherical_unit_vector_LOS_direct_ray(obj,pp)
               % Generates the spherical unit vector r at the receiver from
               % angle of arrival in azimuth and zenith for the direct LOS path
               % See eq. 7.22 in TR 36.873
               %
               % input:     theta_arrival                           ... angle of arrival in zenith of the direct path
               %            phi_arrival                             ... angle of arrival in azimuth of the direct path
               % 
               % output:    RX_spherical_unit_vector_LOS_direct_ray ... a vector with dimensions  3x1 of the corresponding [x; y; z] coordinate
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sin_cos_temp = sind(obj.theta_arrival(pp)).*cosd(obj.phi_arrival(pp));
               sin_sin_temp = sind(obj.theta_arrival(pp)).*sind(obj.phi_arrival(pp));
               cos_temp = cosd(obj.theta_arrival(pp));
               RX_spherical_unit_vector_LOS_direct_ray = [sin_cos_temp;sin_sin_temp;cos_temp];
           end
           
           function TX_spherical_unit_vector_LOS_direct_ray = transmitter_spherical_unit_vector_LOS_direct_ray(obj,pp)
               % Generates the spherical unit vector r at the transmitter from
               % angle of departure in azimuth and zenith for the direct LOS path
               % See eq. 7.23 in TR 36.873
               %
               % input:     theta_departure                           ... angle of departure in zenith of the direct path
               %            phi_departure                             ... angle of departure in azimuth of the direct path
               % 
               % output:    TX_spherical_unit_vector_LOS_direct_ray   ... a vector with dimensions  3x1 of the corresponding [x; y; z] coordinate
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               sin_cos_temp = sind(obj.theta_departure(pp)).*cosd(obj.phi_departure(pp));
               sin_sin_temp = sind(obj.theta_departure(pp)).*sind(obj.phi_departure(pp));
               cos_temp = cosd(obj.theta_departure(pp));
               TX_spherical_unit_vector_LOS_direct_ray = [sin_cos_temp;sin_sin_temp;cos_temp];
           end
           %% Generate UE location and eNodeB location
           function UE_location = UE_position(obj,rx_height,UE_pos,relative_UE_antenna_position) 
               % Generates the actual antenna element position at the receiver in cartesian coordinates [x,y,z]
               % Includes the user position and the relative antenna element position
               % 
               % input:      rx_height                      ... user height in [m]
               %             UE_pos                         ... user position in [x,y]
               %             relative_UE_antenna_position   ... antenna element position in [x,y,z] relative to the antenna array at receiver
               %
               % output:     UE_location                    ... actual antenna element position relative to the network geometry
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               UE_location = [UE_pos(1),UE_pos(2),rx_height]+relative_UE_antenna_position;
           end
           % eNodeB location
           function eNodeB_location = eNodeB_position(obj,LTE_config,attached_site_pos,relative_antenna_position)
               % Generates the actual antenna element position at the transmitter in cartesian coordinates [x,y,z]
               % Includes the base-station position and the relative antenna element position
               % 
               % input:      tx_height                      ... base-station height in [m]
               %             attached_site_pos              ... base-station position in [x,y]
               %             relative_UE_antenna_position   ... antenna element position in [x,y,z] relative to the antenna array at transmitter
               %
               % output:     eNodeB_location                    ... actual antenna element position relative to the network geometry
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               tx_height = LTE_config.tx_height;
               eNodeB_location = [attached_site_pos(1),attached_site_pos(2),tx_height]+relative_antenna_position;
           end
           
 
           % Generate Doppler frequency component Eq(7.24)
           function v_doppler_LOS = calculate_doppler_speed_LOS(obj,LTE_config, direction_of_movement, pp)
               % Generates the Doppler frequency component v
               % See eq. 7.24 in TR 36.873
               %
               % input:           frequency
               %                  UE_speed
               %                  direction_of_movement 
               %                  
               % output:          v_doppler_LOS
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
                frequency           = LTE_config.frequency; 
                UE_speed            = LTE_config.UE_speed.*[cosd(direction_of_movement);sind(direction_of_movement);0];
                receiver_spherical_unit_vector_LOS = obj.receiver_spherical_unit_vector_LOS(pp);
                size_vec=size(receiver_spherical_unit_vector_LOS);
                unit_vector_temp = reshape(receiver_spherical_unit_vector_LOS,3,size_vec(2)*size_vec(3));  
                v_doppler_LOS_temp = unit_vector_temp.'*UE_speed./(299792458./frequency);               
                v_doppler_LOS       = reshape(v_doppler_LOS_temp,size_vec(2),size_vec(3));                       
           end
           
           function v_doppler_NLOS = calculate_doppler_speed_NLOS(obj,LTE_config, direction_of_movement, pp)
               % Generates the Doppler frequency component v
               % See eq. 7.24 in TR 36.873
               %
               % input:           frequency
               %                  UE_speed
               %                  direction_of_movement 
               %                  
               % output:          v_doppler_NLOS
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               frequency            = LTE_config.frequency;
               UE_speed             = LTE_config.UE_speed.*[cosd(direction_of_movement);sind(direction_of_movement);0];
               receiver_spherical_unit_vector_NLOS = obj.receiver_spherical_unit_vector_NLOS(pp); 
               size_vec = size(receiver_spherical_unit_vector_NLOS);
               unit_vector_temp = reshape(receiver_spherical_unit_vector_NLOS,3,size_vec(2)*size_vec(3));
               v_doppler_NLOS_temp = unit_vector_temp.'*UE_speed./(299792458./frequency); 
               v_doppler_NLOS       = reshape(v_doppler_NLOS_temp,size_vec(2),size_vec(3));  
           end
           
           function v_doppler_OTOI = calculate_doppler_speed_OTOI(obj,LTE_config, direction_of_movement, pp)
               % Generates the Doppler frequency component v
               % See eq. 7.24 in TR 36.873
               %
               % input:           frequency
               %                  UE_speed
               %                  direction_of_movement 
               %                  
               % output:          v_doppler_OTOI
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               frequency            = LTE_config.frequency;
               UE_speed             = LTE_config.UE_speed.*[cosd(direction_of_movement);sind(direction_of_movement);0];
               receiver_spherical_unit_vector_OTOI = obj.receiver_spherical_unit_vector_OTOI(pp); 
               size_vec = size(receiver_spherical_unit_vector_OTOI);
               unit_vector_temp = reshape(receiver_spherical_unit_vector_OTOI,3,size_vec(2)*size_vec(3));
               v_doppler_OTOI_temp = unit_vector_temp.'*UE_speed./(299792458./frequency); 
               v_doppler_OTOI       = reshape(v_doppler_OTOI_temp,size_vec(2),size_vec(3));  
           end
           
            function v_doppler_LOS_direct_ray = calculate_doppler_speed_LOS_direct_ray(obj,LTE_config, direction_of_movement, pp)
               % Generates the Doppler frequency component v for the direct LOS path
               % See eq. 7.24 in TR 36.873
               %
               % input:           frequency
               %                  UE_speed
               %                  direction_of_movement 
               %                  
               % output:          v_doppler_LOS_direct_ray
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
                frequency           = LTE_config.frequency;  
                UE_speed            = LTE_config.UE_speed.*[cosd(direction_of_movement);sind(direction_of_movement);0];
                receiver_spherical_unit_vector_LOS_direct_ray = obj.receiver_spherical_unit_vector_LOS_direct_ray(pp);   
                v_doppler_LOS_direct_ray =  receiver_spherical_unit_vector_LOS_direct_ray.'*UE_speed./(299792458./frequency);                           
           end
           
           %% Generate Tx antenna filed pattern in GCS as a function of
           % previously generated ZOD and AOD angles, transformation from
           % LCS (Local Coordinate System) to GCS (Global Coordinate System)
           % LOS 
            function [theta_AntennaField_tx_global_LOS, phi_AntennaField_tx_global_LOS] = antenna_field_pattern_GCS_LOS(obj,LTE_config,attached_eNodeB,slant_angle,pp) 
               % Translates the antenna element field pattern from Local Coordinate System (LCS) in Global Coordinate System (GCS) 
               % See Section 5.1.3 in TR 36.873
               %
               % input:                        theta_antenna_field_tx_LOS        ... polarized antenna element field pattern in zenith in LCS
               %                               phi_antenna_field_tx_LOS          ... polarized antenna element field pattern in azimuth in LCS
               %                               bearing_angle                     ... pointing direction of antenna relative to x axis (0° in +x direction) used for rotation 
               %                               mechanical_downtilt               ... direction of antenna element relative to the y axis (0° in +y direction) used for rotation 
               %                               mechanical_slant                  ... direction of antenna element relative to the z axis (0° in +z direction) used for rotation 
               %                               slant_angle                       ... the angle of slanted antenna element relative to the z axis (0° in +z direction) 
               %                                                                     depends on the plarization mode (defined in the config)
               %
               % output:                       theta_AntennaField_tx_global_LOS  ... polarized antenna element field pattern in zenith in GCS
               %                               phi_AntennaField_tx_global_LOS     ...polarized antenna element field pattern in azimuth in GCS
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               bearing_angle = obj.BS_boresight;% sector boresight angle with respect to x-axis 
               mechanical_downtilt = LTE_config.mechanical_downtilt; % mechanical downtilt between 0° and 180°
               mechanical_slant = LTE_config.mechanical_slant; %mechanical slant angle assumed to be zero
               
               % Psi angle as a functio of mechanical orienation Eq. 5.13
                psi_angle_LOS = angle(sind(mechanical_slant).*cosd(obj.zenith_angle_of_departure_LOS_once(:,:,pp)).*sind(obj.azimuth_angle_of_departure_LOS_once(:,:,pp) - bearing_angle)+...
                                     cosd(mechanical_slant).*(cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_departure_LOS_once(:,:,pp)) - sind(mechanical_downtilt).*...
                                     cosd(obj.zenith_angle_of_departure_LOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_departure_LOS_once(:,:,pp) - bearing_angle))+...
                                         1i.*(sind(mechanical_slant).*cosd(obj.azimuth_angle_of_departure_LOS_once(:,:,pp) - bearing_angle)+...
                                         sind(mechanical_downtilt).*cosd(mechanical_slant).*sind(obj.azimuth_angle_of_departure_LOS_once(:,:,pp) - bearing_angle)));
               psi_angle_LOS = psi_angle_LOS.*(180/pi);
               
               % Convert azimuth and elevation angles from GCS to LCS (as seen from antenna sector)
               theta_global_LOS = acosd(cosd(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.zenith_angle_of_departure_LOS_once(:,:,pp))+...
                                     (sind(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.azimuth_angle_of_departure_LOS_once(:,:,pp) - bearing_angle)-...
                                       sind(mechanical_slant).*sin(obj.azimuth_angle_of_departure_LOS_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_departure_LOS_once(:,:,pp)));
               
               phi_global_LOS = angle((cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_departure_LOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_departure_LOS_once(:,:,pp) - bearing_angle)-...
                                              sind(mechanical_downtilt).*cosd(obj.zenith_angle_of_departure_LOS_once(:,:,pp)))+...
                                              1i.*((cosd(mechanical_downtilt)).*sind(mechanical_slant).*cosd(obj.zenith_angle_of_departure_LOS_once(:,:,pp))+...
                                                 (sind(mechanical_downtilt).*sind(mechanical_slant).*cosd(obj.azimuth_angle_of_departure_LOS_once(:,:,pp) - bearing_angle)+...
                                                  cosd(mechanical_slant).*sind(obj.azimuth_angle_of_departure_LOS_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_departure_LOS_once(:,:,pp))));
               phi_global_LOS = phi_global_LOS.*(180/pi);
               phi_global_LOS = wrapTo180(phi_global_LOS);  
               [theta_antenna_field_tx_LOS,phi_antenna_field_tx_LOS] = attached_eNodeB.antenna.polarization(theta_global_LOS,phi_global_LOS,slant_angle);
               % Antenna field pattern in elevation and azimuth in GCS                                           
               theta_AntennaField_tx_global_LOS = theta_antenna_field_tx_LOS.*cosd(psi_angle_LOS) - phi_antenna_field_tx_LOS.*sind(psi_angle_LOS);
               phi_AntennaField_tx_global_LOS = theta_antenna_field_tx_LOS.*sind(psi_angle_LOS) + phi_antenna_field_tx_LOS.*cosd(psi_angle_LOS);
               
           end
           
           %NLOS
           function [theta_AntennaField_tx_global_NLOS, phi_AntennaField_tx_global_NLOS] = eNodeB_antenna_field_pattern_GCS_NLOS(obj,LTE_config,attached_eNodeB,slant_angle,pp) 
               % Translates the antenna element field pattern from Local Coordinate System (LCS) in Global Coordinate System (GCS) 
               % See Section 5.1.3 in TR 36.873
               %
               % input:                        theta_antenna_field_tx_NLOS        ... polarized antenna element field pattern in zenith in LCS
               %                               phi_antenna_field_tx_NLOS          ... polarized antenna element field pattern in azimuth in LCS
               %                               bearing_angle                      ... pointing direction of antenna relative to x axis (0° in +x direction)
               %                               mechanical_downtilt                ... direction of antenna element relative to the y axis (0° in +y direction) 
               %                               mechanical_slant                  ... direction of antenna element relative to the z axis (0° in +z direction) used for rotation 
               %                               slant_angle                        ... the angle of slanted antenna element relative to the z axis (0° in +z direction) 
               %                                                                      depends on the plarization mode (defined in the config)
               %
               % output:                       theta_AntennaField_tx_global_NLOS   ... polarized antenna element field pattern in zenith in GCS
               %                               phi_AntennaField_tx_global_NLOS     ... polarized antenna element field pattern in azimuth in GCS
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               bearing_angle = obj.BS_boresight;% sector boresight angle with respect to x-axis 
               mechanical_downtilt = LTE_config.mechanical_downtilt; % mechanical downtilt between 0° and 180°
               mechanical_slant = LTE_config.mechanical_slant; %mechanical slant angle assumed to be zero
               
               psi_angle_NLOS = angle(sind(mechanical_slant).*cosd(obj.zenith_angle_of_departure_NLOS_once(:,:,pp)).*sind(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp) - bearing_angle)+...
                                     cosd(mechanical_slant).*(cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_departure_NLOS_once(:,:,pp)) - sind(mechanical_downtilt).*...
                                     cosd(obj.zenith_angle_of_departure_NLOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp) - bearing_angle))+...
                                         1i.*(sind(mechanical_slant).*cosd(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp) - bearing_angle)+...
                                         sind(mechanical_downtilt).*cosd(mechanical_slant).*sind(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp) - bearing_angle)));
               psi_angle_NLOS = psi_angle_NLOS.*(180/pi);   
               
               theta_global_NLOS = acosd(cosd(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.zenith_angle_of_departure_NLOS_once(:,:,pp))+...
                                     (sind(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp) - bearing_angle)-...
                                       sind(mechanical_slant).*sin(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_departure_NLOS_once(:,:,pp)));
                        
               phi_global_NLOS = angle((cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_departure_NLOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp) - bearing_angle)-...
                                              sind(mechanical_downtilt).*cosd(obj.zenith_angle_of_departure_NLOS_once(:,:,pp)))+...
                                              1i.*((cosd(mechanical_downtilt )).*sind(mechanical_slant).*cosd(obj.zenith_angle_of_departure_NLOS_once(:,:,pp))+...
                                                 (sind(mechanical_downtilt ).*sind(mechanical_slant).*cosd(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp) - bearing_angle)+...
                                                  cosd(mechanical_slant).*sind(obj.azimuth_angle_of_departure_NLOS_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_departure_NLOS_once(:,:,pp))));
               phi_global_NLOS = phi_global_NLOS.*(180/pi);                               
               phi_global_NLOS = wrapTo180(phi_global_NLOS);                                                                                         
               [theta_antenna_field_tx_NLOS,phi_antenna_field_tx_NLOS]  = attached_eNodeB.antenna.polarization(theta_global_NLOS,phi_global_NLOS,slant_angle);                     
               % Antenna field pattern in elevation and azimuth in GCS                                             
               theta_AntennaField_tx_global_NLOS = theta_antenna_field_tx_NLOS.*cosd(psi_angle_NLOS) - phi_antenna_field_tx_NLOS.*sind(psi_angle_NLOS);
               phi_AntennaField_tx_global_NLOS = theta_antenna_field_tx_NLOS.*sind(psi_angle_NLOS) + phi_antenna_field_tx_NLOS.*cosd(psi_angle_NLOS);
               
           end
           
           function [theta_AntennaField_tx_global_OTOI, phi_AntennaField_tx_global_OTOI] = eNodeB_antenna_field_pattern_GCS_OTOI(obj,LTE_config,attached_eNodeB,slant_angle,pp) 
               % Translates the antenna element field pattern from Local Coordinate System (LCS) in Global Coordinate System (GCS) 
               % See Section 5.1.3 in TR 36.873
               %
               % input:                        theta_antenna_field_tx_OTOI        ... polarized antenna element field pattern in zenith in LCS
               %                               phi_antenna_field_tx_OTOI          ... polarized antenna element field pattern in azimuth in LCS
               %                               bearing_angle                      ... pointing direction of antenna relative to x axis (0° in +x direction)
               %                               mechanical_downtilt                ... direction of antenna element relative to the y axis (0° in +y direction)  
               %                               mechanical_slant                  ... direction of antenna element relative to the z axis (0° in +z direction) used for rotation 
               %                               slant_angle                        ... the angle of slanted antenna element relative to the z axis (0° in +z direction) 
               %                                                                      depends on the plarization mode (defined in the config)
               %
               % output:                       theta_AntennaField_tx_global_OTOI  ... polarized antenna element field pattern in zenith in GCS
               %                               phi_AntennaField_tx_global_OTOI    ... polarized antenna element field pattern in azimuth in GCS
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               bearing_angle = obj.BS_boresight;% sector boresight angle with respect to x-axis 
               mechanical_downtilt = LTE_config.mechanical_downtilt; % mechanical downtilt between 0° and 180°
               mechanical_slant = LTE_config.mechanical_slant; %mechanical slant angle assumed to be zero
               
               psi_angle_OTOI = angle(sind(mechanical_slant).*cosd(obj.zenith_angle_of_departure_OTOI_once(:,:,pp)).*sind(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp) - bearing_angle)+...
                                     cosd(mechanical_slant).*(cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_departure_OTOI_once(:,:,pp)) - sind(mechanical_downtilt).*...
                                     cosd(obj.zenith_angle_of_departure_OTOI_once(:,:,pp)).*cosd(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp) - bearing_angle))+...
                                         1i.*(sind(mechanical_slant).*cosd(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp) - bearing_angle)+...
                                         sind(mechanical_downtilt).*cosd(mechanical_slant).*sind(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp) - bearing_angle)));
               psi_angle_OTOI = psi_angle_OTOI.*(180/pi);   
               
               theta_global_OTOI = acosd(cosd(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.zenith_angle_of_departure_OTOI_once(:,:,pp))+...
                                     (sind(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp) - bearing_angle)-...
                                       sind(mechanical_slant).*sin(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_departure_OTOI_once(:,:,pp)));
                        
               phi_global_OTOI = angle((cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_departure_OTOI_once(:,:,pp)).*cosd(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp) - bearing_angle)-...
                                              sind(mechanical_downtilt).*cosd(obj.zenith_angle_of_departure_OTOI_once(:,:,pp)))+...
                                              1i.*((cosd(mechanical_downtilt )).*sind(mechanical_slant).*cosd(obj.zenith_angle_of_departure_OTOI_once(:,:,pp))+...
                                                 (sind(mechanical_downtilt ).*sind(mechanical_slant).*cosd(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp) - bearing_angle)+...
                                                  cosd(mechanical_slant).*sind(obj.azimuth_angle_of_departure_OTOI_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_departure_OTOI_once(:,:,pp))));
               phi_global_OTOI = phi_global_OTOI.*(180/pi);                               
               phi_global_OTOI = wrapTo180(phi_global_OTOI);                                                                                         
               [theta_antenna_field_tx_OTOI,phi_antenna_field_tx_OTOI]  = attached_eNodeB.antenna.polarization(theta_global_OTOI,phi_global_OTOI,slant_angle);                     
               % Antenna field pattern in elevation and azimuth in GCS                                             
               theta_AntennaField_tx_global_OTOI = theta_antenna_field_tx_OTOI.*cosd(psi_angle_OTOI) - phi_antenna_field_tx_OTOI.*sind(psi_angle_OTOI);
               phi_AntennaField_tx_global_OTOI = theta_antenna_field_tx_OTOI.*sind(psi_angle_OTOI) + phi_antenna_field_tx_OTOI.*cosd(psi_angle_OTOI);
               
           end
           
           function [theta_AntennaField_tx_global_LOS_direct_ray, phi_AntennaField_tx_global_LOS_direct_ray] = antenna_field_pattern_GCS_LOS_direct_ray(obj,LTE_config,attached_eNodeB,slant_angle,pp)
               % Translates the antenna element field pattern from Local
               % Coordinate System (LCS) in Global Coordinate System (GCS)
               % for the direct LOS path
               % See Section 5.1.3 in TR 36.873
               %
               % input:                   theta_antenna_field_tx_LOS_direct_ray ... polarized antenna element field pattern in zenith in LCS
               %                          phi_antenna_field_tx_LOS_direct_ray   ... polarized antenna element field pattern in azimuth in LCS
               %                          bearing_angle                         ... pointing direction of antenna relative to x axis (0° in +x direction)
               %                          mechanical_downtilt                   ... direction of antenna element relative to the y axis (0° in +y direction)  
               %                          mechanical_slant                      ... direction of antenna element relative to the z axis (0° in +z direction) used for rotation
               %                          slant_angle                           ... the angle of slanted antenna element relative to the z axis (0° in +z direction) 
               %                                                                    depends on the plarization mode (defined in the config)
               %
               % output:                  theta_AntennaField_tx_global_LOS_direct_ray ... polarized antenna element field pattern in zenith in GCS
               %                          phi_AntennaField_tx_global_LOS_direct_ray   ... polarized antenna element field pattern in azimuth in GCS
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               bearing_angle = obj.BS_boresight;% sector boresight angle with respect to x-axis 
               mechanical_downtilt = LTE_config.mechanical_downtilt; % mechanical downtilt between 0° and 180°
               mechanical_slant = LTE_config.mechanical_slant; %mechanical slant angle assumed to be zero
               
               psi_angle_LOS_direct_ray = angle(sind(mechanical_slant)*cosd(obj.theta_departure(pp))*sind(obj.phi_departure(pp) - bearing_angle)+...
                                           cosd(mechanical_slant)*(cosd(mechanical_downtilt)*sind(obj.theta_departure(pp)) - sind(mechanical_downtilt).*...
                                           cosd(obj.theta_departure(pp)).*cosd(obj.phi_departure(pp) - bearing_angle))+...
                                            1i.*(sind(mechanical_slant)*cosd(obj.phi_departure(pp) - bearing_angle)+...
                                                 sind(mechanical_downtilt )*cosd(mechanical_slant)*sind(obj.phi_departure(pp)-bearing_angle)));
               psi_angle_LOS_direct_ray = psi_angle_LOS_direct_ray.*(180/pi);
               
               theta_global_LOS_direct_ray = acosd(cosd(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.theta_departure(pp))+...
                                               (sind(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.phi_departure(pp) - bearing_angle)-...
                                                sind(mechanical_slant).*sin(obj.phi_departure(pp) - bearing_angle)).*sind(obj.theta_departure(pp)));
                                       
               phi_global_LOS_direct_ray = angle((cosd(mechanical_downtilt ).*sind(obj.theta_departure(pp)).*cosd(obj.phi_departure(pp)-bearing_angle)-...
                                              sind(mechanical_downtilt ).*cosd(obj.theta_departure(pp)))+...
                                              1i.*((cosd(mechanical_downtilt)).*sind(mechanical_slant).*cosd(obj.theta_departure(pp))+...
                                                 (sind(mechanical_downtilt ).*sind(mechanical_slant).*cosd(obj.phi_departure(pp)-bearing_angle)+...
                                                  cosd(mechanical_slant).*sind(obj.phi_departure(pp)-bearing_angle)).*sind(obj.theta_departure(pp))));
               phi_global_LOS_direct_ray = phi_global_LOS_direct_ray.*(180/pi);                               
               phi_global_LOS_direct_ray = wrapTo180(phi_global_LOS_direct_ray);                                                             
               [theta_antenna_field_tx_LOS_direct_ray,phi_antenna_field_tx_LOS_direct_ray] = attached_eNodeB.antenna.polarization(theta_global_LOS_direct_ray,phi_global_LOS_direct_ray,slant_angle); 
               % Antenna field pattern in elevation and azimuth in GCS                                              
               theta_AntennaField_tx_global_LOS_direct_ray = theta_antenna_field_tx_LOS_direct_ray.*cosd(psi_angle_LOS_direct_ray) - phi_antenna_field_tx_LOS_direct_ray.*sind(psi_angle_LOS_direct_ray);
               phi_AntennaField_tx_global_LOS_direct_ray = theta_antenna_field_tx_LOS_direct_ray.*sind(psi_angle_LOS_direct_ray) + phi_antenna_field_tx_LOS_direct_ray.*cosd(psi_angle_LOS_direct_ray);
               
           end

  
           function [UE_theta_AntennaField_global_LOS,UE_phi_AntennaField_global_LOS] = UE_antenna_field_pattern_global_LOS(obj,UE_slant_angle,direction_of_movement,pp)  
               % Translates the antenna element field pattern at the receiver from Local Coordinate System (LCS) in Global Coordinate System (GCS) 
               % See Section 5.1.3 in TR 36.873
               %
               % input:                        UE_theta_antenna_field_LOS        ... user polarized antenna element field pattern in zenith in LCS
               %                               UE_phi_antenna_field_LOS          ... user polarized antenna element field pattern in azimuth in LCS
               %                               bearing_angle                     ... pointing direction of antenna at user relative to x axis (0° in +x direction)
               %                                                                     a random angle uniformly distributed between 0° and 360°
               %                               mechanical_downtilt               ... direction of antenna element relative to the y axis (0° in +y direction)   
               %                               mechanical_slant                  ... direction of antenna element relative to the z axis (0° in +z direction) 
               %                               UE_slant_angle                    ... the angle of slanted antenna element at user relative to the z axis (0° in +z direction) 
               %                                                                     depends on the plarization mode (defined in the config)
               %
               % output:                       UE_theta_AntennaField_global_LOS  ... user polarized antenna element field pattern in zenith in GCS
               %                               UE_phi_AntennaField_global_LOS    ... user polarized antenna element field pattern in azimuth in GCS
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               bearing_angle = direction_of_movement; % random angle uniformly distributed between 0° and 360°
               mechanical_downtilt = 90; % downtilt with respect to y-axis
               mechanical_slant    = 0; % 3GPP TR 36.873 Table 8.2-2

               
              UE_psi_angle_LOS = angle(sind(mechanical_slant).*cosd(obj.zenith_angle_of_arrival_LOS_once(:,:,pp)).*sind(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp) - bearing_angle)+...
                   cosd(mechanical_slant).*(cosd(mechanical_downtilt ).*sind(obj.zenith_angle_of_arrival_LOS_once(:,:,pp)) - sind(mechanical_downtilt).*...
                   cosd(obj.zenith_angle_of_arrival_LOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp) - bearing_angle))+...
                   1i.*(sind(mechanical_slant).*cosd(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp) - bearing_angle)+...
                   sind(mechanical_downtilt).*cosd(mechanical_slant).*sind(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp) - bearing_angle)));
               UE_psi_angle_LOS = UE_psi_angle_LOS.*(180/pi);
               
               UE_theta_global_LOS = acosd(cosd(mechanical_downtilt ).*cosd(mechanical_slant).*cosd(obj.zenith_angle_of_arrival_LOS_once(:,:,pp))+...
                   (sind(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp) - bearing_angle)-...
                   sind(mechanical_slant).*sin(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_arrival_LOS_once(:,:,pp)));
               
               UE_phi_global_LOS = angle((cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_arrival_LOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp) - bearing_angle)-...
                   sind(mechanical_downtilt).*cosd(obj.zenith_angle_of_arrival_LOS_once(:,:,pp)))+...
                   1i.*((cosd(mechanical_downtilt)).*sind(mechanical_slant).*cosd(obj.zenith_angle_of_arrival_LOS_once(:,:,pp))+...
                   (sind(mechanical_downtilt).*sind(mechanical_slant).*cosd(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp) -bearing_angle)+...
                   cosd(mechanical_slant).*sind(obj.azimuth_angle_of_arrival_LOS_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_arrival_LOS_once(:,:,pp))));
               UE_phi_global_LOS = UE_phi_global_LOS.*(180/pi);
               
               % Assume omnidirectional antenna at UE
               UE_theta_antenna_field_LOS = ones(size(UE_theta_global_LOS))*cosd(UE_slant_angle);
               UE_phi_antenna_field_LOS = ones(size(UE_phi_global_LOS))*sind(UE_slant_angle);
               UE_theta_AntennaField_global_LOS = UE_theta_antenna_field_LOS.*cosd(UE_psi_angle_LOS) - UE_phi_antenna_field_LOS.*sind(UE_psi_angle_LOS);
               UE_phi_AntennaField_global_LOS = UE_theta_antenna_field_LOS.*sind(UE_psi_angle_LOS) + UE_phi_antenna_field_LOS.*cosd(UE_psi_angle_LOS);
               
           end
           
           %NLOS
           function [UE_theta_AntennaField_global_NLOS, UE_phi_AntennaField_global_NLOS] = UE_antenna_field_pattern_global_NLOS(obj,UE_slant_angle,direction_of_movement,pp)
               % Translates the antenna element field pattern at the receiver from Local Coordinate System (LCS) in Global Coordinate System (GCS) 
               % See Section 5.1.3 in TR 36.873
               %
               % input:                        UE_theta_antenna_field_NLOS        ... user polarized antenna element field pattern in zenith in LCS
               %                               UE_phi_antenna_field_NLOS          ... user polarized antenna element field pattern in azimuth in LCS
               %                               bearing_angle                     ... pointing direction of antenna at user relative to x axis (0° in +x direction)
               %                                                                      a random angle uniformly distributed between 0° and 360°
               %                               mechanical_downtilt               ... direction of antenna element relative to the y axis (0° in +y direction)   
               %                               mechanical_slant                  ... direction of antenna element relative to the z axis (0° in +z direction) 
               %                               UE_slant_angle                    ... the angle of slanted antenna element at user relative to the z axis (0° in +z direction) 
               %                                                                     depends on the plarization mode (defined in the config)
               %
               % output:                       UE_theta_AntennaField_global_NLOS  ... user polarized antenna element field pattern in zenith in GCS
               %                               UE_phi_AntennaField_global_NLOS    ... user polarized antenna element field pattern in azimuth in GCS
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               bearing_angle = direction_of_movement; % random angle uniformly distributed between 0° and 360°
               mechanical_downtilt = 90; % downtilt with respect to y-axis
               mechanical_slant    = 0; % 3GPP TR 36.873 Table 8.2-2
                
               UE_psi_angle_NLOS = angle(sind(mechanical_slant).*cosd(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp)).*sind(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp) - bearing_angle)+...
                   cosd(mechanical_slant).*(cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp)) - sind(mechanical_downtilt).*...
                   cosd(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp) - bearing_angle))+...
                   1i.*(sind(mechanical_slant).*cosd(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp) - bearing_angle)+...
                   sind(mechanical_downtilt ).*cosd(mechanical_slant).*sind(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp) - bearing_angle)));
               UE_psi_angle_NLOS = UE_psi_angle_NLOS.*(180/pi);
               
               UE_theta_global_NLOS = acosd(cosd(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp))+...
                   (sind(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp) - bearing_angle)-...
                   sind(mechanical_slant).*sin(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp)));
               
               UE_phi_global_NLOS = angle((cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp)).*cosd(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp) - bearing_angle)-...
                   sind(mechanical_downtilt).*cosd(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp)))+...
                   1i.*((cosd(mechanical_downtilt)).*sind(mechanical_slant).*cosd(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp))+...
                   (sind(mechanical_downtilt).*sind(mechanical_slant).*cosd(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp) - bearing_angle)+...
                   cosd(mechanical_slant).*sind(obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_arrival_NLOS_once(:,:,pp))));
               UE_phi_global_NLOS = UE_phi_global_NLOS.*180/pi;
               
               % Assume omnidirectional antenna at UE
               UE_theta_antenna_field_NLOS = ones(size(UE_theta_global_NLOS))*cosd(UE_slant_angle);
               UE_phi_antenna_field_NLOS = ones(size(UE_phi_global_NLOS))*sind(UE_slant_angle);
               UE_theta_AntennaField_global_NLOS = UE_theta_antenna_field_NLOS.*cosd(UE_psi_angle_NLOS) - UE_phi_antenna_field_NLOS.*sind(UE_psi_angle_NLOS);
               UE_phi_AntennaField_global_NLOS = UE_theta_antenna_field_NLOS.*sind(UE_psi_angle_NLOS) + UE_phi_antenna_field_NLOS.*cosd(UE_psi_angle_NLOS);
                
           end
           
           
            function [UE_theta_AntennaField_global_OTOI, UE_phi_AntennaField_global_OTOI] = UE_antenna_field_pattern_global_OTOI(obj,UE_slant_angle,direction_of_movement,pp)
               % Translates the antenna element field pattern at the receiver from Local Coordinate System (LCS) in Global Coordinate System (GCS) 
               % See Section 5.1.3 in TR 36.873
               %
               % input:                        UE_theta_antenna_field_OTOI        ... user polarized antenna element field pattern in zenith in LCS
               %                               UE_phi_antenna_field_OTOI          ... user polarized antenna element field pattern in azimuth in LCS
               %                               bearing_angle                     ... pointing direction of antenna at user relative to x axis (0° in +x direction)
               %                                                                      a random angle uniformly distributed between 0° and 360°
               %                               mechanical_downtilt               ... direction of antenna element relative to the y axis (0° in +y direction)
               %                               mechanical_slant                  ... direction of antenna element relative to the z axis (0° in +z direction)
               %                               UE_slant_angle                    ... the angle of slanted antenna element at user relative to the z axis (0° in +z direction) 
               %                                                                     depends on the plarization mode (defined in the config)
               %
               % output:                       UE_theta_AntennaField_global_OTOI  ... user polarized antenna element field pattern in zenith in GCS
               %                               UE_phi_AntennaField_global_OTOI    ... user polarized antenna element field pattern in azimuth in GCS
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
                bearing_angle = direction_of_movement; % random angle uniformly distributed between 0° and 360°
                mechanical_downtilt = 90; % downtilt with respect to y-axis
                mechanical_slant    = 0; % 3GPP TR 36.873 Table 8.2-2
                 
                 UE_psi_angle_OTOI = angle(sind(mechanical_slant).*cosd(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp)).*sind(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp) - bearing_angle)+...
                    cosd(mechanical_slant).*(cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp)) - sind(mechanical_downtilt).*...
                    cosd(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp)).*cosd(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp) - bearing_angle))+...
                    1i.*(sind(mechanical_slant).*cosd(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp) - bearing_angle)+...
                    sind(mechanical_downtilt ).*cosd(mechanical_slant).*sind(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp) - bearing_angle)));
                UE_psi_angle_OTOI = UE_psi_angle_OTOI.*(180/pi);
                
                UE_theta_global_OTOI = acosd(cosd(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp))+...
                    (sind(mechanical_downtilt).*cosd(mechanical_slant).*cosd(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp) - bearing_angle)-...
                    sind(mechanical_slant).*sin(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp)));
                
                UE_phi_global_OTOI = angle((cosd(mechanical_downtilt).*sind(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp)).*cosd(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp) - bearing_angle)-...
                    sind(mechanical_downtilt).*cosd(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp)))+...
                    1i.*((cosd(mechanical_downtilt)).*sind(mechanical_slant).*cosd(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp))+...
                    (sind(mechanical_downtilt).*sind(mechanical_slant).*cosd(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp) - bearing_angle)+...
                    cosd(mechanical_slant).*sind(obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp) - bearing_angle)).*sind(obj.zenith_angle_of_arrival_OTOI_once(:,:,pp))));
                UE_phi_global_OTOI = UE_phi_global_OTOI.*180/pi;
                
                % Assume omnidirectional antenna at UE
                UE_theta_antenna_field_OTOI = ones(size(UE_theta_global_OTOI))*cosd(UE_slant_angle);
                UE_phi_antenna_field_OTOI = ones(size(UE_phi_global_OTOI))*sind(UE_slant_angle);
                UE_theta_AntennaField_global_OTOI = UE_theta_antenna_field_OTOI.*cosd(UE_psi_angle_OTOI) - UE_phi_antenna_field_OTOI.*sind(UE_psi_angle_OTOI);
                UE_phi_AntennaField_global_OTOI = UE_theta_antenna_field_OTOI.*sind(UE_psi_angle_OTOI) + UE_phi_antenna_field_OTOI.*cosd(UE_psi_angle_OTOI);

           end
          
           function [UE_theta_AntennaField_global_LOS_direct_ray,UE_phi_AntennaField_global_LOS_direct_ray] = UE_antenna_field_pattern_global_LOS_direct_ray(obj,UE_slant_angle,direction_of_movement,pp)
               % Translates the antenna element field pattern at the receiver from Local Coordinate System (LCS) in Global Coordinate System (GCS) 
               % for the direct LOS path
               % See Section 5.1.3 in TR 36.873
               %
               % input:           UE_theta_AntennaField_global_LOS_direct_ray  ... user polarized antenna element field pattern in zenith in LCS
               %                  UE_phi_AntennaField_global_LOS_direct_ray    ... user polarized antenna element field pattern in azimuth in LCS
               %                  bearing_angle                                ... pointing direction of antenna at user relative to x axis (0° in +x direction)
               %                                                                   a random angle uniformly distributed between 0° and 360°
               %                  mechanical_downtilt                          ... direction of antenna element relative to the y axis (0° in +y direction)   
               %                  mechanical_slant                             ... direction of antenna element relative to the z axis (0° in +z direction)
               %                  UE_slant_angle                               ... the angle of slanted antenna element at user relative to the z axis (0° in +z direction) 
               %                                                                     depends on the plarization mode (defined in the config)
               %
               % output:          UE_theta_AntennaField_global_LOS_direct_ray  ... user polarized antenna element field pattern in zenith in GCS
               %                  UE_phi_AntennaField_global_LOS_direct_ray    ... user polarized antenna element field pattern in azimuth in GCS
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               bearing_angle = direction_of_movement; % random angle uniformly distributed between 0° and 360°
               mechanical_downtilt = 90; % downtilt with respect to y-axis
               mechanical_slant    = 0; % 3GPP TR 36.873 Table 8.2-2
 
              UE_psi_angle_LOS_direct_ray = angle(sind(mechanical_slant).*cosd(obj.theta_arrival(pp)).*sind(obj.phi_arrival(pp) - bearing_angle)+...
                   cosd(mechanical_slant).*(cosd(mechanical_downtilt).*sind(obj.theta_arrival(pp)) - sind(mechanical_downtilt).*...
                   cosd(obj.theta_arrival(pp)).*cosd(obj.phi_arrival(pp) - bearing_angle))+...
                   1i.*(sind(mechanical_slant).*cosd(obj.phi_arrival(pp) - bearing_angle)+...
                   sind(mechanical_downtilt).*cosd(mechanical_slant ).*sind(obj.phi_arrival(pp)-bearing_angle)));
               UE_psi_angle_LOS_direct_ray = UE_psi_angle_LOS_direct_ray.*(180/pi);
               
               UE_theta_global_LOS_direct_ray = acosd(cosd(mechanical_downtilt).*cosd(mechanical_slant ).*cosd(obj.theta_arrival(pp))+...
                   (sind(mechanical_downtilt).*cosd(mechanical_slant ).*cosd(obj.phi_arrival(pp) - bearing_angle)-...
                   sind(mechanical_slant ).*sin(obj.phi_arrival(pp) - bearing_angle)).*sind(obj.theta_arrival(pp)));
               
               UE_phi_global_LOS_direct_ray = angle((cosd(mechanical_downtilt).*sind(obj.theta_arrival(pp)).*cosd(obj.phi_arrival(pp)-bearing_angle)-...
                   sind(mechanical_downtilt ).*cosd(obj.theta_arrival(pp)))+...
                   1i.*((cosd(mechanical_downtilt)).*sind(mechanical_slant ).*cosd(obj.theta_arrival(pp))+...
                   (sind(mechanical_downtilt).*sind(mechanical_slant ).*cosd(obj.phi_arrival(pp)-bearing_angle)+...
                   cosd(mechanical_slant).*sind(obj.phi_arrival(pp)-bearing_angle)).*sind(obj.theta_arrival(pp))));
               UE_phi_global_LOS_direct_ray = UE_phi_global_LOS_direct_ray.*(180/pi);
                
               % Assume omni-directional antenna at UE
               UE_theta_AntennaField_global_LOS_direct_ray = ones(size(UE_theta_global_LOS_direct_ray))*cosd(UE_slant_angle);
               UE_phi_AntennaField_global_LOS_direct_ray = ones(size(UE_phi_global_LOS_direct_ray))*sind(UE_slant_angle);
               UE_theta_AntennaField_global_LOS_direct_ray = UE_theta_AntennaField_global_LOS_direct_ray.*cosd(UE_psi_angle_LOS_direct_ray) - UE_phi_AntennaField_global_LOS_direct_ray.*sind(UE_psi_angle_LOS_direct_ray);
               UE_phi_AntennaField_global_LOS_direct_ray = UE_theta_AntennaField_global_LOS_direct_ray.*sind(UE_psi_angle_LOS_direct_ray) + UE_phi_AntennaField_global_LOS_direct_ray.*cosd(UE_psi_angle_LOS_direct_ray);
               
           end


           
           %% O- to- I scenario
           % Calculate spherical unit vectors, Doppler and antenna field
           % This function is only called once per antenna array
           
           function field_and_phases_over_rays_and_clusters = calculate_initial_clusters_rays_phases_and_fields_OTOI(obj, LTE_config, rx_height, attached_site_pos, UE_pos, sigmas_OTOI, ZOD_parameters, attached_eNodeB, nTx, nRx, UE_is_indoor, direction_of_movement, pp, dist_indoor)
               % Generates only the field multiplication of the transmitter 
               % and receiver from eq 7.21 in TR 36.873
               % without the Doppler component and exponential terms
               % This function is only called once per antenna array
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               [obj.cl_power_OTOI ,obj.cl_power_OTOI_general, obj.cl_power_OTOI_per_ray] = obj.generate_cluster_power_OTOI(sigmas_OTOI);
               obj.xpr_OTOI = obj.generate_XPR_OTOI; % Cross polarization powers.
               [obj.theta_arrival(pp),obj.theta_departure(pp),obj.phi_arrival(pp),obj.phi_departure(pp)]= obj.eNodeB_UE_LOS_direction_angles_signal_OTOI(LTE_config,rx_height,attached_site_pos,UE_pos, dist_indoor);
               [obj.azimuth_angle_of_departure_OTOI_once(:,:,pp),...
                   obj.azimuth_angle_of_arrival_OTOI_once(:,:,pp),...
                   obj.zenith_angle_of_departure_OTOI_once(:,:,pp),...
                   obj.zenith_angle_of_arrival_OTOI_once(:,:,pp)]       = obj.coupling_rays_OTOI(sigmas_OTOI,ZOD_parameters,obj.cl_power_OTOI,UE_is_indoor,pp); 
               obj.RX_spherical_unit_vector_OTOI(:,:,:,pp)              = obj.receiver_spherical_unit_vector_OTOI(pp);
               obj.TX_spherical_unit_vector_OTOI(:,:,:,pp)              = obj.transmitter_spherical_unit_vector_OTOI(pp);
               obj.Doppler_speed_OTOI(:,:,:,pp)                         = obj.calculate_doppler_speed_OTOI(LTE_config, direction_of_movement, pp);
               initial_phases_OTOI                                      = obj.generate_initial_phases_OTOI;
               initial_phases_OTOI                                      = initial_phases_OTOI *pi/180;  % conversion to rad required
               obj.BS_boresight                                         = attached_eNodeB.boresight;
               % Antenna field
               for uu = 1:nRx
                   % Determine slant angle of UE antenna
                   if strcmp(LTE_config.UE_antenna_polarization,'ULA')
                       UE_slant_angle = LTE_config.UE_antenna_slant_angle;
                   else
                       UE_slant_angle = LTE_config.UE_antenna_slant_angle*mod(uu,2);
                   end
                   
                   for nn = 1:nTx
                       % Determine slant angle of antenna array column: 0...linearly polarized,  +45/-45...x-polarized
                       if strcmp(LTE_config.antenna.antenna_polarization,'XPOL')
                           slant_angle = LTE_config.antenna.slant_angle*(2*mod(nn,2)-1);
                       else
                           slant_angle = LTE_config.antenna.slant_angle;
                       end
                       [obj.UE_theta_AntennaField_global_OTOI(:,:,nn,uu,pp), obj.UE_phi_AntennaField_global_OTOI(:,:,nn,uu,pp)] = ...
                           obj.UE_antenna_field_pattern_global_OTOI(UE_slant_angle,direction_of_movement,pp);
                       [obj.theta_AntennaField_tx_global_OTOI(:,:,nn,pp), obj.phi_AntennaField_tx_global_OTOI(:,:,nn,pp)] = ...
                           obj.eNodeB_antenna_field_pattern_GCS_OTOI(LTE_config,attached_eNodeB,slant_angle,pp);
                       %part form generate_initial_phases_and_fields
                       f_matrix = (exp(1i.*initial_phases_OTOI(:,:,1)).* obj.theta_AntennaField_tx_global_OTOI(:,:,nn,pp) + sqrt(obj.xpr_OTOI.^(-1)).*exp(1i.*initial_phases_OTOI(:,:,2)).*obj.phi_AntennaField_tx_global_OTOI(:,:,nn,pp)) .* obj.UE_theta_AntennaField_global_OTOI(:,:,nn,uu,pp) +...
                           (sqrt(obj.xpr_OTOI.^(-1)).*exp(1i.*initial_phases_OTOI(:,:,3)).*obj.theta_AntennaField_tx_global_OTOI(:,:,nn,pp) + exp(1i.*initial_phases_OTOI(:,:,4)).*obj.phi_AntennaField_tx_global_OTOI(:,:,nn,pp)).* obj.UE_phi_AntennaField_global_OTOI(:,:,nn,uu,pp);
                       
                       field_and_phases_over_rays_and_clusters(:,:,nn,uu) = sqrt(obj.cl_power_OTOI_per_ray).*f_matrix;
                   end
               end
           end
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           %%%% USED FOR CALIBRATION ONLY
            % Save NLOS angles 
            function azimuth_of_departure_OTOI = azimuth_departure_angles_OTOI(obj)
               azimuth_of_departure_OTOI = obj.azimuth_angle_of_departure_OTOI_once(:,:,1);
            end
            
           function zenith_of_departure_OTOI = zenith_departure_angles_OTOI(obj)
               zenith_of_departure_OTOI  = obj.zenith_angle_of_departure_OTOI_once(:,:,1);
           end
            
           function azimuth_of_arrival_OTOI = azimuth_arrival_angles_OTOI(obj)
               azimuth_of_arrival_OTOI   = obj.azimuth_angle_of_arrival_OTOI_once(:,:,1);
           end
           function zenith_of_arrival_OTOI = zenith_arrival_angles_OTOI(obj)
               zenith_of_arrival_OTOI    = obj.zenith_angle_of_arrival_OTOI_once(:,:,1);
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
           function exp_part_without_Doppler_OTOI = calculate_exp_part_without_Doppler_OTOI(obj,LTE_config,rx_height,attached_site_pos,UE_pos,relative_antenna_position,relative_UE_antenna_position, pp)
               % Generates the exponential term in eq. 7.21 without the
               % Doppler component
               % A product of exponential terms at the transmitter and
               % receiver, where each exponential term is a product 
               % between receiver/transmitter spherical unit vector and 
               % the receiver/transmitter location vector
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               
               wavelength                         = 299792458./LTE_config.frequency;
               UE_location                        = obj.UE_position(rx_height,UE_pos,relative_UE_antenna_position);
               eNodeB_location                    = obj.eNodeB_position(LTE_config,attached_site_pos,relative_antenna_position);
               % Product between receiver spherical unit vector and the
               % receiver location vector
               size_vec = size(obj.RX_spherical_unit_vector_OTOI);
               RX_vector_temp = reshape(obj.RX_spherical_unit_vector_OTOI(:,:,:,pp),3,size_vec(2)*size_vec(3));
               calculate_r_rx_and_location_vector = reshape(RX_vector_temp.'*UE_location.',size_vec(2),size_vec(3));
               % Product between transmitter spherical unit vector and the
               % transmitter location vector
               TX_vector_temp = reshape(obj.TX_spherical_unit_vector_OTOI(:,:,:,pp),3,size_vec(2)*size_vec(3));
               calculate_r_tx_and_location_vector = reshape(TX_vector_temp.'*eNodeB_location.',size_vec(2),size_vec(3));

               % Calculate exponential part of channel matrix
               exp_part_without_Doppler_OTOI   = exp(1i*2*pi/wavelength.*calculate_r_rx_and_location_vector).*exp(1i*2*pi/wavelength.*calculate_r_tx_and_location_vector);  
           end
         
           
           function Doppler_over_time_OTOI = calculate_Doppler_over_time_OTOI(obj, current_TTI, time_length, tSubframe, pp)
               % Generates the Doppler frequency component
               % The Doppler component is calculated for each time instant
               Doppler_over_time_OTOI = zeros(size(obj.Doppler_speed_OTOI,1),size(obj.Doppler_speed_OTOI,2), time_length);
               % Calculate Doopler effect for each time instant    
               for tt = current_TTI:time_length
                   Doppler_over_time_OTOI(:,:,tt-current_TTI+1)           = exp(1i*2*pi.*obj.Doppler_speed_OTOI(:,:,pp) * tt*tSubframe);
                   % Calibration, t=0, no Doppler
                   % Doppler_over_time_OTOI(:,:,tt-current_TTI+1)           = exp(0*obj.Doppler_speed_OTOI(:,:,pp) * tt*tSubframe);
               end
           end
           
       
           
           function channel_coefficients_over_clusters_OTOI = calculate_channel_coefficient_OTOI(obj,LTE_config,rx_height,attached_site_pos,UE_pos,sigmas_OTOI,ZOD_parameters,attached_eNodeB, nTx, nRx, UE_is_indoor, pp, current_TTI, time_length, tSubframe, direction_of_movement, dist_indoor)
               % Generates the channel impulse response for each antenna
               % element pair between tramsitter and receiver
               % 
               % output: channel_coefficients_over_clusters_OTOI  ... a matrix of dimensions:
               %                                                       - 1st dimension is number of receive antenna ports
               %                                                       - 2nd dimension is number of transmit antenna ports
               %                                                       - 3rd dimension is 1 (no reason, can be removed)
               %                                                       - 4th dimension is number of time blocks
               %                                                       - 5th dimension is number of clusters
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               wavelength                   = 299792458./LTE_config.frequency;
               % Check polarization of transmit and receive antennas - for horizontal spacing of antenna elements.
               if strcmp(LTE_config.antenna.antenna_polarization,'XPOL') tx_spacing_factor = 2; else  tx_spacing_factor = 1; end
               if strcmp(LTE_config.UE_antenna_polarization,'XPOL')      rx_spacing_factor = 2; else  rx_spacing_factor = 1; end
               %Generate channel parameters for antenna fields and phases
               %over rays and clusters
               obj.field_and_phases_over_rays_and_clusters_OTOI = obj.calculate_initial_clusters_rays_phases_and_fields_OTOI(LTE_config, rx_height, attached_site_pos, UE_pos, sigmas_OTOI, ZOD_parameters, attached_eNodeB, nTx, nRx, UE_is_indoor, direction_of_movement, pp, dist_indoor);
               obj.Doppler_over_time_OTOI = obj.calculate_Doppler_over_time_OTOI(current_TTI, time_length, tSubframe, pp);
               for uu = 1:nRx
                   relative_UE_antenna_position  = [0, LTE_config.UE_antenna_element_horizontal_spacing*(ceil(uu/rx_spacing_factor)-1), 0];
                   element_channel_coefficients = zeros(obj.PerClusterRays,obj.NumClusters_OTOI, nTx, 1, LTE_config.nr_of_antenna_elements_in_each_column, time_length-current_TTI+1); % use nr of clusters instead numbers rays, clusters
                   element_weight               = zeros(nTx, LTE_config.nr_of_antenna_elements_in_each_column);
                   
                   for nn = 1:nTx
                       % Generate channel coefficients for each antenna
                       % element
                       for mm = 1:LTE_config.nr_of_antenna_elements_in_each_column
                           relative_antenna_position = [0, LTE_config.antenna_element_horizontal_spacing*(ceil(nn/tx_spacing_factor)-1), LTE_config.antenna_element_vertical_spacing*(mm - 1)];
                           obj.exponential_part_without_Doppler_OTOI = obj.calculate_exp_part_without_Doppler_OTOI(LTE_config,rx_height,attached_site_pos,UE_pos,relative_antenna_position,relative_UE_antenna_position, pp);
                           element_weight(nn,mm) = 1/sqrt(LTE_config.nr_of_antenna_elements_in_each_column)*exp(-1i*2*pi/wavelength*(mm-1)*LTE_config.antenna_element_vertical_spacing*cosd(LTE_config.electrical_downtilt));
                           % perform simulation from current TTI, for all
                           % subsequent TTIs
                           for tt=current_TTI:time_length
                               element_channel_coefficients(:,:,nn, 1 ,mm,tt-current_TTI+1)= obj.field_and_phases_over_rays_and_clusters_OTOI(:,:,nn,uu)...
                                   .* obj.exponential_part_without_Doppler_OTOI .* obj.Doppler_over_time_OTOI(:,:,tt-current_TTI+1);
                               
                               element_channel_coeff_over_clusters = sum(element_channel_coefficients,1); % sum ch-coeff. over rays within each cluster
                               rearranged_element_channel_coefficent= permute(element_channel_coeff_over_clusters,[3 5 6 4 2 1]);
                               
                               for kk=1:obj.NumClusters_OTOI
                                   portwise_channel_coefficients(:,1,uu,tt-current_TTI+1,kk) = sum((rearranged_element_channel_coefficent(:,:,tt-current_TTI+1,1,kk).*element_weight),2);
                                   %obj.chan_norm_NLOS(nn,1,uu,tt) = sum(abs(portwise_channel_coefficients(nn,1,uu,tt,:)).^2);
                               end
                               chanel_coeff_norm_OTOI(nn,1,uu,tt-current_TTI+1) = sum(abs(portwise_channel_coefficients(nn,1,uu,tt-current_TTI+1,:)).^2);
                           end
                       end
                   end         
               end
               channel_coefficients_over_clusters_OTOI = permute(portwise_channel_coefficients, [3 1 2 4 5]);
               obj.channel_coefficients_over_clusters_OTOI = permute(portwise_channel_coefficients, [3 1 2 4 5]); % channel matrix with dimensions [nRx, nTx, 1  nTTI, nCluster]
               obj.chan_norm_OTOI = permute(chanel_coeff_norm_OTOI, [3 1 2 4]);
           end
           
           % Save channel matrix
           function channel_matrix_OTOI = saved_channel_matrix_OTOI(obj)
               channel_matrix_OTOI = obj.channel_coefficients_over_clusters_OTOI;
           end
%% NLOS 
           % Calculate spherical unit vectors, Doppler and antenna field
           % This function is only called once per antenna array
           
           function field_and_phases_over_rays_and_clusters = calculate_initial_clusters_rays_phases_and_fields_NLOS(obj, LTE_config, rx_height, attached_site_pos, UE_pos, sigmas_NLOS, ZOD_parameters, attached_eNodeB, nTx, nRx, direction_of_movement, pp)
               % Generates only the field multiplication of the transmitter 
               % and receiver from eq 7.21 in TR 36.873
               % without the Doppler component and exponential terms
               % This function is only called once per antenna array
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               [obj.cl_power_NLOS ,obj.cl_power_NLOS_general, obj.cl_power_NLOS_per_ray] = obj.generate_cluster_power_NLOS(sigmas_NLOS);
               obj.xpr_NLOS = obj.generate_XPR_NLOS; % Cross polarization powers.
               [obj.theta_arrival(pp),obj.theta_departure(pp),obj.phi_arrival(pp),obj.phi_departure(pp)]= obj.eNodeB_UE_LOS_direction_angles_signal(LTE_config,rx_height,attached_site_pos,UE_pos);
               [obj.azimuth_angle_of_departure_NLOS_once(:,:,pp),...
                   obj.azimuth_angle_of_arrival_NLOS_once(:,:,pp),...
                   obj.zenith_angle_of_departure_NLOS_once(:,:,pp),...
                   obj.zenith_angle_of_arrival_NLOS_once(:,:,pp)]       = obj.coupling_rays_NLOS(sigmas_NLOS,ZOD_parameters,obj.cl_power_NLOS,pp); 
               obj.RX_spherical_unit_vector_NLOS(:,:,:,pp)              = obj.receiver_spherical_unit_vector_NLOS(pp);
               obj.TX_spherical_unit_vector_NLOS(:,:,:,pp)              = obj.transmitter_spherical_unit_vector_NLOS(pp);
               obj.Doppler_speed_NLOS(:,:,:,pp)                         = obj.calculate_doppler_speed_NLOS(LTE_config, direction_of_movement, pp);
               initial_phases_NLOS                                      = obj.generate_initial_phases_NLOS;
               initial_phases_NLOS                                      = initial_phases_NLOS *pi/180;  % conversion to rad required
               obj.BS_boresight                                         = attached_eNodeB.boresight;
               % Antenna field
               for uu = 1:nRx
                   % Determine slant angle of UE antenna
                   if strcmp(LTE_config.UE_antenna_polarization,'ULA')
                       UE_slant_angle = LTE_config.UE_antenna_slant_angle;
                   else
                       UE_slant_angle = LTE_config.UE_antenna_slant_angle*mod(uu,2);
                   end
                   
                   for nn = 1:nTx
                       % Determine slant angle of antenna array column: 0...linearly polarized,  +45/-45...x-polarized
                       if strcmp(LTE_config.antenna.antenna_polarization,'XPOL')
                           slant_angle = LTE_config.antenna.slant_angle*(2*mod(nn,2)-1);
                       else
                           slant_angle = LTE_config.antenna.slant_angle;
                       end
                       [obj.UE_theta_AntennaField_global_NLOS(:,:,nn,uu,pp), obj.UE_phi_AntennaField_global_NLOS(:,:,nn,uu,pp)] = ...
                           obj.UE_antenna_field_pattern_global_NLOS(UE_slant_angle,direction_of_movement,pp);
                       [obj.theta_AntennaField_tx_global_NLOS(:,:,nn,pp), obj.phi_AntennaField_tx_global_NLOS(:,:,nn,pp)] = ...
                           obj.eNodeB_antenna_field_pattern_GCS_NLOS(LTE_config,attached_eNodeB,slant_angle,pp);
                       %part form generate_initial_phases_and_fields
                       f_matrix = (exp(1i.*initial_phases_NLOS(:,:,1)).* obj.theta_AntennaField_tx_global_NLOS(:,:,nn,pp) + sqrt(obj.xpr_NLOS.^(-1)).*exp(1i.*initial_phases_NLOS(:,:,2)).*obj.phi_AntennaField_tx_global_NLOS(:,:,nn,pp)) .* obj.UE_theta_AntennaField_global_NLOS(:,:,nn,uu,pp) +...
                           (sqrt(obj.xpr_NLOS.^(-1)).*exp(1i.*initial_phases_NLOS(:,:,3)).*obj.theta_AntennaField_tx_global_NLOS(:,:,nn,pp) + exp(1i.*initial_phases_NLOS(:,:,4)).*obj.phi_AntennaField_tx_global_NLOS(:,:,nn,pp)).* obj.UE_phi_AntennaField_global_NLOS(:,:,nn,uu,pp);
                       
                       field_and_phases_over_rays_and_clusters(:,:,nn,uu) = sqrt(obj.cl_power_NLOS_per_ray).*f_matrix;
                   end
               end
           end
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           %%%% USED FOR CALIBRATION ONLY
            % Save NLOS angles 
            function azimuth_of_departure_NLOS = azimuth_departure_angles_NLOS(obj)
               azimuth_of_departure_NLOS = obj.azimuth_angle_of_departure_NLOS_once(:,:,1);
            end
            
           function zenith_of_departure_NLOS = zenith_departure_angles_NLOS(obj)
               zenith_of_departure_NLOS  = obj.zenith_angle_of_departure_NLOS_once(:,:,1);
           end
            
           function azimuth_of_arrival_NLOS = azimuth_arrival_angles_NLOS(obj)
               azimuth_of_arrival_NLOS   = obj.azimuth_angle_of_arrival_NLOS_once(:,:,1);
           end
           function zenith_of_arrival_NLOS = zenith_arrival_angles_NLOS(obj)
               zenith_of_arrival_NLOS    = obj.zenith_angle_of_arrival_NLOS_once(:,:,1);
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
           function exp_part_without_Doppler_NLOS = calculate_exp_part_without_Doppler_NLOS(obj,LTE_config,rx_height,attached_site_pos,UE_pos,relative_antenna_position,relative_UE_antenna_position, pp)
               % Generates the exponential term in eq. 7.21 without the
               % Doppler component
               % A product of exponential terms at the transmitter and
               % receiver, where each exponential term is a product 
               % between receiver/transmitter spherical unit vector and 
               % the receiver/transmitter location vector
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               wavelength                         = 299792458./LTE_config.frequency;
               UE_location                        = obj.UE_position(rx_height,UE_pos,relative_UE_antenna_position);
               eNodeB_location                    = obj.eNodeB_position(LTE_config,attached_site_pos,relative_antenna_position);
               % Product between receiver spherical unit vector and the
               % receiver location vector
               size_vec = size(obj.RX_spherical_unit_vector_NLOS); 
               RX_vector_temp = reshape(obj.RX_spherical_unit_vector_NLOS(:,:,:,pp),3,size_vec(2)*size_vec(3));
               calculate_r_rx_and_location_vector = reshape(RX_vector_temp.'*UE_location.',size_vec(2),size_vec(3));
               % Product between transmitter spherical unit vector and the
               % transmitter location vector
               TX_vector_temp = reshape(obj.TX_spherical_unit_vector_NLOS(:,:,:,pp),3,size_vec(2)*size_vec(3));
               calculate_r_tx_and_location_vector = reshape(TX_vector_temp.'*eNodeB_location.',size_vec(2),size_vec(3));

               % Calculate exponential part of channel matrix
               exp_part_without_Doppler_NLOS   = exp(1i*2*pi/wavelength.*calculate_r_rx_and_location_vector).*exp(1i*2*pi/wavelength.*calculate_r_tx_and_location_vector);  
           end
         
           
           function Doppler_over_time_NLOS = calculate_Doppler_over_time_NLOS(obj, current_TTI, time_length, tSubframe, pp)
               % Generates the Doppler frequency component
               % The Doppler component is calculated for each time instant
               Doppler_over_time_NLOS = zeros(size(obj.Doppler_speed_NLOS,1),size(obj.Doppler_speed_NLOS,2), time_length);
               % Calculate Doopler effect for each time instant
               for tt = current_TTI:time_length
                   Doppler_over_time_NLOS(:,:,tt-current_TTI+1)           = exp(1i*2*pi.*obj.Doppler_speed_NLOS(:,:,pp) * tt*tSubframe);%check
                   % Calibration
                   % Doppler_over_time_NLOS(:,:,tt-current_TTI+1)           = exp(0*obj.Doppler_speed_NLOS(:,:,pp) * tt*tSubframe);%check
               end
           end
           
       
           
           function channel_coefficients_over_clusters_NLOS = calculate_channel_coefficient_NLOS(obj,LTE_config,rx_height,attached_site_pos,UE_pos,sigmas_NLOS,ZOD_parameters,attached_eNodeB, nTx, nRx, pp, current_TTI, time_length, tSubframe, direction_of_movement)
               % Generates the channel impulse response for each antenna
               % element pair between tramsitter and receiver
               % 
               % output: channel_coefficients_over_clusters_NLOS    ... a matrix of dimensions:
               %                                                       - 1st dimension is number of receive antenna ports
               %                                                       - 2nd dimension is number of transmit antenna ports
               %                                                       - 3rd dimension is 1 (no reason, can be removed)
               %                                                       - 4th dimension is number of time blocks
               %                                                       - 5th dimension is number of clusters
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               wavelength                   = 299792458./LTE_config.frequency;
               % Check polarization of transmit and receive antennas - for horizontal spacing of antenna elements.
               if strcmp(LTE_config.antenna.antenna_polarization,'XPOL') tx_spacing_factor = 2; else  tx_spacing_factor = 1; end
               if strcmp(LTE_config.UE_antenna_polarization,'XPOL')      rx_spacing_factor = 2; else  rx_spacing_factor = 1; end
               %Generate channel parameters for antenna fields and phases
               %over rays and clusters
               obj.field_and_phases_over_rays_and_clusters_NLOS = obj.calculate_initial_clusters_rays_phases_and_fields_NLOS(LTE_config, rx_height, attached_site_pos, UE_pos, sigmas_NLOS, ZOD_parameters, attached_eNodeB, nTx, nRx, direction_of_movement, pp);
               obj.Doppler_over_time_NLOS = obj.calculate_Doppler_over_time_NLOS(current_TTI, time_length, tSubframe, pp);
               for uu = 1:nRx
                   relative_UE_antenna_position  = [0, LTE_config.UE_antenna_element_horizontal_spacing*(ceil(uu/rx_spacing_factor)-1), 0];
                   element_channel_coefficients = zeros(obj.PerClusterRays,obj.NumClusters_NLOS, nTx, 1, LTE_config.nr_of_antenna_elements_in_each_column, time_length-current_TTI+1); % use nr of clusters instead numbers rays, clusters
                   element_weight               = zeros(nTx, LTE_config.nr_of_antenna_elements_in_each_column);
                   
                   for nn = 1:nTx
                       % Generate channel coefficients for each antenna
                       % element
                       for mm = 1:LTE_config.nr_of_antenna_elements_in_each_column
                           relative_antenna_position = [0, LTE_config.antenna_element_horizontal_spacing*(ceil(nn/tx_spacing_factor)-1), LTE_config.antenna_element_vertical_spacing*(mm - 1)];
                           obj.exponential_part_without_Doppler_NLOS = obj.calculate_exp_part_without_Doppler_NLOS(LTE_config,rx_height,attached_site_pos,UE_pos,relative_antenna_position,relative_UE_antenna_position, pp);
                           element_weight(nn,mm) = 1/sqrt(LTE_config.nr_of_antenna_elements_in_each_column)*exp(-1i*2*pi/wavelength*(mm-1)*LTE_config.antenna_element_vertical_spacing*cosd(LTE_config.electrical_downtilt));
                           % perform simulation from current TTI, for all
                           % subsequent TTIs
                           for tt=current_TTI:time_length
                               element_channel_coefficients(:,:,nn, 1 ,mm,tt-current_TTI+1)= obj.field_and_phases_over_rays_and_clusters_NLOS(:,:,nn,uu)...
                                   .* obj.exponential_part_without_Doppler_NLOS .* obj.Doppler_over_time_NLOS(:,:,tt-current_TTI+1);
                               
                               element_channel_coeff_over_clusters = sum(element_channel_coefficients,1); % sum ch-coeff. over rays within each cluster
                               rearranged_element_channel_coefficent= permute(element_channel_coeff_over_clusters,[3 5 6 4 2 1]);
                               
                               for kk=1:obj.NumClusters_NLOS
                                   portwise_channel_coefficients(:,1,uu,tt-current_TTI+1,kk) = sum((rearranged_element_channel_coefficent(:,:,tt-current_TTI+1,1,kk).*element_weight),2);
                               end
                               chanel_coeff_norm_NLOS(nn,1,uu,tt-current_TTI+1) = sum(abs(portwise_channel_coefficients(nn,1,uu,tt-current_TTI+1,:)).^2);
                           end
                       end
                   end 
               end
               channel_coefficients_over_clusters_NLOS = permute(portwise_channel_coefficients, [3 1 2 4 5]);
               obj.channel_coefficients_over_clusters_NLOS = permute(portwise_channel_coefficients, [3 1 2 4 5]); % channel matrix with dimensions [nRx, nTx, 1  nTTI, nCluster]
               obj.chan_norm_NLOS = permute(chanel_coeff_norm_NLOS, [3 1 2 4]);
           end
           
           % Save channel matrix
           function channel_matrix_NLOS = saved_channel_matrix_NLOS(obj)
               channel_matrix_NLOS = obj.channel_coefficients_over_clusters_NLOS;
           end

%% LOS
           % Calculate spherical unit vectors, Doppler and antenna field
           % This function is only called once per antenna array
           function [field_and_phases_over_rays_and_clusters_LOS, field_and_phases_single_LOS_ray] = calculate_initial_clusters_rays_phases_and_fields_LOS(obj,LTE_config, rx_height, attached_site_pos, UE_pos, sigmas_LOS, ZOD_parameters, attached_eNodeB, nTx, nRx, direction_of_movement, pp)
               % Generates only the field multiplication of the transmitter 
               % and receiver from eq 7.21 in TR 36.873
               % without the Doppler component and exponential terms
               % This function is only called once per antenna array
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               obj.xpr_LOS   = obj.generate_XPR_LOS;
               [obj.cl_power_LOS,obj.cl_power_LOS_general,obj.cl_power_LOS_per_ray, obj.ray_power_without_LOS_final] = obj.generate_cluster_power_LOS(sigmas_LOS);
               [obj.theta_arrival(pp),obj.theta_departure(pp),obj.phi_arrival(pp),obj.phi_departure(pp)]= obj.eNodeB_UE_LOS_direction_angles_signal(LTE_config,rx_height,attached_site_pos,UE_pos);
               [obj.azimuth_angle_of_departure_LOS_once(:,:,pp),...
                   obj.azimuth_angle_of_arrival_LOS_once(:,:,pp),...
                   obj.zenith_angle_of_departure_LOS_once(:,:,pp),...
                   obj.zenith_angle_of_arrival_LOS_once(:,:,pp)]      = obj.coupling_rays_LOS(sigmas_LOS,ZOD_parameters,obj.cl_power_LOS,pp);
               obj.Doppler_speed_LOS(:,:,pp)                          = obj.calculate_doppler_speed_LOS(LTE_config, direction_of_movement, pp);
               obj.Doppler_speed_LOS_direct_ray(:,:,pp)                   = obj.calculate_doppler_speed_LOS_direct_ray(LTE_config, direction_of_movement, pp);
               obj.RX_spherical_unit_vector_LOS(:,:,:,pp)             = obj.receiver_spherical_unit_vector_LOS(pp);
               obj.TX_spherical_unit_vector_LOS(:,:,:,pp)             = obj.transmitter_spherical_unit_vector_LOS(pp);
               obj.RX_spherical_unit_vector_LOS_direct_ray(:,:,:,pp)  = obj.receiver_spherical_unit_vector_LOS_direct_ray(pp);
               obj.TX_spherical_unit_vector_LOS_direct_ray(:,:,:,pp)  = obj.transmitter_spherical_unit_vector_LOS_direct_ray(pp);
               initial_phases_LOS            = obj.generate_initial_phases_LOS;
               initial_phases_LOS            = initial_phases_LOS *pi/180;  % conversion to rad required
               initial_phases_LOS_direct_ray = obj.generate_initial_phases_LOS_direct_ray;
               initial_phases_LOS_direct_ray = initial_phases_LOS_direct_ray *pi/180;
               % Antenna field
                obj.BS_boresight                                         = attached_eNodeB.boresight;
               for uu = 1:nRx
                   % Determine slant angle of UE antenna
                   if strcmp(LTE_config.UE_antenna_polarization,'ULA')
                       UE_slant_angle = LTE_config.UE_antenna_slant_angle;
                   else
                       UE_slant_angle = LTE_config.UE_antenna_slant_angle*mod(uu,2);
                   end
                   
                   for nn = 1:nTx
                       % Determine slant angle of antenna array column: 0...linearly polarized,  +45/-45...x-polarized
                       if strcmp(LTE_config.antenna.antenna_polarization,'XPOL')
                           slant_angle = LTE_config.antenna.slant_angle*(2*mod(nn,2)-1);
                       else
                           slant_angle = LTE_config.antenna.slant_angle;
                       end
                       [obj.UE_theta_AntennaField_global_LOS(:,:,nn,uu,pp), obj.UE_phi_AntennaField_global_LOS(:,:,nn,uu,pp)] = ...
                           obj.UE_antenna_field_pattern_global_LOS(UE_slant_angle,direction_of_movement,pp);
                       [obj.theta_AntennaField_tx_global_LOS(:,:,nn,pp), obj.phi_AntennaField_tx_global_LOS(:,:,nn,pp)] = ...
                           obj.antenna_field_pattern_GCS_LOS(LTE_config,attached_eNodeB,slant_angle,pp);
                       [obj.UE_theta_AntennaField_global_LOS_direct_ray(:,:,nn,uu,pp), obj.UE_phi_AntennaField_global_LOS_direct_ray(:,:,nn,uu,pp)] = ...
                           obj.UE_antenna_field_pattern_global_LOS_direct_ray(UE_slant_angle,direction_of_movement,pp);
                       [obj.theta_AntennaField_tx_global_LOS_direct_ray(:,:,nn,pp), obj.phi_AntennaField_tx_global_LOS_direct_ray(:,:,nn,pp)] = ...
                           obj.antenna_field_pattern_GCS_LOS_direct_ray(LTE_config,attached_eNodeB,slant_angle,pp);
                       
                       % Calculate UE antenna field pattern multiplied by random initial phasis matrix and eNodeB antenna field pattern
                       f_matrix = (exp(1i.*initial_phases_LOS(:,:,1)).* obj.theta_AntennaField_tx_global_LOS(:,:,nn,pp) + sqrt(obj.xpr_LOS.^(-1)).*exp(1i.*initial_phases_LOS(:,:,2)).*obj.phi_AntennaField_tx_global_LOS(:,:,nn,pp)) .* obj.UE_theta_AntennaField_global_LOS(:,:,nn,uu,pp) +...
                           (sqrt(obj.xpr_LOS.^(-1)).*exp(1i.*initial_phases_LOS(:,:,3)).*obj.theta_AntennaField_tx_global_LOS(:,:,nn,pp) + exp(1i.*initial_phases_LOS(:,:,4)).*obj.phi_AntennaField_tx_global_LOS(:,:,nn,pp)).* obj.UE_phi_AntennaField_global_LOS(:,:,nn,uu,pp);
                       %Generate channel_coefficients using Eq(7.21) [TR. 36.873]
                       field_and_phases_over_rays_and_clusters_LOS(:,:,nn,uu) = sqrt(obj.cl_power_LOS_per_ray).*f_matrix;
                       
                       % Calculate UE antenna field pattern multiplied by random initial phasis matrix and eNodeB antenna field pattern for a single LOS ray
                       f_matrix_LOS_ray = (exp(1i.*initial_phases_LOS_direct_ray(1)).*obj.theta_AntennaField_tx_global_LOS_direct_ray(:,:,nn,pp)).*obj.UE_theta_AntennaField_global_LOS_direct_ray(:,:,nn,uu,pp) +...
                           (-exp(1i.*initial_phases_LOS_direct_ray(2)).*obj.phi_AntennaField_tx_global_LOS_direct_ray(:,:,nn,pp)).*obj.UE_phi_AntennaField_global_LOS_direct_ray(:,:,nn,uu,pp);
                       field_and_phases_single_LOS_ray(:,:,nn,uu)        = f_matrix_LOS_ray;
                   end
               end
           end
           
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           %%%% USED FOR CALIBRATION ONLY
            % Save NLOS angles 
           % Save LOS angles
           function azimuth_of_departure_LOS = azimuth_departure_angles_LOS(obj)
               azimuth_of_departure_LOS = obj.azimuth_angle_of_departure_LOS_once(:,:,1);
           end
           
           function zenith_of_departure_LOS = zenith_departure_angles_LOS(obj)
               zenith_of_departure_LOS  = obj.zenith_angle_of_departure_LOS_once(:,:,1);
           end
           
           function azimuth_of_arrival_LOS = azimuth_arrival_angles_LOS(obj)
               azimuth_of_arrival_LOS   = obj.azimuth_angle_of_arrival_LOS_once(:,:,1);
           end
           
           function zenith_of_arrival_LOS = zenith_arrival_angles_LOS(obj)
               zenith_of_arrival_LOS    = obj.zenith_angle_of_arrival_LOS_once(:,:,1);
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
           function [exp_part_without_Doppler_LOS, exp_part_without_Doppler_direct_ray] = calculate_exp_part_without_Doppler_LOS(obj,LTE_config,rx_height,attached_site_pos,UE_pos,relative_antenna_position,relative_UE_antenna_position, pp)
               % Generates the exponential term in eq. 7.21 without the
               % Doppler component
               % A product of exponential terms at the transmitter and
               % receiver, where each exponential term is a product 
               % between receiver/transmitter spherical unit vector and 
               % the receiver/transmitter location vector
               %
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               wavelength                         = 299792458./LTE_config.frequency;
               UE_location                        = obj.UE_position(rx_height,UE_pos,relative_UE_antenna_position);
               eNodeB_location                    = obj.eNodeB_position(LTE_config,attached_site_pos,relative_antenna_position);
               % Product between receiver spherical unit vector and the
               % receiver location vector
               size_vec = size(obj.RX_spherical_unit_vector_LOS);
               RX_vector_temp = reshape(obj.RX_spherical_unit_vector_LOS(:,:,:,pp),3,size_vec(2)*size_vec(3));
               calculate_r_rx_and_location_vector = reshape(RX_vector_temp.'*UE_location.',size_vec(2),size_vec(3));
               % Product between transmitter spherical unit vector and the
               % transmitter location vector
               TX_vector_temp = reshape(obj.TX_spherical_unit_vector_LOS(:,:,:,pp),3,size_vec(2)*size_vec(3));
               calculate_r_tx_and_location_vector = reshape(TX_vector_temp.'*eNodeB_location.',size_vec(2),size_vec(3));
               % Calculate exponential part of channel matrix
               exp_part_without_Doppler_LOS = exp(1i*2*pi/wavelength.*calculate_r_rx_and_location_vector).*exp(1i*2*pi/wavelength.*calculate_r_tx_and_location_vector);
              
               % LOS ray
               % Calculate exponential part of channel matrix for LOS ray
               calculate_r_rx_and_location_vector_direct_ray = obj.RX_spherical_unit_vector_LOS_direct_ray(:,:,:,pp).'*UE_location.'; 
               calculate_r_tx_and_location_vector_direct_ray = obj.TX_spherical_unit_vector_LOS_direct_ray(:,:,:,pp).'*eNodeB_location.';
               % Calculate exponential part of channel matrix
               exp_part_without_Doppler_direct_ray = exp(1i*2*pi/wavelength.*calculate_r_rx_and_location_vector_direct_ray).*exp(1i*2*pi/wavelength.*calculate_r_tx_and_location_vector_direct_ray);    
           end
       
          
            function [Doppler_over_time_LOS, Doopler_over_time_LOS_direct_ray] = calculate_Doppler_over_time_LOS(obj, current_TTI, time_length, tSubframe, pp)
               % Generates the Doppler frequency component
               % The Doppler component is calculated for each time instant
               Doppler_over_time_LOS = zeros(size(obj.Doppler_speed_LOS,1),size(obj.Doppler_speed_LOS,2), time_length);
               Doopler_over_time_LOS_direct_ray = zeros(size(obj.Doppler_speed_LOS_direct_ray,1),size(obj.Doppler_speed_LOS_direct_ray,2),time_length);
               % Calculate Doopler effect for each time instant    
               for tt = current_TTI:time_length
                   Doppler_over_time_LOS(:,:,tt-current_TTI+1)           = exp(1i*2*pi.*obj.Doppler_speed_LOS(:,:,pp) * tt*tSubframe);
                   Doopler_over_time_LOS_direct_ray(:,:,tt-current_TTI+1) = exp(1i*2*pi.*obj.Doppler_speed_LOS_direct_ray(:,:,pp) * tt*tSubframe);
                   % Calibration
                   % Doppler_over_time_LOS(:,:,tt-current_TTI+1)           = exp(0*pi.*obj.Doppler_speed_LOS(:,:,pp) * tt*tSubframe);
                   % Doopler_over_time_LOS_direct_ray(:,:,tt-current_TTI+1) = exp(0*pi.*obj.Doppler_speed_LOS_direct_ray(:,:,pp) * tt*tSubframe);
               end
           end
           
           function channel_coefficients_LOS_final_matrix = calculate_channel_coefficient_LOS(obj,LTE_config,rx_height,attached_site_pos,UE_pos,sigmas_LOS,ZOD_parameters,attached_eNodeB, nTx, nRx,  pp, current_TTI, time_length, tSubframe, direction_of_movement)
               % Generates the channel impulse response for each antenna
               % element pair between tramsitter and receiver
               % 
               % output: channel_coefficients_over_clusters_LOS    ... a matrix of dimensions:
               %                                                       - 1st dimension is number of receive antenna ports
               %                                                       - 2nd dimension is number of transmit antenna ports
               %                                                       - 3rd dimension is 1 (no reason, can be removed)
               %                                                       - 4th dimension is number of time blocks
               %                                                       - 5th dimension is number of clusters
               % (c) Fjolla Ademaj, Martin Taranetz, ITC, 2016
               wavelength                   = 299792458./LTE_config.frequency;
               sigmas_KF_LOS = obj.KF_mu_LOS;
               % Check polarization of transmit and receive antennas - for horizontal spacing of antenna elements.
               if strcmp(LTE_config.antenna.antenna_polarization,'XPOL') tx_spacing_factor = 2; else  tx_spacing_factor = 1; end
               if strcmp(LTE_config.UE_antenna_polarization,'XPOL')      rx_spacing_factor = 2; else  rx_spacing_factor = 1; end
               %Generate channel parameters for antenna fields and phases
               %over rays and clusters
               [obj.field_and_phases_over_rays_and_clusters_LOS, obj.field_and_phases_single_LOS_ray] = obj.calculate_initial_clusters_rays_phases_and_fields_LOS(LTE_config, rx_height, attached_site_pos, UE_pos, sigmas_LOS, ZOD_parameters, attached_eNodeB, nTx, nRx, direction_of_movement, pp);
               [obj.Doppler_over_time_LOS, obj.Doopler_LOS_direct_ray] = obj.calculate_Doppler_over_time_LOS(current_TTI, time_length, tSubframe, pp);
               for uu = 1:nRx
                   relative_UE_antenna_position  = [0, LTE_config.UE_antenna_element_horizontal_spacing*(ceil(uu/rx_spacing_factor)-1), 0];
                   element_channel_coefficients = zeros(obj.PerClusterRays,obj.NumClusters_LOS,nTx,1, LTE_config.nr_of_antenna_elements_in_each_column, time_length-current_TTI+1);
                   element_channel_coefficients_LOS_ray = zeros(nTx,1, LTE_config.nr_of_antenna_elements_in_each_column);
                   element_weight               = zeros(nTx, LTE_config.nr_of_antenna_elements_in_each_column);
                   
                   for nn = 1:nTx
                       % Generate channel coefficients for each antenna
                       % element
                       for mm = 1:LTE_config.nr_of_antenna_elements_in_each_column
                           relative_antenna_position = [0, LTE_config.antenna_element_horizontal_spacing*(ceil(nn/tx_spacing_factor)-1), LTE_config.antenna_element_vertical_spacing*(mm - 1)];
                           [obj.exp_part_without_Doppler_LOS, obj.exp_part_without_Doppler_direct_ray] = obj.calculate_exp_part_without_Doppler_LOS(LTE_config,rx_height,attached_site_pos,UE_pos,relative_antenna_position,relative_UE_antenna_position, pp);
                           element_weight(nn,mm) = 1/sqrt(LTE_config.nr_of_antenna_elements_in_each_column)*exp(-1i*2*pi/wavelength*(mm-1)*LTE_config.antenna_element_vertical_spacing*cosd(LTE_config.electrical_downtilt));
                           % perform simulation from current TTI, for all
                           % subsequent TTIs
                           for tt=current_TTI:time_length
                               element_channel_coefficients(:,:,nn,1,mm,tt-current_TTI+1)= obj.field_and_phases_over_rays_and_clusters_LOS(:,:,nn,uu)...
                                   .* obj.exp_part_without_Doppler_LOS .* obj.Doppler_over_time_LOS(:,:,tt-current_TTI+1);
                               
                               element_channel_coefficients_LOS_ray(nn,1,mm,tt-current_TTI+1)= obj.field_and_phases_single_LOS_ray(:,:,nn,uu)...
                                   .* obj.exp_part_without_Doppler_direct_ray .* obj.Doopler_LOS_direct_ray(tt-current_TTI+1);
                               
                               element_channel_coeff_over_clusters = sum(element_channel_coefficients,1); % sum ch-coeff. over rays within each cluster
                               rearranged_element_channel_coefficent= permute(element_channel_coeff_over_clusters,[3 5 6 4 2 1]);
                               rearranges_element_channel_coefficients_LOS_ray = permute(element_channel_coefficients_LOS_ray,[1 3 2 4]);
                               
                               for kk=1:obj.NumClusters_LOS
                                   portwise_channel_coefficients(:,1,uu,tt-current_TTI+1,kk) = sum((rearranged_element_channel_coefficent(:,:,tt-current_TTI+1,1,kk).*element_weight),2);
                                   portwise_channel_coefficients_LOS_ray(:,1,uu,tt-current_TTI+1,kk)= sum((rearranges_element_channel_coefficients_LOS_ray(:,:,:,tt-current_TTI+1).*element_weight),2);
                               end
                               chanel_coeff_norm_LOS(nn,1,uu,tt-current_TTI+1) = sum(abs(portwise_channel_coefficients(nn,1,uu,tt-current_TTI+1,:)).^2);
                               chanel_coeff_norm_LOS_ray(nn,1,uu,tt-current_TTI+1) = sum(abs(portwise_channel_coefficients_LOS_ray(nn,1,uu,tt-current_TTI+1,:)).^2);
                           end
                       end
                   end  
               end
               % Determine the single LOS ray channel coefficients
               % to be added to other multipath components
               dirac_vector = zeros(size(portwise_channel_coefficients_LOS_ray,1),size(portwise_channel_coefficients_LOS_ray,2),size(portwise_channel_coefficients_LOS_ray,3),size(portwise_channel_coefficients_LOS_ray,4),obj.NumClusters_LOS);
               dirac_vector(:,:,:,:,1) = 1;
               dirac_vector = permute(dirac_vector, [3 1 2 4 5]);

                channel_coefficients_over_clusters_LOS = permute(portwise_channel_coefficients, [3 1 2 4 5]); % channel matrix with dimensions [nRx, nTx, 1  nTTI, nCluster] 
                channel_coefficients_over_clusters_LOS_ray = permute(portwise_channel_coefficients_LOS_ray, [3 1 2 4 5]); % channel matrix with dimensions [nRx, nTx, 1  nTTI, nCluster]
                portwize_channel_coeff_LOS_ray_transformed = dirac_vector.*sqrt(10^(sigmas_KF_LOS/10)/(10^(sigmas_KF_LOS/10) +1)).* channel_coefficients_over_clusters_LOS_ray;
               obj.chan_norm_LOS = permute(chanel_coeff_norm_LOS, [3 1 2 4]);
               obj.chan_norm_LOS_ray = permute(chanel_coeff_norm_LOS_ray, [3 1 2 4]);

               channel_coefficients_LOS_final_matrix = sqrt(1/(10^(sigmas_KF_LOS/10) +1)).*channel_coefficients_over_clusters_LOS + portwize_channel_coeff_LOS_ray_transformed;
               obj.channel_coefficients_over_clusters_LOS = channel_coefficients_LOS_final_matrix;
           end
           
           % Save channel matrix
           function channel_matrix_LOS = saved_channel_matrix_LOS(obj)
               channel_matrix_LOS = obj.channel_coefficients_over_clusters_LOS;
           end
           
           %% Sample channel coefficients based on delay taps chenerated per each cluster
           function channel_coefficients_powers_LOS_norm = sampled_channel_LOS(obj, nRx, nTx, time_length, current_TTI, channel_coefficient)
               % Samples the channel impulse response based on the delay
               % taps
               % The sampled channel impulse response is normalized 
               
               delay_tap = round(obj.delay_LOS*obj.fs);
               obj.delay_tap_size_LOS = max(delay_tap+1);
               tap_position = [1, find(diff(delay_tap)) + 1];
               % Sum up all of the taps that merge (sum power!): nearest neighbor interpolation
               channel_coefficients_powers_LOS = zeros(nRx,nTx,1,time_length-current_TTI+1,obj.delay_tap_size_LOS);
               norm_sampled_channel = zeros(nRx,nTx,1,time_length-current_TTI+1);
               for tap_idx = unique(delay_tap)
                   taps_to_sum = delay_tap==tap_idx;
                   channel_coefficients_powers_LOS(:,:,:,:,tap_idx+1) = sum(channel_coefficient(:,:,:,:,taps_to_sum),5); 
               end
               for uu= 1:nRx
                   for nn=1:nTx
                       for tt=1:time_length-current_TTI+1
                           norm_sampled_channel(uu,nn,1,tt) = sum(abs(channel_coefficients_powers_LOS(uu,nn,1,tt,:)).^2);
                       end
                   end
               end

               for ss = 1: size(channel_coefficients_powers_LOS,5)
                   channel_coefficients_powers_LOS_norm(:,:,:,:,ss)= channel_coefficients_powers_LOS(:,:,:,:,ss).*sqrt(obj.chan_norm_LOS./norm_sampled_channel);
               end
               
               obj.int_sampled_channel_LOS= channel_coefficients_powers_LOS_norm;
           end
           
           function delay_tap_dimension = saved_delay_tap_dimension_LOS(obj)
               delay_tap_dimension = obj.delay_tap_size_LOS;
           end
           
           function sampled_channel_NLOS = saved_sampled_channel_LOS(obj)
               sampled_channel_NLOS = obj.int_sampled_channel_LOS;
           end
                     
           function channel_coefficients_powers_NLOS_norm = sampled_channel_NLOS(obj, nRx, nTx, time_length, current_TTI, channel_coefficient)
               % Samples the channel impulse response based on the delay
               % taps
               % The sampled channel impulse response is normalized 
               delay_tap = round(obj.delay_NLOS*obj.fs);
               obj.delay_tap_size_NLOS = max(delay_tap+1);
               tap_position = [1, find(diff(delay_tap)) + 1];
               % Sum up all of the taps that merge (sum power!): nearest neighbor interpolation
               channel_coefficients_powers_NLOS = zeros(nRx,nTx,1,time_length-current_TTI+1, obj.delay_tap_size_NLOS);
               norm_sampled_channel = zeros(nRx,nTx,1,time_length-current_TTI+1);
               for tap_idx = unique(delay_tap)
                   taps_to_sum = delay_tap==tap_idx;
                   channel_coefficients_powers_NLOS(:,:,:,:,tap_idx+1) = sum(channel_coefficient(:,:,:,:,taps_to_sum),5); % 
               end
               for uu= 1:nRx
                   for nn=1:nTx
                       for tt=1:time_length-current_TTI+1
                           norm_sampled_channel(uu,nn,1,tt) = sum(abs(channel_coefficients_powers_NLOS(uu,nn,1,tt,:)).^2);
                       end
                   end
               end
               for ss = 1: size(channel_coefficients_powers_NLOS,5)
                   channel_coefficients_powers_NLOS_norm(:,:,:,:,ss)= channel_coefficients_powers_NLOS(:,:,:,:,ss).*sqrt(obj.chan_norm_NLOS./norm_sampled_channel);
               end
               
               obj.int_sampled_channel_NLOS = channel_coefficients_powers_NLOS_norm;
           end
           
           function delay_tap_dimension = saved_delay_tap_dimension_NLOS(obj)  
                delay_tap_dimension = obj.delay_tap_size_NLOS;
           end
            
          function sampled_channel_NLOS = saved_sampled_channel_NLOS(obj)  
                sampled_channel_NLOS = obj.int_sampled_channel_NLOS;
          end
          
          
          % OTOI
           function channel_coefficients_powers_OTOI_norm = sampled_channel_OTOI(obj, nRx, nTx, time_length, current_TTI, channel_coefficient)
               % Samples the channel impulse response based on the delay
               % taps
               % The sampled channel impulse response is normalized 
               delay_tap = round(obj.delay_OTOI*obj.fs);
               obj.delay_tap_size_OTOI = max(delay_tap+1);
               tap_position = [1, find(diff(delay_tap)) + 1];
               % Sum up all of the taps that merge (sum power!): nearest neighbor interpolation
               channel_coefficients_powers_OTOI = zeros(nRx,nTx,1,time_length-current_TTI+1, obj.delay_tap_size_OTOI);
               norm_sampled_channel = zeros(nRx,nTx,1,time_length-current_TTI+1);
               for tap_idx = unique(delay_tap)
                   taps_to_sum = delay_tap==tap_idx;
                   channel_coefficients_powers_OTOI(:,:,:,:,tap_idx+1) = sum(channel_coefficient(:,:,:,:,taps_to_sum),5); 
               end
               for uu= 1:nRx
                   for nn=1:nTx
                       for tt=1:time_length-current_TTI+1
                           norm_sampled_channel(uu,nn,1,tt) = sum(abs(channel_coefficients_powers_OTOI(uu,nn,1,tt,:)).^2);
                       end
                   end
               end
               for ss = 1: size(channel_coefficients_powers_OTOI,5)
                   channel_coefficients_powers_OTOI_norm(:,:,:,:,ss)= channel_coefficients_powers_OTOI(:,:,:,:,ss).*sqrt(obj.chan_norm_OTOI./norm_sampled_channel);
               end
               
               obj.int_sampled_channel_OTOI = channel_coefficients_powers_OTOI_norm;
           end
           
           function delay_tap_dimension = saved_delay_tap_dimension_OTOI(obj)  
                delay_tap_dimension = obj.delay_tap_size_OTOI;
           end
            
          function sampled_channel_OTOI = saved_sampled_channel_OTOI(obj)  
                sampled_channel_OTOI = obj.int_sampled_channel_OTOI;
          end


           %% FFT over the sampled channel coefficients
           
           function H_fft_RB = get_RB_trace(obj, channel)
               % Returns back a frequency channel trace jumping each FFT_sampling_interval subcarriers
               H_fft_large = fft(channel,obj.Nfft,5);
               % Eliminate guardband
               H_fft       = H_fft_large(:,:,:,:,[obj.Nfft-obj.Ntot/2+1:obj.Nfft 2:obj.Ntot/2+1]);
               % Do not return the channel for all subcarriers, but just a subset of it
               H_fft_RB    = H_fft(:,:,:,:,1:obj.FFT_sampling_interval:end);
           end

           
            function H_fft_RB = get_RB_trace_per_subcarrier(obj, channel)
               % Returns back a frequency channel trace jumping each FFT_sampling_interval subcarriers
               H_fft_large = fft(channel,obj.Nfft,5);
               % Eliminate guardband
               H_fft       = H_fft_large(:,:,:,:,[obj.Nfft-obj.Ntot/2+1:obj.Nfft 2:obj.Ntot/2+1]);
               H_fft_RB    = H_fft(:,:,:,:,1:1:end); % calculate the channel per each subcarrier---Debugging/Calibration 
           end
    end
end 

      
