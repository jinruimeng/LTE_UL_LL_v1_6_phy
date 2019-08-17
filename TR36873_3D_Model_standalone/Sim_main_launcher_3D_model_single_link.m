%% 3GPP TR 36.973: 3D Channel Model - Main launcher file
% This is the launcher file that generates the 3D model channel trace based on 3GPP TR 36.873. 

% (c) Fjolla Ademaj, Martin Taranetz, ITC 2016

clc;
clear all;
close all;
config.parallel_network = false;

%% Input parameters 
% * In this part are given the following settings
% * Basic settings
% * User and eNodeB settings
% * Network geometry
% * UE antenna parameters
% * eNodeB antenna parameters
% * Load specific channel model parameters from Table 7.3-6 and system parameters for the FFT of the channel impulse response


% Basic settings
config.frequency            = 2e9;
config.wavelength           = 299792458/config.frequency;
config.debug_level          = 1;
config.bandwidth            = 10e6; 
config.BF_length            = 1e-3/14;       % Denotes the Block Fading length
config.show_network_plots   = true;
config.CP_length            = 'normal';

% User and eNodeB settings
config.UE_per_eNodeB        = 1;
config.UE_speed             = 50/3.6; % m/s
config.eNodeB_tx_power      = 40;
config.eNodeB_nTX           = 4;
config.UE_nRX               = 2;
config.channel_model.type   = '3D_UMa_fading';  % '3D_UMa_fading'; '3D_UMi_fading'
config.min_floor_number       = 4; % number of floors in buildings
config.max_floor_number       = 8; % number of floors in buildings
config.indoor_UE_fraction     = 0.8; % fraction of indoor users

% Network geometry
config.tx_height                                       = 25;     %in [m];   3D_UMa_fading: 25[m]; 3D_UMi_fading: 10[m]; 
config.min_UE_eNodeB_distance                          = 35;     %in [m];   3D_UMa_fading: 35[m]; 3D_UMi_fading: 10[m]; 
config.inter_eNodeB_distance                           = 500;    %in [m];   3D_UMa_fading: 500[m]; 3D_UMi_fading: 200[m]; 
config.nr_eNodeB_rings                                 = 0;
config.nr_sectors                                      = 1;
config.boresight_direction                             = 45;     % in [deg];  The boresight direction of the eNodeB; This is a desired parameter, can be changed

% UE antenna parameters
config.UE_antenna_polarization                       = 'ULA';
config.UE_antenna_slant_angle                        = 0;     % if ULA --> slant angle is 0deg; if 'XPOL' --> slant_angle is 90deg
config.UE_antenna_element_horizontal_spacing         = 0.5 * config.wavelength;

% eNodeB antenna parameters
config.antenna.antenna_gain_pattern                  = 'TR36.873 3D antenna';
config.antenna.antenna_polarization                  = 'ULA'; % 'XPOL' or 'COPOL'.
config.antenna.slant_angle                           = 0; % If 'COPOL' --> slant_angle is 0deg, 'XPOL' --> slant_angle is 45deg
config.antenna_element_vertical_spacing              = 0.5 * config.wavelength;     % 0.8*wavelength % 
config.antenna_element_horizontal_spacing            = 0.5 * config.wavelength;   % 0.5*wavelength can also be 0.8*wavelength
config.nr_of_antenna_elements_in_each_column         = 10;
config.electrical_downtilt                           = 90;
config.mechanical_downtilt                           = 0;
config.mechanical_slant                              = 0;
config.antenna.max_antenna_gain                      = 8; % As defined in TR 36.873 (2000 MHz)

% Load specific channel model parameters from Table 7.3-6 and system
% parameters for the FFT of the channel impulse response
config = load_specific_params_link_level(config);

config.simulation_time_blocks                         = 100; % Design parameter: the simulation length as multiple of Block Fading length (config.BF_length)
config.results_folder                                 = './results';
config.results_file                                   = 'auto'; %NOTE: 'auto' assigns a filename automatically


%% Generate network layout 
% Generates the eNodeB and UE position
%
% # Generates eNodeBs properties
%
% * |eNodeBs.pos|          - The eNodeB position [x, y]
% * |eNodeBs.tx_height|    - The eNodeB height in [m]
% * |eNodeBs.nTX|          - The number of antenna ports at eNodeB
% * |eNodeBs.antenna_type| - Antenna model type
% * |eNodeBs.boresight|    - The boresight direction of antenna pointing in [deg]
%
% # Generate UE properties 
%
% * |UEs.pos|                  - The UE position in [x, y]
% * |UEs.rx_height|            - The UE height in [m]
% * |UEs.direction|            - The Ue direction in [deg]     
% * |UEs.nRX|                  - The number of antenna ports at UE
% * |UEs.is_LOS|               - Boolean operator that indicates the LOS/NLOS condition of the UE; 1->LOS; 0->NLOS
% * |UEs.is_indoor|            - Boolean operator that indicates the indoor/outdoor condition of the UE; 1->indoors; 0->outdoors
% * |UEs.dist_indoor|          - The indoor distance in [m] from the UE postion to the outer wall of the building where the Ue is located
% * |UEs.H_0_channel_trace|    - The channel impulse response from the channel model
% * |UEs.sampled_channel_H_0|  - The sampled channel impulse response
% * |UEs.H_0_final|            - The channel transfer function 

fprintf('Generating eNodeB and UE position');
eNodeB_site_positions = [0 0];
UE_positions = randi([config.min_UE_eNodeB_distance,config.inter_eNodeB_distance],1,2); %

%
% Generate eNodeBs properties
fprintf('Generating eNodeB properites');
eNodeB_sites = network_elements.eNodeB_site;
eNodeBs = network_elements.eNodeB;
eNodeB_sites.id = 1;
eNodeB_sites.pos = eNodeB_site_positions;
eNodeB_sites.site_type = 'macro';
eNodeBs = network_elements.eNodeB;
eNodeBs.id = 1;
eNodeBs.eNodeB_id = 1;
eNodeBs.parent_eNodeB = eNodeB_sites;
eNodeBs.pos = eNodeB_site_positions;
eNodeBs.max_power = config.eNodeB_tx_power;
eNodeBs.nTX = config.eNodeB_nTX;
eNodeBs.tx_height = config.tx_height;
eNodeBs.boresight = config.boresight_direction;
eNodeBs.antenna_type  = config.antenna.antenna_gain_pattern;
eNodeBs = antennas.antenna.attach_antenna_to_eNodeB(eNodeBs,config);

% Generate UE properties
fprintf('Generating UE properties');
UEs                     = network_elements.UE;
UEs.id                  = 1;
UEs.pos                 = UE_positions;
UEs.nRX                 = config.UE_nRX;
UEs.direction           = floor(random('unif',0,359));
UEs.H_0_channel_trace   = [];  % denotes the channel impulse response from the channel model
UEs.sampled_channel_H_0 = [];  % denotes the sampled channel impulse response
UEs.H_0_final           = [];  % denotes the channel transfer function after the FFT
eNodeBs.attachUser(UEs);
UEs.is_indoor           = binornd(1,config.indoor_UE_fraction); % determines whether the UE is indoors or outdoors
UEs.dist_indoor         = UEs.is_indoor.*(25*rand); % generates the indoor distance from the UE position to the outer wall of the building where the UE located
indoor_UE_height            = 3*(randi(randi([config.min_floor_number,config.max_floor_number])) - 1) + 1.5; % the height for an indoor UE as specified in the standard
UEs.rx_height           = UEs.is_indoor.*indoor_UE_height + ~UEs.is_indoor.*1.5; % determines the UE height: if outdoor UE the default value is 1.5 as given in the standard
UE_distance_to_site         = abs(sqrt(sum((UEs.attached_site.pos - UEs.pos).^2)) - UEs.dist_indoor);
UEs.is_LOS              = macroscopic_pathloss_models.generate_LOS_probability(UE_distance_to_site, UEs.rx_height, config.channel_model.type);
UEs.channels = channel_gain_wrappers.TR36873_Fading_3D_Channel(config);

% initialize clock
networkClock = network_elements.clock(config.BF_length);
eNodeBs.clock = networkClock;

%% Generate Correlated Large Scale Parameters
% Generates the correlated large scale parameters 
%
% |UEs.all_large_scale_params|     - Denotes a vector of large scale parameters in the following order [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
%
% * SF - Shadow fading
% * K-factor - Ricean K-factor
% * DS - Delay spread
% * ASD - Azimuth spread of departure
% * ASA - Azimuth spread of arrival
% * ZSD - Zenith spread of departure
% * ZSA - Zenith spread of arrival

fprintf('Calculating large scale parameters for eNodeb %d \n', eNodeBs.eNodeB_id);

if eNodeBs.attached_UEs_vector.is_indoor
    UE_distance_to_site = abs(sqrt(sum((eNodeBs.parent_eNodeB.pos - eNodeBs.attached_UEs_vector.pos).^2)) - eNodeBs.attached_UEs_vector.dist_indoor);
    if eNodeBs.attached_UEs_vector.is_LOS
        [ZS_D_mu_OTOI_LOS,ZS_D_sigma_OTOI_LOS,ZOD_mu_offset_OTOI_LOS] = eNodeBs.generate_ZSD_ZoD_offset_parameters_OTOI_LOS(config, eNodeBs.attached_UEs_vector.rx_height, UE_distance_to_site);
        eNodeBs.attached_UEs_vector.all_ZOD_params =  [ZS_D_mu_OTOI_LOS,ZS_D_sigma_OTOI_LOS,ZOD_mu_offset_OTOI_LOS];
    else
        [ZS_D_mu_OTOI_NLOS,ZS_D_sigma_OTOI_NLOS,ZOD_mu_offset_OTOI_NLOS] = eNodeBs.generate_ZSD_ZoD_offset_parameters_OTOI_NLOS(config, eNodeBs.attached_UEs_vector.rx_height, UE_distance_to_site);
        eNodeBs.attached_UEs_vector.all_ZOD_params =  [ZS_D_mu_OTOI_NLOS,ZS_D_sigma_OTOI_NLOS,ZOD_mu_offset_OTOI_NLOS];
    end
    [sigmas_SF_OTOI, non_value_param, sigmas_DS_OTOI, sigmas_ASD_OTOI, sigmas_ASA_OTOI, sigmas_ZSD_OTOI, sigmas_ZSA_OTOI] = eNodeBs.sigmas_OTOI(config, randn(7,1),eNodeBs.attached_UEs_vector.all_ZOD_params);
    eNodeBs.attached_UEs_vector.all_large_scale_params = [sigmas_SF_OTOI, non_value_param, sigmas_DS_OTOI, sigmas_ASD_OTOI, sigmas_ASA_OTOI, sigmas_ZSD_OTOI, sigmas_ZSA_OTOI];
else
    UE_distance_to_site = sqrt(sum((eNodeBs.parent_eNodeB.pos - eNodeBs.attached_UEs_vector.pos).^2));
    if eNodeBs.attached_UEs_vector.is_LOS
        [ZS_D_mu_LOS,ZS_D_sigma_LOS,ZOD_mu_offset_LOS] = eNodeBs.generate_ZSD_ZoD_offset_parameters_LOS(config, eNodeBs.attached_UEs_vector.rx_height, UE_distance_to_site);
        eNodeBs.attached_UEs_vector.all_ZOD_params = [ZS_D_mu_LOS,ZS_D_sigma_LOS,ZOD_mu_offset_LOS];
        [sigmas_SF_LOS, sigmas_KF_LOS, sigmas_DS_LOS, sigmas_ASD_LOS, sigmas_ASA_LOS, sigmas_ZSD_LOS, sigmas_ZSA_LOS] = eNodeBs.sigmas_LOS(config, randn(7,1), eNodeBs.attached_UEs_vector.all_ZOD_params);
        eNodeBs.attached_UEs_vector.all_large_scale_params = [sigmas_SF_LOS, sigmas_KF_LOS, sigmas_DS_LOS, sigmas_ASD_LOS, sigmas_ASA_LOS, sigmas_ZSD_LOS, sigmas_ZSA_LOS];
        
    else
        [ZS_D_mu_NLOS,ZS_D_sigma_NLOS,ZOD_mu_offset_NLOS] = eNodeBs.generate_ZSD_ZoD_offset_parameters_NLOS(config, eNodeBs.attached_UEs_vector.rx_height, UE_distance_to_site);
        eNodeBs.attached_UEs_vector.all_ZOD_params =  [ZS_D_mu_NLOS,ZS_D_sigma_NLOS,ZOD_mu_offset_NLOS];
        [sigmas_SF_NLOS, non_value_param,  sigmas_DS_NLOS, sigmas_ASD_NLOS, sigmas_ASA_NLOS, sigmas_ZSD_NLOS, sigmas_ZSA_NLOS] = eNodeBs.sigmas_NLOS(config, randn(7,1),eNodeBs.attached_UEs_vector.all_ZOD_params);
        eNodeBs.attached_UEs_vector.all_large_scale_params = [sigmas_SF_NLOS, non_value_param, sigmas_DS_NLOS, sigmas_ASD_NLOS, sigmas_ASA_NLOS, sigmas_ZSD_NLOS, sigmas_ZSA_NLOS];
    end
end



%% Generate channel impulse response
% Generates the channel impulse response
% * UEs.H_0_channel_trace   - Denotes the channel impulse response 
% * UEs.sampled_channel_H_0   - Denotes the sampled channel impulse response 

fprintf('Generating the channel impulse response for UE %d \n', eNodeBs.attached_UEs_vector.id);

if eNodeBs.attached_UEs_vector.is_indoor
    eNodeBs.attached_UEs_vector.H_0_channel_trace = zeros (eNodeBs.attached_UEs_vector.nRX,  eNodeBs.nTX, 1, config.simulation_time_blocks, config.NumClusters_OTOI);
    eNodeBs.attached_UEs_vector.channels.calculate_channel_coefficient_OTOI(config, eNodeBs.attached_UEs_vector.rx_height, eNodeBs.pos, eNodeBs.attached_UEs_vector.pos, eNodeBs.attached_UEs_vector.all_large_scale_params, eNodeBs.attached_UEs_vector.all_ZOD_params, eNodeBs, eNodeBs.nTX, eNodeBs.attached_UEs_vector.nRX, eNodeBs.attached_UEs_vector.is_indoor, 1, networkClock.current_block, config.simulation_time_blocks, networkClock.BF_time, eNodeBs.attached_UEs_vector.direction, eNodeBs.attached_UEs_vector.dist_indoor);
    eNodeBs.attached_UEs_vector.H_0_channel_trace = eNodeBs.attached_UEs_vector.channels.saved_channel_matrix_OTOI;
    eNodeBs.attached_UEs_vector.sampled_channel_H_0 = eNodeBs.attached_UEs_vector.channels.sampled_channel_OTOI(eNodeBs.attached_UEs_vector.nRX, eNodeBs.nTX, config.simulation_time_blocks, networkClock.current_block, eNodeBs.attached_UEs_vector.H_0_channel_trace);
    
else
    if eNodeBs.attached_UEs_vector.is_LOS
        eNodeBs.attached_UEs_vector.H_0_channel_trace = zeros(eNodeBs.attached_UEs_vector.nRX, eNodeBs.nTX, 1, config.simulation_time_blocks, config.NumClusters_LOS);
        eNodeBs.attached_UEs_vector.channels.calculate_channel_coefficient_LOS(config, eNodeBs.attached_UEs_vector.rx_height, eNodeBs.pos, eNodeBs.attached_UEs_vector.pos, eNodeBs.attached_UEs_vector.all_large_scale_params, eNodeBs.attached_UEs_vector.all_ZOD_params, eNodeBs, eNodeBs.nTX, eNodeBs.attached_UEs_vector.nRX, 1, networkClock.current_block, config.simulation_time_blocks, networkClock.BF_time, eNodeBs.attached_UEs_vector.direction);
        eNodeBs.attached_UEs_vector.H_0_channel_trace = eNodeBs.attached_UEs_vector.channels.saved_channel_matrix_LOS;
        eNodeBs.attached_UEs_vector.sampled_channel_H_0 = eNodeBs.attached_UEs_vector.channels.sampled_channel_LOS(eNodeBs.attached_UEs_vector.nRX, eNodeBs.nTX, config.simulation_time_blocks, networkClock.current_block, eNodeBs.attached_UEs_vector.H_0_channel_trace);
    else
        eNodeBs.attached_UEs_vector.H_0_channel_trace = zeros(eNodeBs.attached_UEs_vector.nRX, eNodeBs.nTX, 1, config.simulation_time_blocks, config.NumClusters_NLOS);
        eNodeBs.attached_UEs_vector.channels.calculate_channel_coefficient_NLOS(config, eNodeBs.attached_UEs_vector.rx_height, eNodeBs.pos, eNodeBs.attached_UEs_vector.pos, eNodeBs.attached_UEs_vector.all_large_scale_params, eNodeBs.attached_UEs_vector.all_ZOD_params, eNodeBs, eNodeBs.nTX, eNodeBs.attached_UEs_vector.nRX, 1, networkClock.current_block, config.simulation_time_blocks, networkClock.BF_time, eNodeBs.attached_UEs_vector.direction);
        eNodeBs.attached_UEs_vector.H_0_channel_trace = eNodeBs.attached_UEs_vector.channels.saved_channel_matrix_NLOS;
        eNodeBs.attached_UEs_vector.sampled_channel_H_0 = eNodeBs.attached_UEs_vector.channels.sampled_channel_NLOS(eNodeBs.attached_UEs_vector.nRX, eNodeBs.nTX, config.simulation_time_blocks, networkClock.current_block, eNodeBs.attached_UEs_vector.H_0_channel_trace);
    end
end

%% Save the results
%check if the folder exists
if ~exist(config.results_folder,'dir')
    mkdir(config.results_folder);
end

the_date = clock;
date_time_string = sprintf('%04d%02d%02d_%02d%02d%02d',...
    the_date(1),...                     %Date: year
    the_date(2),...                     %Date: month
    the_date(3),...                     %Date: day
    the_date(4),...                     %Date: hour
    the_date(5),...                     %Date: minutes
    floor(the_date(6)));

if exist('labindex','builtin')   %does not exist in older Matlab versions
    lab_ind = labindex;
else
    lab_ind = 1;
end

if config.parallel_network
    t = getCurrentTask();
    if ~isempty(t)
        par_ID = num2str(t.ID);
    else
        par_ID = '';
    end
else
    par_ID = '';
end
if config.frequency/1e9>=1
    this_freq  = sprintf('%3.2fGHz',config.frequency/1e9);
else
    this_freq = sprintf('%3.0fMHz',config.frequency/1e6);
end

this_bw = sprintf('%gfMHz',config.bandwidth/1e6);

if strcmp(config.results_file,'auto')
    results_filename = sprintf('%s_freq_%s_bw_%s_%.1fKmph_%s_lab%02.0f_%s.mat',...
        this_freq,...                               % Frequency
        this_bw,...                                 % System Bandwidth
        config.channel_model.type,...         % Channel model userd
        config.UE_speed*3.6,...               % Channel speed (Km/h)
        date_time_string,...                  % Date string
        lab_ind,...                                 % Lab index (in order to avoid multiple labs writing to the same file)
        par_ID);                                    % ID of the current worker (for parallel simulations only)
end
save('-v7.3',fullfile(config.results_folder,results_filename),'config','UEs','eNodeBs','eNodeB_sites');

