%% 3GPP TR 36.973: 3D Channel Model - Main launcher file
% - Generates a hexagonal grid of theree-sectorized eNodeb sites 
% - Generates the channel impulse response and the channel transfer function for all UEs served by each eNodeB
%   sector: desired channel and interfering channels using the 3GPP 3D channel model

% (c) Fjolla Ademaj, Martin Taranetz, ITC 2016

clc;
clear all;
close all;
config.parallel_network = false;

%% Input parameters 
% * In this part are given the following settings
% * Basic settings
% * User and eNodeB settings
% * Small- and large scale fading
% * Network geometry
% * UE antenna parameters
% * eNodeB antenna parameters
% * Load specific channel model parameters from Table 7.3-6 and system parameters for the FFT of the channel impulse response


% Basic settings
config.frequency            = 2e9;
config.wavelength           = 299792458/config.frequency;
config.debug_level          = 1;
config.bandwidth            = 10e6; 
config.BF_length            = 1e-3;   % Denotes the Block Fading length: e.g. For system level the TTI length of 1e-3; For link level the symbol length 1e-3/14
config.show_network_plots   = true;
config.CP_length            = 'normal';

% User and eNodeB settings
config.UE_per_eNodeB             = 3;
config.UE_speed                  = 50/3.6; % m/s
config.eNodeB_tx_power           = 40;
config.eNodeB_nTX                = 4;
config.UE_nRX                    = 2;
config.UE.thermal_noise_density  = -174;
config.UE.receiver_noise_figure  = 9;
config.min_floor_number          = 4; % number of floors in buildings
config.max_floor_number          = 8; % number of floors in buildings
config.indoor_UE_fraction        = 0.8; % fraction of indoor users

% Small- and large scale fading
config.channel_model.type                              = '3D_UMa_fading';  % '3D_UMa_fading'; '3D_UMi_fading'
config.macroscopic_pathloss_model                      = 'TR 36.873';
config.macroscopic_pathloss_model_settings.environment = '3D_UMa';         % '3D_UMa'; '3D_UMi'
config.shadow_fading_type                              = 'claussen_3D_UMa';  % 'none'; 'claussen_3D_UMa'; 'claussen_3D_UMi'
config.macroscopic_pathloss_is_model                   = true;
config.decouple_site_shadow_fading_maps                = false;
config.claussen_stream                                 = RandStream('mt19937ar','Seed',0);
config.deactivate_claussen_spatial_correlation         = false;
config.network_source                                  = 'generated';
config.minimum_coupling_loss                           = 70;
config.parallel_toolbox_installed                      =[];%Added

% Network geometry
config.nr_eNodeB_rings                                 = 1;
config.nr_sectors                                      = 1;
config.min_UE_eNodeB_distance                          = 35;     % in [m];   3D_UMa_fading: 35[m]; 3D_UMi_fading: 10[m]; 
config.inter_eNodeB_distance                           = 500;    % in [m];   3D_UMa_fading: 500[m]; 3D_UMi_fading: 200[m];
config.antenna_azimuth_offsett                         = 60;     % % in [deg];  The boresight direction of the eNodeB; This is a desired parameter, three-sector eNodeBs assumed 
config.map_resolution                                  = 5;
config.tx_height                                       = 25;     % in [m];   3D_UMa_fading: 25[m]; 3D_UMi_fading: 10[m];
config.rx_height                                       = 1.5;

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
config.electrical_downtilt                           = 102;
config.mechanical_downtilt                           = 0;
config.mechanical_slant                              = 0;
config.antenna.max_antenna_gain                      = 8; % As defined in TR 36.873 (2000 MHz)

% Load specific channel model parameters from Table 7.3-6 and system
% parameters for the FFT of the channel impulse response
config = load_specific_params(config);

config.simulation_time_blocks                        = 50;           % Design parameter: the simulation length as multiple of Block Fading length (config.BF_length)
config.results_folder                                = './results';
config.results_file                                  = 'auto';       % NOTE: 'auto' assigns a filename automatically

%% Generate network layout
% Generate a hexagonal grid of eNodeB sites
fprintf('Generating eNodeB positions');
eNodeB_site_positions = network_geometry.hexagonal_eNodeB_grid(config);
if config.nr_eNodeB_rings == 0;
    roi_x = [-config.inter_eNodeB_distance,config.inter_eNodeB_distance];
    roi_y = [-config.inter_eNodeB_distance,config.inter_eNodeB_distance];
else
    ROI_increase_factor = 0.5;
    tx_pos = eNodeB_site_positions;
    roi_x = [min(tx_pos(:,1)),max(tx_pos(:,1))];
    roi_y = [min(tx_pos(:,2)),max(tx_pos(:,2))];
     roi_x = roi_x + ROI_increase_factor*abs(roi_x(2)-roi_x(1))*[-1,1];
    roi_y = roi_y + ROI_increase_factor*abs(roi_y(2)-roi_y(1))*[-1,1];
%     roi_x = roi_x + ROI_increase_factor*abs(roi_x(2)-roi_x(1))*[-1,1];
%     roi_y = roi_y + ROI_increase_factor*abs(roi_y(2)-roi_y(1))*[-1,1];
end

if config.show_network_plots
    figure(1)
    plot(eNodeB_site_positions(:,1),eNodeB_site_positions(:,2),'bo','linewidth',2,'Markersize',12)
    xlim(roi_x);
    ylim(roi_y);
    hold on
    grid on
end

config.sector_azimuths = 0:360/config.nr_sectors:359; 
s_idx = 1;

%% Generate eNodeBs object and its properties
% The eNodeB site properties are: 
% * eNodeB_sites.id
% * eNodeB_sites.pos
% * eNodeB_sites.site_type
% * eNodeB_sites.sectors
% 
% The eNodeBs properties are:
% * eNodeBs.id
% * eNodeBs.parent_eNodeB
% * eNodeBs.pos
% * eNodeBs.nTX
% * eNodeBs.tx_height
% * eNodeBs.boresight
% * eNodeBs.antenna_type
eNodeB_sites = network_elements.eNodeB_site;
eNodeBs = network_elements.eNodeB;
for b_ = 1:size(eNodeB_site_positions,1)
    eNodeB_sites(b_).id = b_;
    eNodeB_sites(b_).pos = eNodeB_site_positions(b_,:);
    eNodeB_sites(b_).site_type = 'macro';
    eNodeBs(s_idx) = network_elements.eNodeB;
    for s_ = 1:length(config.sector_azimuths)
        eNodeBs(s_idx).id = s_;
        eNodeBs(s_idx).eNodeB_id = s_idx;
        eNodeBs(s_idx).parent_eNodeB = eNodeB_sites(b_);
        eNodeBs(s_idx).pos = eNodeB_site_positions(b_,:);
        eNodeBs(s_idx).max_power = config.eNodeB_tx_power;
        eNodeBs(s_idx).nTX = config.eNodeB_nTX;      
        eNodeBs(s_idx).tx_height = config.tx_height;
        eNodeBs(s_idx).boresight = mod(config.antenna_azimuth_offsett + config.sector_azimuths(s_),360);%Changed 
%         eNodeBs(s_idx).boresight = mod(config.antenna_azimuth_offsett + config.sector_azimuths(s_),360);%Changed 
        eNodeBs(s_idx).antenna_type  = config.antenna.antenna_gain_pattern;        
        eNodeBs(s_idx) = antennas.antenna.attach_antenna_to_eNodeB(eNodeBs(s_idx),config); 

        eNodeBs(s_idx).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
            config,...
            config.macroscopic_pathloss_model,...
            config.frequency,...
            config.macroscopic_pathloss_model_settings);          
        s_idx = s_idx + 1;        
    end
    eNodeB_sites(b_).sectors = eNodeBs((b_-1)*length(config.sector_azimuths) + (1:length(config.sector_azimuths)));
end
N_BS = s_idx-1;

% initialize clock
networkClock = network_elements.clock(config.BF_length);

for bb = 1:N_BS 
    eNodeBs(bb).clock = networkClock;
    eNodeBs(bb).neighbors_eNodeB = eNodeBs([1:(bb-1) (bb+1):length(eNodeBs)]);
end


%% Generate the pathloss map for each eNodeB sector
% Generate new maps that determine:
% * UE_indoor_map   - Pregenerates a map of the same size as pathloss map that indicates the indoor/outdoor condition for each pixel in the map: Boolean operator: 1->indoors; 0->outdoors
% * UE_height_map   - Pregenerates a map of the same size as pathloss map that indicates the height in [m] for each pixel
% * UE_indoor_distance_map  - Pregenerates a map of the same size as pathloss map that indicates the indoor distance in [m] for the pixels that are determined to be indoors

networkMacroscopicPathlossMap = channel_gain_wrappers.macroscopicPathlossMap;
networkMacroscopicPathlossMap.data_res = config.map_resolution;
networkMacroscopicPathlossMap.roi_x = roi_x;
networkMacroscopicPathlossMap.roi_y = roi_y;

networkMacroscopicPathlossMap.UE_indoor_map = macroscopic_pathloss_models.generalPathlossModel.calculate_UE_indoor_map(config,networkMacroscopicPathlossMap);
networkMacroscopicPathlossMap.UE_height_map = macroscopic_pathloss_models.generalPathlossModel.calculate_UE_height_map(config,networkMacroscopicPathlossMap);
networkMacroscopicPathlossMap.UE_indoor_distance_map = macroscopic_pathloss_models.generalPathlossModel.calculate_UE_indoor_distance(networkMacroscopicPathlossMap);
macroscopic_pathloss_models.generalPathlossModel.calculate_pathloss_maps(config,eNodeB_sites,networkMacroscopicPathlossMap);
networkPathlossMap = networkMacroscopicPathlossMap;


%% Generate shadow fading map
switch config.shadow_fading_type
    case 'claussen_3D_UMa'
        [config.roi_x,config.roi_y] = networkPathlossMap.valid_range;
        if config.debug_level>=1
            fprintf('Generating Claussen space-correlated shadow fading maps for 3D UMa (combining O-to-I, LOS and NLOS maps), ');
        end
        fprintf('One map per site');
        nShadowFadingMaps = size(eNodeB_site_positions,1);
        
        networkShadowFadingMap = channel_gain_wrappers.shadowFadingMapClaussen(...
            config.shadow_fading_map_resolution,...
            config.roi_x,...
            config.roi_y,...
            config.shadow_fading_n_neighbors,...
            nShadowFadingMaps,...
            config.shadow_fading_mean,...
            config.shadow_fading_sd_LOS,...
            config.r_eNodeBs,...
            config.claussen_stream,...
            config.deactivate_claussen_spatial_correlation,...
            config.parallel_toolbox_installed);
        
        networkShadowFadingMap_NLOS =  channel_gain_wrappers.shadowFadingMapClaussen(...
            config.shadow_fading_map_resolution,...
            config.roi_x,...
            config.roi_y,...
            config.shadow_fading_n_neighbors,...
            nShadowFadingMaps,...
            config.shadow_fading_mean,...
            config.shadow_fading_sd_NLOS,...
            config.r_eNodeBs,...
            config.claussen_stream,...
            config.deactivate_claussen_spatial_correlation,...
            config.parallel_toolbox_installed);
        
        networkShadowFadingMap_OTOI =  channel_gain_wrappers.shadowFadingMapClaussen(...
            config.shadow_fading_map_resolution,...
            config.roi_x,...
            config.roi_y,...
            config.shadow_fading_n_neighbors,...
            nShadowFadingMaps,...
            config.shadow_fading_mean,...
            config.shadow_fading_sd_OTOI,...
            config.r_eNodeBs,...
            config.claussen_stream,...
            config.deactivate_claussen_spatial_correlation,...
            config.parallel_toolbox_installed);
        
        % Merge shadowing maps for OTOI/LOS / NLOS:
        for i = 1:size(eNodeB_site_positions,1)
            UE_indoor_map_(:,:,i)= networkPathlossMap.UE_indoor_map;
        end
        networkShadowFadingMap.pathloss = UE_indoor_map_.*networkShadowFadingMap_OTOI.pathloss + ~UE_indoor_map_.*(...
            networkPathlossMap.LOS_map.*networkShadowFadingMap.pathloss + ~networkPathlossMap.LOS_map.*networkShadowFadingMap_NLOS.pathloss);
        
        networkShadowFadingMap.oneMapPerSite = ~config.decouple_site_shadow_fading_maps;
        

    case 'claussen_3D_UMi'
        [config.roi_x,config.roi_y] = networkPathlossMap.valid_range; %Changed from 'LTE_init_network_generation' line 213
        if config.debug_level>=1
            fprintf('Generating Claussen space-correlated shadow fading map for 3D UMi (combining O-to-I, LOS and NLOS maps), '); %Added O-to-I
        end
        
        fprintf('One map per site');
        nShadowFadingMaps = size(eNodeB_site_positions,1)+nRRHs;
        
        networkShadowFadingMap =  channel_gain_wrappers.shadowFadingMapClaussen(...
            config.shadow_fading_map_resolution,...
            config.roi_x,...
            config.roi_y,...
            config.shadow_fading_n_neighbors,...
            nShadowFadingMaps,...
            config.shadow_fading_mean,...
            config.shadow_fading_sd_LOS,...
            config.r_eNodeBs,...
            config.claussen_stream,...
            config.deactivate_claussen_spatial_correlation,...
            config.parallel_toolbox_installed);
        
        networkShadowFadingMap_NLOS = channel_gain_wrappers.shadowFadingMapClaussen(...
            config.shadow_fading_map_resolution,...
            config.roi_x,...
            config.roi_y,...
            config.shadow_fading_n_neighbors,...
            nShadowFadingMaps,...
            config.shadow_fading_mean,...
            config.shadow_fading_sd_NLOS,...
            config.r_eNodeBs,...
            config.claussen_stream,...
            config.deactivate_claussen_spatial_correlation,...
            config.parallel_toolbox_installed);
        
        networkShadowFadingMap_OTOI =  channel_gain_wrappers.shadowFadingMapClaussen(...
            config.shadow_fading_map_resolution,...
            config.roi_x,...
            config.roi_y,...
            config.shadow_fading_n_neighbors,...
            nShadowFadingMaps,...
            config.shadow_fading_mean,...
            config.shadow_fading_sd_OTOI,...
            config.r_eNodeBs,...
            config.claussen_stream,...
            config.deactivate_claussen_spatial_correlation,...
            config.parallel_toolbox_installed);
        
        % Merge shadowing maps for OTOI / LOS / NLOS:
        for i = 1:size(eNodeB_site_positions,1)
            UE_indoor_map_(:,:,i)= networkPathlossMap.UE_indoor_map;
        end
        
        networkShadowFadingMap.pathloss = UE_indoor_map_.*networkShadowFadingMap_OTOI.pathloss + ~UE_indoor_map_.*(...
            networkPathlossMap.LOS_map.*networkShadowFadingMap.pathloss + ~networkPathlossMap.LOS_map.*networkShadowFadingMap_NLOS.pathloss);
        networkShadowFadingMap.oneMapPerSite = ~config.decouple_site_shadow_fading_maps;
        
        for rrh_=1:nRRHs
            RRHs(rrh_).site_id = rrh_+length(sites);
        end
        
    case 'none'
        [config.roi_x, config.roi_y] = networkPathlossMap.valid_range; %Changed from 'LTE_init_network_generation' line 278
        if config.debug_level>=1
            fprintf('Generating dummy shadow fading map (i.e. no shadow fading');
        end
        networkShadowFadingMap = channel_gain_wrappers.shadowFadingDummyMap(...
            config.roi_x,...
            config.roi_y,...
            length(eNodeBs));
    otherwise
        error('%s shadow fading type not supported. Only "claussen" and "none" supported',config.shadow_fading_type);
end

networkPathlossMap.apply_MCL(config);

if config.show_network_plots
    figure(2)
    temp_gain = networkPathlossMap.pathloss(:,:,1);
    imagesc(networkShadowFadingMap.roi_x,networkShadowFadingMap.roi_y,temp_gain/max(max(temp_gain))*255)
    set(gca,'YDir','normal');
    axis equal
end

% With shadow fading
[networkPathlossMap.capacity,...
    networkPathlossMap.SINR,...
    networkPathlossMap.sector_assignment,...
    networkPathlossMap.maxSINR_assignment,...
    networkPathlossMap.diff_SINR_dB,...
    networkPathlossMap.sector_sizes,...
    networkPathlossMap.sector_centers] = LTE_common_calculate_cell_capacity(config,networkPathlossMap,eNodeB_sites,eNodeBs,networkShadowFadingMap);
% Without shadow fading
[networkPathlossMap.capacity2,...
    networkPathlossMap.SINR2,...
    networkPathlossMap.sector_assignment2,...
    networkPathlossMap.maxSINR_assignment,...
    networkPathlossMap.diff_SINR_dB2,...
    networkPathlossMap.sector_sizes2,...
    networkPathlossMap.sector_centers2] = LTE_common_calculate_cell_capacity(config,networkPathlossMap,eNodeB_sites,eNodeBs);

if config.show_network_plots
    figure(3)
    temp_gain = networkShadowFadingMap.pathloss(:,:,1);
    image(networkShadowFadingMap.roi_x,networkShadowFadingMap.roi_y,temp_gain/max(max(temp_gain))*255)
    set(gca,'YDir','normal');
    axis equal
end

if N_BS>1
    if config.show_network_plots
        figure(4)
        imagesc(networkShadowFadingMap.roi_x,networkShadowFadingMap.roi_y,networkPathlossMap.sector_assignment);
        set(gca,'YDir','normal');
        axis equal
    end
end
%% Generate UE positions
% Generates the UE properties:
% * |UEs.pos|                  - The UE position in [x, y]
% * |UEs.rx_height|            - The UE height in [m]
% * |UEs.direction|            - The Ue direction in [deg]     
% * |UEs.nRX|                  - The number of antenna ports at UE
% * |UEs.is_LOS|               - Boolean operator that indicates the LOS/NLOS condition of the UE; 1->LOS; 0->NLOS
% * |UEs.is_indoor|            - Boolean operator that indicates the indoor/outdoor condition of the UE; 1->indoors; 0->outdoors
% * |UEs.dist_indoor|          - The indoor distance in [m] from the UE postion to the outer wall of the building where the UE is located
% * |UEs.H_0_channel_trace|    - The channel impulse response from the channel model
% * |UEs.sampled_channel_H_0|  - The sampled channel impulse response
% * |UEs.H_0_final|            - The channel transfer function 

fprintf('Generating UE track');
UE_spatial_distribution = spatial_distributions.constantElementsPerCellSpatialDistribution(networkPathlossMap, config.UE_per_eNodeB);
UEs = network_elements.UE;
UE_positions = UE_spatial_distribution.generate_positions;
N_UE = size(UE_positions,1);
data_res = config.map_resolution;
for u_ = 1:N_UE 
    UEs(u_)     = network_elements.UE;
    UEs(u_).id  = u_;
    UEs(u_).pos = UE_positions(u_,:);
    UEs(u_).pos_pixel = LTE_common_pos_to_pixel(UEs(u_).pos,[roi_x(1) roi_y(1)], data_res);
    UEs(u_).rx_height = config.rx_height;
    UEs(u_).direction = floor(random('unif',0,359)); 
    UEs(u_).nRX = config.UE_nRX;
    UEs(u_).is_indoor = false; 
    UEs(u_).dist_indoor =false; 
    UEs(u_).H_0_channel_trace = [];
    UEs(u_).sampled_channel_H_0 = []; 
    UEs(u_).H_0_final = []; 
    UEs(u_).H_i_channel_trace = [];
    UEs(u_).sampled_channel_H_i = [];
    UEs(u_).H_i_after_fft = [];
    UEs(u_).H_i_full_final =[];
    UEs(u_).TTI_of_smallscale_fading_recalculation =[];
    UEs(u_).recalculate_3D_smallscale_fading = true;
    
    [site_id, sector_num, eNodeB_id ] = networkPathlossMap.cell_assignment(UEs(u_).pos);
    eNodeBs(eNodeB_id).attachUser(UEs(u_));  
   
    if config.show_network_plots
        figure(1)
        plot(UEs(u_).pos(1),UEs(u_).pos(2),'rx','linewidth',1,'Markersize',7);
        xlim(roi_x);
        axis equal;
    end
    UEs(u_).channels = channel_gain_wrappers.TR36873_Fading_3D_Channel(config);
end

      


%% Start main simulation loop that outputs the channel impulse response
% Calculate correlated Large Scale Parameters
% Generate the channel impulse response and the channel transfer function

while networkClock.current_block < config.simulation_time_blocks 
    % Advance the network clock
    networkClock.advance_1_TTI;
    if config.debug_level>=1
        fprintf(' TTI %5.0f/%d: ',networkClock.current_block,config.simulation_time_blocks);
    end
    
    UE_positions = [UEs.pos];
    for bb = 1:N_BS
        %% Generate correlated Large Scale Parameters
        if networkClock.current_block==1
            fprintf('Calculating large scale parameters for eNodeb %d \n', eNodeBs(bb).eNodeB_id);
            UE_positions = zeros(eNodeBs(bb).attached_UEs, 2);
            data_res = config.map_resolution;
            % Get UEs that are attached to this eNodeB and check whether they are in LOS or NLOS
            % Generate LSPs for interfering links
            if eNodeBs(bb).attached_UEs > 0
                for jj=1:eNodeBs(bb).attached_UEs
                    UE_distance_to_site = zeros(length(eNodeBs(bb).neighbors_eNodeB)+1,1);
                    eNodeBs(bb).attached_UEs_vector(jj).is_LOS = zeros(length(eNodeBs(bb).neighbors_eNodeB)+1,1);
                    eNodeBs(bb).attached_UEs_vector(jj).all_large_scale_params = zeros(length(eNodeBs(bb).neighbors_eNodeB)+1,7);
                    UE_positions(jj,:)                    = eNodeBs(bb).attached_UEs_vector(jj).pos;
                    UE_pos_pixel                          = LTE_common_pos_to_pixel(UE_positions(jj,:), [networkPathlossMap.roi_x(1), networkPathlossMap.roi_y(1)], data_res);
                    eNodeBs(bb).attached_UEs_vector(jj).is_LOS(1)    = networkPathlossMap.LOS_map(UE_pos_pixel(2),UE_pos_pixel(1), eNodeBs(bb).parent_eNodeB.id);
                    eNodeBs(bb).attached_UEs_vector(jj).rx_height = networkPathlossMap.UE_height_map(UE_pos_pixel(2),UE_pos_pixel(1));
                    eNodeBs(bb).attached_UEs_vector(jj).is_indoor = networkPathlossMap.UE_indoor_map(UE_pos_pixel(2),UE_pos_pixel(1));
                    eNodeBs(bb).attached_UEs_vector(jj).dist_indoor = networkPathlossMap.UE_indoor_distance_map(UE_pos_pixel(2),UE_pos_pixel(1));
                    eNodeBs(bb).attached_UEs_vector(jj).recalculate_3D_smallscale_fading = true; % In case the large scale parameters have changed, the fast scale fading is recalculated.
                    if ~isempty(eNodeBs(bb).neighbors_eNodeB)
                        for ii= 1:length(eNodeBs(bb).neighbors_eNodeB)
                            if config.nr_eNodeB_rings == 0
                                eNodeBs(bb).attached_UEs_vector(jj).is_LOS(ii+1)    = networkPathlossMap.LOS_map(UE_pos_pixel(2),UE_pos_pixel(1));
                                eNodeBs(bb).attached_UEs_vector(jj).recalculate_3D_smallscale_fading_i(ii) = true;
                            else
                                eNodeBs(bb).attached_UEs_vector(jj).is_LOS(ii+1)    = networkPathlossMap.LOS_map(UE_pos_pixel(2),UE_pos_pixel(1), eNodeBs(bb).neighbors_eNodeB(ii).id);
                                eNodeBs(bb).attached_UEs_vector(jj).recalculate_3D_smallscale_fading_i(ii) = true;
                            end
                            % Calculate large scale parameters for
                            % interfering links
                            % O-to-I
                            if eNodeBs(bb).attached_UEs_vector(jj).is_indoor
                                UE_distance_to_site(ii+1) = abs(sqrt(sum((eNodeBs(bb).neighbors_eNodeB(ii).parent_eNodeB.pos - eNodeBs(bb).attached_UEs_vector(jj).pos).^2)) - eNodeBs(bb).attached_UEs_vector(jj).dist_indoor); % Consider only outdoor distance for indoor UEs
                                if eNodeBs(bb).attached_UEs_vector(jj).is_LOS(1)  % indoor UE case doesn't change (interf. vs. desired channel)
                                    [ZS_D_mu_OTOI_LOS,ZS_D_sigma_OTOI_LOS,ZOD_mu_offset_OTOI_LOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_OTOI_LOS(config, eNodeBs(bb).attached_UEs_vector(jj).rx_height, UE_distance_to_site(ii+1));
                                    eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(ii+1,:) =  [ZS_D_mu_OTOI_LOS,ZS_D_sigma_OTOI_LOS,ZOD_mu_offset_OTOI_LOS];
                                else
                                    [ZS_D_mu_OTOI_NLOS,ZS_D_sigma_OTOI_NLOS,ZOD_mu_offset_OTOI_NLOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_OTOI_NLOS(config, eNodeBs(bb).attached_UEs_vector(jj).rx_height, UE_distance_to_site(ii+1));
                                    eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(ii+1,:) =  [ZS_D_mu_OTOI_NLOS,ZS_D_sigma_OTOI_NLOS,ZOD_mu_offset_OTOI_NLOS];
                                end
                                [sigmas_SF_OTOI, non_value_param, sigmas_DS_OTOI, sigmas_ASD_OTOI, sigmas_ASA_OTOI, sigmas_ZSD_OTOI, sigmas_ZSA_OTOI] = eNodeBs(bb).sigmas_OTOI(config, randn(7,1), eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(ii+1,:));
                                eNodeBs(bb).attached_UEs_vector(jj).all_large_scale_params(ii+1,:) = [sigmas_SF_OTOI, non_value_param, sigmas_DS_OTOI, sigmas_ASD_OTOI, sigmas_ASA_OTOI, sigmas_ZSD_OTOI, sigmas_ZSA_OTOI];
                                
                            else
                                % LOS_i
                                if eNodeBs(bb).attached_UEs_vector(jj).is_LOS(ii+1)
                                    UE_distance_to_site(ii+1) = sqrt(sum((eNodeBs(bb).neighbors_eNodeB(ii).parent_eNodeB.pos - eNodeBs(bb).attached_UEs_vector(jj).pos).^2));
                                    [ZS_D_mu_LOS,ZS_D_sigma_LOS,ZOD_mu_offset_LOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_LOS(config, eNodeBs(bb).attached_UEs_vector(jj).rx_height, UE_distance_to_site(ii+1));
                                    eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(ii+1,:) = [ZS_D_mu_LOS,ZS_D_sigma_LOS,ZOD_mu_offset_LOS];
                                    [sigmas_SF_LOS, sigmas_KF_LOS, sigmas_DS_LOS, sigmas_ASD_LOS, sigmas_ASA_LOS, sigmas_ZSD_LOS, sigmas_ZSA_LOS] = eNodeBs(bb).sigmas_LOS(config, randn(7,1), eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(ii+1,:));
                                    eNodeBs(bb).attached_UEs_vector(jj).all_large_scale_params(ii+1,:) = [sigmas_SF_LOS, sigmas_KF_LOS, sigmas_DS_LOS, sigmas_ASD_LOS, sigmas_ASA_LOS, sigmas_ZSD_LOS, sigmas_ZSA_LOS];
                                    
                                    % NLOS_i
                                else
                                    UE_distance_to_site(ii+1) = sqrt(sum((eNodeBs(bb).neighbors_eNodeB(ii).parent_eNodeB.pos - eNodeBs(bb).attached_UEs_vector(jj).pos).^2));
                                    [ZS_D_mu_NLOS,ZS_D_sigma_NLOS,ZOD_mu_offset_NLOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_NLOS(config, eNodeBs(bb).attached_UEs_vector(jj).rx_height, UE_distance_to_site(ii+1));
                                    eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(ii+1,:) =  [ZS_D_mu_NLOS,ZS_D_sigma_NLOS,ZOD_mu_offset_NLOS];
                                    [sigmas_SF_NLOS, non_value_param, sigmas_DS_NLOS, sigmas_ASD_NLOS, sigmas_ASA_NLOS, sigmas_ZSD_NLOS, sigmas_ZSA_NLOS] = eNodeBs(bb).sigmas_NLOS(config, randn(7,1), eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(ii+1,:));
                                    eNodeBs(bb).attached_UEs_vector(jj).all_large_scale_params(ii+1,:) = [sigmas_SF_NLOS, non_value_param, sigmas_DS_NLOS, sigmas_ASD_NLOS, sigmas_ASA_NLOS, sigmas_ZSD_NLOS, sigmas_ZSA_NLOS];
                                    
                                end
                            end
                        end
                    end
                end
            end
            
            % Generate LSP for desired channel: if more than one UE is attached to the BS sector perform cross-correlation
            if eNodeBs(bb).attached_UEs > 1
                [filtered_rand_grid_OTOI, filtered_rand_grid_LOS, filtered_rand_grid_NLOS, min_roi_position, ~] = ...
                    calculate_correlations(UE_positions, data_res, config);
                
                % Apply spatial filter at each UE.
                for jj=1:eNodeBs(bb).attached_UEs
                    % Generate UE indices from grid
                    [UE_index] = LTE_common_pos_to_pixel(eNodeBs(bb).attached_UEs_vector(jj).pos, min_roi_position, data_res);
                    UE_index_column = UE_index(2);
                    UE_index_row    = UE_index(1);
                    % Calculate large scale parameters for signal link
                    % OTOI
                    if eNodeBs(bb).attached_UEs_vector(jj).is_indoor
                        filtered_matrix                 = squeeze(filtered_rand_grid_OTOI(UE_index_row, UE_index_column,:));
                        correlated_large_scale_params   = eNodeBs(bb).calculate_correlated_large_scale_param_OTOI(config, filtered_matrix);
                        UE_distance_to_site(1) = abs(sqrt(sum((eNodeBs(bb).parent_eNodeB.pos - eNodeBs(bb).attached_UEs_vector(jj).pos).^2)) - eNodeBs(bb).attached_UEs_vector(jj).dist_indoor); % consider outdoor distance
                        if eNodeBs(bb).attached_UEs_vector(jj).is_LOS(1)
                            [ZS_D_mu_OTOI_LOS,ZS_D_sigma_OTOI_LOS,ZOD_mu_offset_OTOI_LOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_OTOI_LOS(config, eNodeBs(bb).attached_UEs_vector(jj).rx_height, UE_distance_to_site(1));
                            eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(1,:) =  [ZS_D_mu_OTOI_LOS, ZS_D_sigma_OTOI_LOS, ZOD_mu_offset_OTOI_LOS];
                        else
                            [ZS_D_mu_OTOI_NLOS,ZS_D_sigma_OTOI_NLOS,ZOD_mu_offset_OTOI_NLOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_OTOI_NLOS(config, eNodeBs(bb).attached_UEs_vector(jj).rx_height, UE_distance_to_site(1));
                            eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(1,:) =  [ZS_D_mu_OTOI_NLOS, ZS_D_sigma_OTOI_NLOS, ZOD_mu_offset_OTOI_NLOS];
                        end
                        [sigmas_SF_OTOI, non_value_param, sigmas_DS_OTOI, sigmas_ASD_OTOI, sigmas_ASA_OTOI, sigmas_ZSD_OTOI, sigmas_ZSA_OTOI] = eNodeBs(bb).sigmas_OTOI(config, correlated_large_scale_params, eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(1,:));
                        eNodeBs(bb).attached_UEs_vector(jj).all_large_scale_params(1,:) = [sigmas_SF_OTOI, non_value_param, sigmas_DS_OTOI, sigmas_ASD_OTOI, sigmas_ASA_OTOI, sigmas_ZSD_OTOI, sigmas_ZSA_OTOI];
                        
                    else
                        % LOS
                        if eNodeBs(bb).attached_UEs_vector(jj).is_LOS(1)
                            filtered_matrix                 = squeeze(filtered_rand_grid_LOS(UE_index_row, UE_index_column,:));
                            correlated_large_scale_params   = eNodeBs(bb).calculate_correlated_large_scale_param_LOS(config, filtered_matrix);
                            UE_distance_to_site = sqrt(sum((eNodeBs(bb).parent_eNodeB.pos - eNodeBs(bb).attached_UEs_vector(jj).pos).^2));
                            [ZS_D_mu_LOS,ZS_D_sigma_LOS,ZOD_mu_offset_LOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_LOS(config, eNodeBs(bb).attached_UEs_vector(jj).rx_height, UE_distance_to_site(1));
                            eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(1,:) = [ZS_D_mu_LOS,ZS_D_sigma_LOS,ZOD_mu_offset_LOS];
                            [sigmas_SF_LOS, sigmas_KF_LOS, sigmas_DS_LOS, sigmas_ASD_LOS, sigmas_ASA_LOS, sigmas_ZSD_LOS, sigmas_ZSA_LOS] = eNodeBs(bb).sigmas_LOS(config, correlated_large_scale_params, eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(1,:));
                            eNodeBs(bb).attached_UEs_vector(jj).all_large_scale_params(1,:) = [sigmas_SF_LOS, sigmas_KF_LOS, sigmas_DS_LOS, sigmas_ASD_LOS, sigmas_ASA_LOS, sigmas_ZSD_LOS, sigmas_ZSA_LOS];
                            % NLOS
                        else
                            filtered_matrix                 = squeeze(filtered_rand_grid_NLOS(UE_index_row, UE_index_column,:));
                            correlated_large_scale_params   = eNodeBs(bb).calculate_correlated_large_scale_param_NLOS(config, filtered_matrix);
                            UE_distance_to_site = sqrt(sum((eNodeBs(bb).parent_eNodeB.pos - eNodeBs(bb).attached_UEs_vector(jj).pos).^2));
                            [ZS_D_mu_NLOS,ZS_D_sigma_NLOS,ZOD_mu_offset_NLOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_NLOS(config, eNodeBs(bb).attached_UEs_vector(jj).rx_height, UE_distance_to_site(1));
                            eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(1,:) =  [ZS_D_mu_NLOS,ZS_D_sigma_NLOS,ZOD_mu_offset_NLOS];
                            [sigmas_SF_NLOS, non_value_param, sigmas_DS_NLOS, sigmas_ASD_NLOS, sigmas_ASA_NLOS, sigmas_ZSD_NLOS, sigmas_ZSA_NLOS] = eNodeBs(bb).sigmas_NLOS(config, correlated_large_scale_params, eNodeBs(bb).attached_UEs_vector(jj).all_ZOD_params(1,:));
                            eNodeBs(bb).attached_UEs_vector(jj).all_large_scale_params(1,:) = [sigmas_SF_NLOS, non_value_param, sigmas_DS_NLOS, sigmas_ASD_NLOS, sigmas_ASA_NLOS, sigmas_ZSD_NLOS, sigmas_ZSA_NLOS];
                        end
                    end
                end
                
            % If only one UE is attached to the eNodeB than no cross-correlation is assumed (based on TR 36.873)
            else
                if ~isempty(eNodeBs(bb).attached_UEs_vector)
                    if eNodeBs(bb).attached_UEs_vector(1).is_indoor
                        UE_distance_to_site(1) = abs(sqrt(sum((eNodeBs(bb).parent_eNodeB.pos - eNodeBs(bb).attached_UEs_vector(1).pos).^2)) - eNodeBs(bb).attached_UEs_vector(1).dist_indoor);
                        if eNodeBs(bb).attached_UEs_vector(1).is_LOS(1)
                            [ZS_D_mu_OTOI_LOS,ZS_D_sigma_OTOI_LOS,ZOD_mu_offset_OTOI_LOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_OTOI_LOS(config, eNodeBs(bb).attached_UEs_vector(1).rx_height, UE_distance_to_site(1));
                            eNodeBs(bb).attached_UEs_vector(1).all_ZOD_params(1,:) =  [ZS_D_mu_OTOI_LOS,ZS_D_sigma_OTOI_LOS,ZOD_mu_offset_OTOI_LOS];
                        else
                            [ZS_D_mu_OTOI_NLOS,ZS_D_sigma_OTOI_NLOS,ZOD_mu_offset_OTOI_NLOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_OTOI_NLOS(config, eNodeBs(bb).attached_UEs_vector(1).rx_height, UE_distance_to_site(1));
                            eNodeBs(bb).attached_UEs_vector(1).all_ZOD_params(1,:) =  [ZS_D_mu_OTOI_NLOS,ZS_D_sigma_OTOI_NLOS,ZOD_mu_offset_OTOI_NLOS];
                        end
                        [sigmas_SF_OTOI, non_value_param, sigmas_DS_OTOI, sigmas_ASD_OTOI, sigmas_ASA_OTOI, sigmas_ZSD_OTOI, sigmas_ZSA_OTOI] = eNodeBs(bb).sigmas_OTOI(config, randn(7,1),eNodeBs(bb).attached_UEs_vector(1).all_ZOD_params(1,:) );
                        eNodeBs(bb).attached_UEs_vector(1).all_large_scale_params(1,:) = [sigmas_SF_OTOI, non_value_param, sigmas_DS_OTOI, sigmas_ASD_OTOI, sigmas_ASA_OTOI, sigmas_ZSD_OTOI, sigmas_ZSA_OTOI];
                    else
                        UE_distance_to_site = sqrt(sum((eNodeBs(bb).parent_eNodeB.pos - eNodeBs(bb).attached_UEs_vector(jj).pos).^2));
                        if eNodeBs(bb).attached_UEs_vector(1).is_LOS(1)
                            [ZS_D_mu_LOS,ZS_D_sigma_LOS,ZOD_mu_offset_LOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_LOS(config, eNodeBs(bb).attached_UEs_vector(1).rx_height, UE_distance_to_site);
                            eNodeBs(bb).attached_UEs_vector(1).all_ZOD_params(1,:) = [ZS_D_mu_LOS,ZS_D_sigma_LOS,ZOD_mu_offset_LOS];
                            [sigmas_SF_LOS, sigmas_KF_LOS, sigmas_DS_LOS, sigmas_ASD_LOS, sigmas_ASA_LOS, sigmas_ZSD_LOS, sigmas_ZSA_LOS] = eNodeBs(bb).sigmas_LOS(config, randn(7,1), eNodeBs(bb).attached_UEs_vector(1).all_ZOD_params(1,:));
                            eNodeBs(bb).attached_UEs_vector(1).all_large_scale_params(1,:) = [sigmas_SF_LOS, sigmas_KF_LOS, sigmas_DS_LOS, sigmas_ASD_LOS, sigmas_ASA_LOS, sigmas_ZSD_LOS, sigmas_ZSA_LOS];
                            
                        else
                            [ZS_D_mu_NLOS,ZS_D_sigma_NLOS,ZOD_mu_offset_NLOS] = eNodeBs(bb).generate_ZSD_ZoD_offset_parameters_NLOS(config, eNodeBs(bb).attached_UEs_vector(1).rx_height, UE_distance_to_site);
                            eNodeBs(bb).attached_UEs_vector(1).all_ZOD_params(1,:) =  [ZS_D_mu_NLOS,ZS_D_sigma_NLOS,ZOD_mu_offset_NLOS];
                            [sigmas_SF_NLOS, non_value_param,  sigmas_DS_NLOS, sigmas_ASD_NLOS, sigmas_ASA_NLOS, sigmas_ZSD_NLOS, sigmas_ZSA_NLOS] = eNodeBs(bb).sigmas_NLOS(config, randn(7,1),eNodeBs(bb).attached_UEs_vector(1).all_ZOD_params(1,:) );
                            eNodeBs(bb).attached_UEs_vector(1).all_large_scale_params(1,:) = [sigmas_SF_NLOS, non_value_param, sigmas_DS_NLOS, sigmas_ASD_NLOS, sigmas_ASA_NLOS, sigmas_ZSD_NLOS, sigmas_ZSA_NLOS];
                        end
                    end
                end
                
            end
        end
        
        %% Generate channel coefficients for desired- and interfering links
        % DESIRED CHANNEL
        
     
        
        
        for uu = 1:eNodeBs(bb).attached_UEs
        fprintf('Generating desired and interfering channels for UE %d \n', eNodeBs(bb).attached_UEs_vector(uu).id);
            if eNodeBs(bb).attached_UEs_vector(uu).recalculate_3D_smallscale_fading
                % Generate the channel matrix for each UE
                if eNodeBs(bb).attached_UEs_vector(uu).is_indoor
                    eNodeBs(bb).attached_UEs_vector(uu).H_0_channel_trace = zeros (eNodeBs(bb).attached_UEs_vector(uu).nRX,  eNodeBs(bb).nTX, 1, config.simulation_time_blocks, config.NumClusters_OTOI);
                    eNodeBs(bb).attached_UEs_vector(uu).channels.calculate_channel_coefficient_OTOI(config, eNodeBs(bb).attached_UEs_vector(uu).rx_height, eNodeBs(bb).pos, eNodeBs(bb).attached_UEs_vector(uu).pos, eNodeBs(bb).attached_UEs_vector(uu).all_large_scale_params(1,:), eNodeBs(bb).attached_UEs_vector(uu).all_ZOD_params(1,:), eNodeBs(bb), eNodeBs(bb).nTX, eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).attached_UEs_vector(uu).is_indoor, 1, networkClock.current_block, config.simulation_time_blocks, networkClock.BF_time, eNodeBs(bb).attached_UEs_vector(uu).direction, eNodeBs(bb).attached_UEs_vector(uu).dist_indoor);
                    eNodeBs(bb).attached_UEs_vector(uu).H_0_channel_trace = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_channel_matrix_OTOI;
                    eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_0 = eNodeBs(bb).attached_UEs_vector(uu).channels.sampled_channel_OTOI(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).nTX, config.simulation_time_blocks, networkClock.current_block, eNodeBs(bb).attached_UEs_vector(uu).H_0_channel_trace);
                    
                else
                    if eNodeBs(bb).attached_UEs_vector(uu).is_LOS(1)
                        eNodeBs(bb).attached_UEs_vector(uu).H_0_channel_trace = zeros(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).nTX, 1, config.simulation_time_blocks, config.NumClusters_LOS);
                        eNodeBs(bb).attached_UEs_vector(uu).channels.calculate_channel_coefficient_LOS(config, eNodeBs(bb).attached_UEs_vector(uu).rx_height, eNodeBs(bb).pos, eNodeBs(bb).attached_UEs_vector(uu).pos, eNodeBs(bb).attached_UEs_vector(uu).all_large_scale_params(1,:), eNodeBs(bb).attached_UEs_vector(uu).all_ZOD_params(1,:), eNodeBs(bb), eNodeBs(bb).nTX, eNodeBs(bb).attached_UEs_vector(uu).nRX, 1, networkClock.current_block, config.simulation_time_blocks, networkClock.BF_time, eNodeBs(bb).attached_UEs_vector(uu).direction);
                        eNodeBs(bb).attached_UEs_vector(uu).H_0_channel_trace = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_channel_matrix_LOS;
                        eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_0 = eNodeBs(bb).attached_UEs_vector(uu).channels.sampled_channel_LOS(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).nTX, config.simulation_time_blocks, networkClock.current_block, eNodeBs(bb).attached_UEs_vector(uu).H_0_channel_trace);
                    else
                        eNodeBs(bb).attached_UEs_vector(uu).H_0_channel_trace = zeros(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).nTX, 1, config.simulation_time_blocks, config.NumClusters_NLOS);
                        eNodeBs(bb).attached_UEs_vector(uu).channels.calculate_channel_coefficient_NLOS(config, eNodeBs(bb).attached_UEs_vector(uu).rx_height, eNodeBs(bb).pos, eNodeBs(bb).attached_UEs_vector(uu).pos, eNodeBs(bb).attached_UEs_vector(uu).all_large_scale_params(1,:), eNodeBs(bb).attached_UEs_vector(uu).all_ZOD_params(1,:), eNodeBs(bb), eNodeBs(bb).nTX, eNodeBs(bb).attached_UEs_vector(uu).nRX, 1, networkClock.current_block, config.simulation_time_blocks, networkClock.BF_time, eNodeBs(bb).attached_UEs_vector(uu).direction);
                        eNodeBs(bb).attached_UEs_vector(uu).H_0_channel_trace = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_channel_matrix_NLOS;
                        eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_0 = eNodeBs(bb).attached_UEs_vector(uu).channels.sampled_channel_NLOS(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).nTX, config.simulation_time_blocks, networkClock.current_block, eNodeBs(bb).attached_UEs_vector(uu).H_0_channel_trace);
                    end
                end
                
                eNodeBs(bb).attached_UEs_vector(uu).recalculate_3D_smallscale_fading = false;
                eNodeBs(bb).attached_UEs_vector(uu).TTI_of_smallscale_fading_recalculation = networkClock.current_block;
                
           
                % get the proper value from pathloss map
                 pathloss_linear_0 = networkMacroscopicPathlossMap.get_pathloss_eNodeB(eNodeBs(bb).attached_UEs_vector(uu).pos_pixel, bb); 

                % Here goes the channel matrix after FFT with size [nRx nTx 1 nTTI nRB]
                H_0_after_fft = eNodeBs(bb).attached_UEs_vector(uu).channels.get_RB_trace(eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_0);
                eNodeBs(bb).attached_UEs_vector(uu).H_0_final = H_0_after_fft./ sqrt(pathloss_linear_0.*eNodeBs(bb).attached_UEs_vector(uu).all_large_scale_params(1,1));
            end
            H_0 = permute(eNodeBs(bb).attached_UEs_vector(uu).H_0_final(:,:,:,networkClock.current_block-eNodeBs(bb).attached_UEs_vector(uu).TTI_of_smallscale_fading_recalculation+1,:),[1,2,5,4,3]); % final H_0 with size [nRx nTx nRB nTTI]
            
            % INTERFERING CHANNELS
            for kk = 1:length(eNodeBs(bb).neighbors_eNodeB)     % indicates the interfering eNodeB index
                if eNodeBs(bb).attached_UEs_vector(uu).recalculate_3D_smallscale_fading_i(kk)
                    eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_i = zeros(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).neighbors_eNodeB(kk).nTX, 1, config.simulation_time_blocks-networkClock.current_block+1, 200);
                    eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft = zeros(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).neighbors_eNodeB(kk).nTX, 1, config.simulation_time_blocks-networkClock.current_block+1, 2*config.N_RB , length(eNodeBs(bb).neighbors_eNodeB));
                    % preallocate memory
                    if networkClock.current_block == 1
                        if kk ~= 1
                            eNodeBs(bb).attached_UEs_vector(uu).H_i_full_final(:,:,:,:,:,kk)=0;
                        end
                    else
                        if kk == 1
                            eNodeBs(bb).attached_UEs_vector(uu).H_i_full_final = zeros(size(eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft,1), size(eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft,2), size(eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft,3), size(eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft,4), size(eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft,5), size(eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft,6));
                        end
                    end
                    if eNodeBs(bb).attached_UEs_vector(uu).is_indoor
                        eNodeBs(bb).attached_UEs_vector(uu).H_i_channel_trace = zeros(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).neighbors_eNodeB(kk).nTX, 1, config.simulation_time_blocks-networkClock.current_block+1, config.NumClusters_OTOI);
                        eNodeBs(bb).attached_UEs_vector(uu).channels.calculate_channel_coefficient_OTOI(config, eNodeBs(bb).attached_UEs_vector(uu).rx_height, eNodeBs(bb).neighbors_eNodeB(kk).parent_eNodeB.pos, eNodeBs(bb).attached_UEs_vector(uu).pos, eNodeBs(bb).attached_UEs_vector(uu).all_large_scale_params(kk+1,:), eNodeBs(bb).attached_UEs_vector(uu).all_ZOD_params(kk+1,:), eNodeBs(bb).neighbors_eNodeB(kk), eNodeBs(bb).neighbors_eNodeB(kk).nTX, eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).attached_UEs_vector(uu).is_indoor, kk+1, networkClock.current_block, config.simulation_time_blocks, networkClock.BF_time, eNodeBs(bb).attached_UEs_vector(uu).direction, eNodeBs(bb).attached_UEs_vector(uu).dist_indoor);
                        eNodeBs(bb).attached_UEs_vector(uu).H_i_channel_trace(:,:,:,:,:,kk) = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_channel_matrix_OTOI;
                        eNodeBs(bb).attached_UEs_vector(uu).channels.sampled_channel_OTOI(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).neighbors_eNodeB(kk).nTX, config.simulation_time_blocks, networkClock.current_block, eNodeBs(bb).attached_UEs_vector(uu).H_i_channel_trace(:,:,:,:,:,kk));
                        delay_tap_dimension = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_delay_tap_dimension_OTOI;
                        eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_i(:,:,:,:,1:delay_tap_dimension,kk) = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_sampled_channel_OTOI;
                        eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft(:,:,:,:,:,kk) = eNodeBs(bb).attached_UEs_vector(uu).channels.get_RB_trace(eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_i(:,:,:,:,:,kk));
                    else
                        if eNodeBs(bb).attached_UEs_vector(uu).is_LOS(kk+1)
                            eNodeBs(bb).attached_UEs_vector(uu).H_i_channel_trace = zeros(eNodeBs(bb).attached_UEs_vector(uu).nRX,  eNodeBs(bb).neighbors_eNodeB(kk).nTX, 1, config.simulation_time_blocks-networkClock.current_block+1, config.NumClusters_LOS);
                            eNodeBs(bb).attached_UEs_vector(uu).channels.calculate_channel_coefficient_LOS(config, eNodeBs(bb).attached_UEs_vector(uu).rx_height,  eNodeBs(bb).neighbors_eNodeB(kk).parent_eNodeB.pos, eNodeBs(bb).attached_UEs_vector(uu).pos, eNodeBs(bb).attached_UEs_vector(uu).all_large_scale_params(kk+1,:), eNodeBs(bb).attached_UEs_vector(uu).all_ZOD_params(kk+1,:),  eNodeBs(bb).neighbors_eNodeB(kk),  eNodeBs(bb).neighbors_eNodeB(kk).nTX, eNodeBs(bb).attached_UEs_vector(uu).nRX, kk+1, networkClock.current_block, config.simulation_time_blocks, networkClock.BF_time, eNodeBs(bb).attached_UEs_vector(uu).direction);
                            eNodeBs(bb).attached_UEs_vector(uu).H_i_channel_trace(:,:,:,:,:,kk) = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_channel_matrix_LOS;
                            eNodeBs(bb).attached_UEs_vector(uu).channels.sampled_channel_LOS(eNodeBs(bb).attached_UEs_vector(uu).nRX,  eNodeBs(bb).neighbors_eNodeB(kk).nTX, config.simulation_time_blocks, networkClock.current_block, eNodeBs(bb).attached_UEs_vector(uu).H_i_channel_trace(:,:,:,:,:,kk));
                            delay_tap_dimension = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_delay_tap_dimension_LOS;
                            eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_i(:,:,:,:,1:delay_tap_dimension,kk) = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_sampled_channel_LOS;
                            eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft(:,:,:,:,:,kk) = eNodeBs(bb).attached_UEs_vector(uu).channels.get_RB_trace(eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_i(:,:,:,:,:,kk));
                        else
                            eNodeBs(bb).attached_UEs_vector(uu).H_i_channel_trace = zeros(eNodeBs(bb).attached_UEs_vector(uu).nRX,  eNodeBs(bb).neighbors_eNodeB(kk).nTX, 1, config.simulation_time_blocks-networkClock.current_block+1, config.NumClusters_NLOS);
                            eNodeBs(bb).attached_UEs_vector(uu).channels.calculate_channel_coefficient_NLOS(config, eNodeBs(bb).attached_UEs_vector(uu).rx_height,  eNodeBs(bb).neighbors_eNodeB(kk).parent_eNodeB.pos, eNodeBs(bb).attached_UEs_vector(uu).pos, eNodeBs(bb).attached_UEs_vector(uu).all_large_scale_params(kk+1,:), eNodeBs(bb).attached_UEs_vector(uu).all_ZOD_params(kk+1,:),  eNodeBs(bb).neighbors_eNodeB(kk),  eNodeBs(bb).neighbors_eNodeB(kk).nTX, eNodeBs(bb).attached_UEs_vector(uu).nRX, kk+1, networkClock.current_block, config.simulation_time_blocks, networkClock.BF_time, eNodeBs(bb).attached_UEs_vector(uu).direction);
                            eNodeBs(bb).attached_UEs_vector(uu).H_i_channel_trace(:,:,:,:,:,kk) = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_channel_matrix_NLOS;
                            eNodeBs(bb).attached_UEs_vector(uu).channels.sampled_channel_NLOS(eNodeBs(bb).attached_UEs_vector(uu).nRX,  eNodeBs(bb).neighbors_eNodeB(kk).nTX, config.simulation_time_blocks, networkClock.current_block, eNodeBs(bb).attached_UEs_vector(uu).H_i_channel_trace(:,:,:,:,:,kk));
                            delay_tap_dimension = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_delay_tap_dimension_NLOS;
                            eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_i(:,:,:,:,1:delay_tap_dimension,kk) = eNodeBs(bb).attached_UEs_vector(uu).channels.saved_sampled_channel_NLOS;
                            eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft(:,:,:,:,:,kk) = eNodeBs(bb).attached_UEs_vector(uu).channels.get_RB_trace(eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_i(:,:,:,:,:,kk));
                        end
                    end
                    eNodeBs(bb).attached_UEs_vector(uu).recalculate_3D_smallscale_fading_i(kk) = false;
                    eNodeBs(bb).attached_UEs_vector(uu).TTI_of_smallscale_fading_recalculation_i(kk) = networkClock.current_block;
                    
                     pathloss_linear_i = networkMacroscopicPathlossMap.get_pathloss_eNodeB(eNodeBs(bb).attached_UEs_vector(uu).pos_pixel, eNodeBs(bb).neighbors_eNodeB(kk).eNodeB_id); 

%                     pathloss_linear_i = networkPathlossMap.pathloss(eNodeBs(bb).attached_UEs_vector(uu).pos_pixel(1,1),eNodeBs(bb).attached_UEs_vector(uu).pos_pixel(1,2),eNodeBs(bb).neighbors_eNodeB(kk).eNodeB_id); 
                    % apply path-loss and shadow fading
                    eNodeBs(bb).attached_UEs_vector(uu).H_i_full_final(:,:,:,:,:,kk) = eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft(:,:,:,:,:,kk)./ sqrt(pathloss_linear_i*eNodeBs(bb).attached_UEs_vector(uu).all_large_scale_params(kk+1,1));
                end
                fprintf('%d ', kk+1);
                H_i_current_block = eNodeBs(bb).attached_UEs_vector(uu).H_i_full_final(:,:,:,networkClock.current_block-eNodeBs(bb).attached_UEs_vector(uu).TTI_of_smallscale_fading_recalculation_i(kk)+1,:,:);
            end
            fprintf('\n');
            H_i_full = permute(H_i_current_block,[1,2,5,3,6,4]);  % Size of final H_i [nRx nTx nRB 1 nInterf_eNodeB nTTI]
            
        end
    end
    fprintf('\n');
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
        this_bw,...                                 % System bandwidth
        config.channel_model.type,...               % Channel model used
        config.UE_speed*3.6,...                     % Channel speed (Km/h)
        date_time_string,...                        % Date string
        lab_ind,...                                 % Lab index (in order to avoid multiple labs writing to the same file)
        par_ID);                                    % ID of the current worker (for parallel simulations only)
end
save('-v7.3',fullfile(config.results_folder,results_filename),'config','networkPathlossMap','networkShadowFadingMap','UEs','eNodeBs','eNodeB_sites');
