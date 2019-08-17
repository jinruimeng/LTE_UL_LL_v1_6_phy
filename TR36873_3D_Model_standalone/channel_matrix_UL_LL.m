function [channel_matrices, pathlosses] = channel_matrix_UL_LL(N_subframes, params)
% Generates 3D model TR 36.873 channels for Vienna LTE-A Uplink Linklevel
% simulator

%% Parameter checks
if params.TR_36_873.LOS_according_model
    if isfield(params.TR_36_873, 'UE_is_LOS')
        warning('Ignoring explicitly defined LOS conditions because random setting is active');
    end
else
    if ~isfield(params.TR_36_873, 'UE_is_LOS')
        error('LOS conditions not defined');
    end
    if ~isequal(size(params.TR_36_873.UE_is_LOS), [params.nBS, params.nBS*params.nUE])
        error('LOS conditions variable has not size [params.nBS, params.nBS*params.nUE]');
    end
end

if params.TR_36_873.indoor_according_model
    if isfield(params.TR_36_873, 'UE_is_indoor')
        warning('Ignoring explicitly defined indoor states because random setting is active');
    end
    if isfield(params.TR_36_873, 'UE_dist_indoor')
        warning('Ignoring explicitly defined indoor distances because random setting is active');
    end
else
    if ~isfield(params.TR_36_873, 'UE_is_indoor')
        error('Indoor conditions not defined');
    end
    if ~isfield(params.TR_36_873, 'UE_dist_indoor')
        error('Indoor distances not defined');
    end    
    if ~isequal(size(params.TR_36_873.UE_is_indoor), [1, params.nBS*params.nUE])
        error('Indoor state variable has not size [1, params.nBS*params.nUE]');
    end
    if ~isequal(size(params.TR_36_873.UE_dist_indoor), [params.nBS, params.nBS*params.nUE])
        error('Indoor distance variable has not size [params.nBS, params.nBS*params.nUE]');
    end
end

if params.TR_36_873.eNodeB_hex_grid
    if isfield(params.TR_36_873, 'eNodeB_pos')
        warning('Ignoring explicitly defined eNodeB positions because hex grid option is active');
    end
else
    if ~isfield(params.TR_36_873, 'eNodeB_pos')
        error('eNodeB positions not defined');
    end
    if ~isequal(size(params.TR_36_873.eNodeB_pos), [params.nBS, 2])
        error('eNodeB positions variable has not size [params.nBS, 2]');
    end
end

if params.TR_36_873.UE_pos_rand
    if isfield(params.TR_36_873, 'UE_pos')
        warning('Ignoring explicitly defined UE positions because random pacement option is active');
    end
else
    if ~isfield(params.TR_36_873, 'UE_pos')
        error('UE positions not defined');
    end
    if ~isequal(size(params.TR_36_873.UE_pos), [params.nBS*params.nUE, 2])
        error('UE positions variable has not size [params.nBS*params.nUE, 2]');
    end
end

if params.TR_36_873.UE_heights_according_to_model
    if isfield(params.TR_36_873, 'UE_heights')
        warning('Ignoring explicitly defined UE heights because random placement is active');
    end
else
    if ~isfield(params.TR_36_873, 'UE_heights')
        error('UE heights not defined');
    end
    if ~isequal(size(params.TR_36_873.UE_heights), [1, params.nBS*params.nUE])
        error('UE height variable has not size [1, params.nBS*params.nUE]');
    end
end

if params.TR_36_873.eNodeB_heights_according_to_model
    if isfield(params.TR_36_873, 'eNodeB_heights')
    	warning('Ignoring explicitly defined eNodeB heights because setting according to model is active');
    end
else
    if ~isfield(params.TR_36_873, 'eNodeB_heights')
    	error('eNodeB heights not defined');
    end
    if ~isequal(size(params.TR_36_873.eNodeB_heights), [1, params.nBS])
        error('eNodeB height variable has not size [1, params.nBS]');
    end
end

if ~isequal(size(params.TR_36_873.antenna_azimuth_offset), [1, params.nBS])
    error('Antenna azimuth offset variable has not size [1, params.nBS]');
end

if ~ismember(params.TR_36_873.antenna.antenna_gain_pattern, {'TR36.873 3D antenna omnidirectional', 'TR36.873 3D antenna'})
    error('Specified antenna gain pattern not supported');
end

%% Basic settings
global DEBUG_LEVEL;
config.generate_uplink = true; % true means uplink, false means downlink
config.frequency = params.carrier_freq_UP;
config.wavelength = params.speed_of_light/config.frequency;

config.debug_level = DEBUG_LEVEL;
config.bandwidth = params.Bandwidth;
config.CP_length = params.CyclicPrefix;

switch params.ChanMod_config.filtering
    case 'BlockFading'
        config.BF_length = params.Tsubframe;
        config.simulation_time_blocks = N_subframes;
        config.fading_interval_name = 'Subframe';
    case 'FastFading'
        config.BF_length = params.Tsubframe/params.Nsub;
        config.simulation_time_blocks = params.Nsub*N_subframes;
        config.fading_interval_name = 'Symbol';
    otherwise
        error('Channel filtering type invalid');
end

% General user and eNodeB settings
config.UE_per_eNodeB = params.nUE;
config.number_eNodeBs = params.nBS;
config.UE_speed = params.UE_config.user_speed;
config.eNodeB_nTX = params.BS_config.nRX;
config.UE_nRX = params.UE_config.nTX;

% Fading models
if strcmp(params.TR_36_873.environment, 'UMa')
    config.channel_model.type = '3D_UMa_fading';
    config.macroscopic_pathloss_model_settings.environment = '3D_UMa';
    config.inter_site_dist = 500; % meters
elseif strcmp(params.TR_36_873.environment, 'UMi')
    config.channel_model.type = '3D_UMi_fading';
    config.macroscopic_pathloss_model_settings.environment = '3D_UMi';
    config.inter_site_dist = 200; % meters
else
    error('Environment not supported');
end

config.macroscopic_pathloss_model = 'TR 36.873';
config.shadow_fading_type = 'none';
config.map_resolution = params.TR_36_873.map_resolution;

% Determine eNodeBs azimuth and height
config.antenna_azimuth_offset = params.TR_36_873.antenna_azimuth_offset;

if params.TR_36_873.eNodeB_heights_according_to_model
    if strcmp(params.TR_36_873.environment, 'UMa')
        config.eNodeB_heights = repmat(25, [1, params.nBS]);
    elseif strcmp(params.TR_36_873.environment, 'UMi')
        config.eNodeB_heights = repmat(10, [1, params.nBS]);
    end
else
    config.eNodeB_heights = params.TR_36_873.eNodeB_heights;
end

% Determine UE antenna parameters
config.UE_antenna_polarization = params.TR_36_873.UE_antenna_polarization;
if strcmp(config.UE_antenna_polarization, 'XPOL')
    config.UE_antenna_slant_angle = 90;
elseif strcmp(config.UE_antenna_polarization, 'ULA')
    config.UE_antenna_slant_angle = 0;
else
    error('Antenna polarization not supported');
end
config.UE_antenna_element_horizontal_spacing = params.TR_36_873.UE_antenna_element_horizontal_spacing * config.wavelength;

% Determine eNodeB antenna parameters
config.antenna.antenna_gain_pattern = params.TR_36_873.antenna.antenna_gain_pattern;
config.antenna.antenna_polarization = params.TR_36_873.antenna.antenna_polarization;
if strcmp(config.antenna.antenna_polarization, 'XPOL')
    config.antenna.slant_angle = 45;
elseif strcmp(config.antenna.antenna_polarization, 'ULA')
    config.antenna.slant_angle = 0;
else
    error('Antenna polarization not supported');
end
config.antenna_element_vertical_spacing = params.TR_36_873.antenna_element_vertical_spacing * config.wavelength;
config.antenna_element_horizontal_spacing = params.TR_36_873.antenna_element_horizontal_spacing * config.wavelength;
config.nr_of_antenna_elements_in_each_column = params.TR_36_873.nr_of_antenna_elements_in_each_column;
config.electrical_downtilt = params.TR_36_873.electrical_downtilt;
config.mechanical_downtilt = params.TR_36_873.mechanical_downtilt;
config.mechanical_slant = params.TR_36_873.mechanical_slant;
config.antenna.max_antenna_gain = params.TR_36_873.antenna.max_antenna_gain;

% Load additional parameters
config = load_specific_params(config);

%% generate eNodeBs and eNodeB_sites
if config.debug_level > 0
    fprintf('TR 36.873 - 3D Channel Model\n');
    fprintf('\tGenerating eNodeBs\n');
end

% Determine the eNodeBs positions
if params.TR_36_873.eNodeB_hex_grid
    n_eNodeB_rings = ceil(-0.5+sqrt((4*params.nBS-1)/12)) + (params.nBS == 1);
    eNodeB_site_positions = linklevel.network_geometry.hexagonal_eNodeB_grid(config.inter_site_dist, n_eNodeB_rings);
    eNodeB_site_positions = eNodeB_site_positions(1:params.nBS, :);
else
    eNodeB_site_positions = params.TR_36_873.eNodeB_pos;
end
eNodeB_sites = network_elements.eNodeB_site;
eNodeBs = network_elements.eNodeB;

% Set most of the eNodeB parameters
for bb = 1:params.nBS
    % Site settings
    eNodeB_sites(bb).id = bb;
    eNodeB_sites(bb).pos = eNodeB_site_positions(bb,:);
    eNodeB_sites(bb).site_type = 'macro';
    % Sector settings
    eNodeBs(bb) = network_elements.eNodeB;
    eNodeBs(bb).id = 1;
    eNodeBs(bb).eNodeB_id = bb;
    eNodeBs(bb).parent_eNodeB = eNodeB_sites(bb);
    eNodeBs(bb).pos = eNodeB_site_positions(bb,:);
    eNodeBs(bb).nTX = config.eNodeB_nTX;
    eNodeBs(bb).tx_height = config.eNodeB_heights(bb);
    eNodeBs(bb).boresight = config.antenna_azimuth_offset(bb);
    eNodeBs(bb).antenna_type  = config.antenna.antenna_gain_pattern;
    % NOTE: needed so antenna instantiation does not fail
    config.tx_height = eNodeBs(bb).tx_height;
    eNodeBs(bb) = antennas.antenna.attach_antenna_to_eNodeB(eNodeBs(bb),config);
    eNodeB_sites(bb).sectors = eNodeBs(bb);
end

% Initialize clock
networkClock = network_elements.clock(config.BF_length);
for bb = 1:params.nBS
    eNodeBs(bb).clock = networkClock;
    eNodeBs(bb).neighbors_eNodeB = eNodeBs([1:(bb-1) (bb+1):length(eNodeBs)]);
end

%% Generate UEs
if config.debug_level > 0
    fprintf('\tGenerating UEs\n');
end
UEs = network_elements.UE;

% Determine the UE positions
if params.TR_36_873.UE_pos_rand
    UE_positions = zeros(params.nBS*params.nUE, 2);
    for bb = 1:params.nBS
        UE_positions(1+(bb-1)*params.nUE : params.nUE+(bb-1)*params.nUE, :) = ...
            linklevel.network_geometry.rand_circle(params.nUE,...
            eNodeBs(bb).pos(1), eNodeBs(bb).pos(2), config.inter_site_dist);
    end
else
    UE_positions = params.TR_36_873.UE_pos;
end

%% Determine UE indoor state and distance
if params.TR_36_873.indoor_according_model
    UE_is_indoor = rand(1, params.nBS*params.nUE) < 0.8;
    UE_dist_indoor = 25 * rand(params.nBS, params.nBS*params.nUE, 1);
else
    UE_is_indoor = params.TR_36_873.UE_is_indoor;
    UE_dist_indoor = params.TR_36_873.UE_dist_indoor;
end

% Determine UE height
if params.TR_36_873.UE_heights_according_to_model
    config.UE_heights = linklevel.network_geometry.UE_heights(params.nBS*params.nUE, UE_is_indoor);
else
    config.UE_heights = params.TR_36_873.UE_heights;
end

% Set most UE parameters
for uu = 1:params.nUE*params.nBS
    UEs(uu)= network_elements.UE;
    UEs(uu).id = uu;
    UEs(uu).pos = UE_positions(uu,:);
    UEs(uu).rx_height = config.UE_heights(uu);
    UEs(uu).nRX = config.UE_nRX;
    UEs(uu).direction = floor(random('unif',0,359));
    UEs(uu).H_0_channel_trace = [];
    UEs(uu).sampled_channel_H_0 = [];
    UEs(uu).H_i_channel_trace = [];
    UEs(uu).sampled_channel_H_i = [];
    UEs(uu).H_i_after_fft = [];
    UEs(uu).H_0_final = [];
    UEs(uu).H_i_full_final =[];
    UEs(uu).TTI_of_smallscale_fading_recalculation =[];
    UEs(uu).recalculate_3D_smallscale_fading = true;
    UEs(uu).channels = channel_gain_wrappers.TR36873_Fading_3D_Channel_LL(config);
end

%% Determine LOS state
if params.TR_36_873.LOS_according_model
    UE_height = zeros(1, params.nBS*params.nUE);
    UE_BS_dist_out = zeros(params.nBS, params.nBS*params.nUE);
    for uu = 1:params.nBS*params.nUE
        UE_height(uu) = UEs(uu).rx_height;
        for bb = 1:params.nBS
            UE_BS_dist_out(bb, uu) = norm(eNodeB_sites(bb).pos - UEs(uu).pos) - UE_is_indoor(uu)*UE_dist_indoor(bb,uu);
        end
    end
    UE_height = repmat(UE_height, [params.nBS, 1]);
    if strcmp(config.channel_model.type, '3D_UMa_fading')
        LOS_probabilities = linklevel.propagation.LOS_probability_UMa(UE_BS_dist_out, UE_height);
    elseif strcmp(config.channel_model.type, '3D_UMi_fading')
        LOS_probabilities = linklevel.propagation.LOS_probability_UMi(UE_BS_dist_out, UE_height);
    end
    UE_is_LOS = rand(params.nBS, params.nBS*params.nUE) < LOS_probabilities;
else
    UE_is_LOS = params.TR_36_873.UE_is_LOS;
end

% Set LOS state for desired connections
for uu = 1:params.nBS*params.nUE
    bb = find(params.connection_table(:, uu) == true);
    UEs(uu).is_indoor = UE_is_indoor(uu);
    UEs(uu).is_LOS(1) = UE_is_LOS(bb, uu);
    UEs(uu).dist_indoor = UE_dist_indoor(bb, uu);
end

% Set LOS state for interfering connections
for uu = 1:params.nBS*params.nUE  
    bb_interf = 1;
    for bb = 1:params.nBS
        if ~params.connection_table(bb, uu)
            UEs(uu).is_LOS(bb_interf + 1) = UE_is_LOS(bb, uu);
            bb_interf = bb_interf + 1;
        end
    end
end

% Attach UEs to eNodeBs according to connection table
for uu = 1:params.nBS*params.nUE
    bb = params.connection_table(:, uu) == true;
    eNodeBs(bb).attachUser(UEs(uu));
end

%% Generate pathlosses
if isfield(params, 'TR_36_873.pathloss_enabled') && params.TR_36_873.pathloss_enabled == true
    pathlosses = NaN(params.nBS, params.nBS*params.nUE);
    if strcmp(params.TR_36_873.environment, 'UMa')
        for bb = 1:params.nBS
            for uu = 1:params.nBS*params.nUE
                d2D = norm(UEs(uu).pos - eNodeB_sites(bb).pos);
                if params.TR_36_873.UE_is_indoor(uu)
                    d_2D_in = UE_dist_indoor(bb, uu);
                    is_LOS = UE_is_LOS(bb, uu);
                    pathlosses(bb, uu) = linklevel.pathloss.pathloss_UMa_OtoI(config.frequency, eNodeBs(bb).tx_height, UEs(uu).rx_height,...
                        d2D, d_2D_in, is_LOS);
                elseif params.TR_36_873.UE_is_LOS(bb, uu)
                    pathlosses(bb, uu) = linklevel.pathloss.pathloss_UMa_LOS(config.frequency, eNodeBs(bb).tx_height, UEs(uu).rx_height, d2D);
                else
                    pathlosses(bb, uu) = linklevel.pathloss.pathloss_UMa_NLOS(config.frequency, eNodeBs(bb).tx_height, UEs(uu).rx_height, d2D);
                end
            end
        end
    elseif strcmp(params.TR_36_873.environment, 'UMi')
        for bb = 1:params.nBS
            for uu = 1:params.nBS*params.nUE
                d2D = norm(UEs(uu).pos - eNodeB_sites(bb).pos);
                if params.TR_36_873.UE_is_indoor(uu)
                    d_2D_in = UE_dist_indoor(bb, uu);
                    is_LOS = UE_is_LOS(bb, uu);
                    pathlosses(bb, uu) = linklevel.pathloss.pathloss_UMi_OtoI(config.frequency, eNodeBs(bb).tx_height, UEs(uu).rx_height,...
                        d2D, d_2D_in, is_LOS);
                elseif params.TR_36_873.UE_is_LOS(bb, uu)
                    pathlosses(bb, uu) = linklevel.pathloss.pathloss_UMi_LOS(config.frequency, eNodeBs(bb).tx_height, UEs(uu).rx_height, d2D);
                else
                    pathlosses(bb, uu) = linklevel.pathloss.pathloss_UMi_NLOS(config.frequency, eNodeBs(bb).tx_height, UEs(uu).rx_height, d2D);
                end
            end
        end
    end
else
    pathlosses = NaN;
end

%% Start main simulation loop: calculate correlated LSPs and generate channel coefficients
while networkClock.current_block < config.simulation_time_blocks
    % Advance the network clock
    networkClock.advance_1_TTI;
    if config.debug_level > 0
        fprintf('\t%s %d/%d:\n', config.fading_interval_name, networkClock.current_block, config.simulation_time_blocks);
    end
    
    for bb = 1:params.nBS
        %% Generate correlated Large Scale Parameters
        if networkClock.current_block == 1
            if config.debug_level > 0
                fprintf('\t\tLarge scale parameter generation: eNodeB %d\n', eNodeBs(bb).eNodeB_id);
            end
            UE_positions = zeros(eNodeBs(bb).attached_UEs, 2);
            data_res = config.map_resolution;
            % Get UEs that are attached to this eNodeB and check whether they are in LOS or NLOS
            % Generate LSPs for interfering links
            if eNodeBs(bb).attached_UEs > 0
                for jj=1:eNodeBs(bb).attached_UEs
                    UE_distance_to_site = zeros(length(eNodeBs(bb).neighbors_eNodeB)+1,1);
                    eNodeBs(bb).attached_UEs_vector(jj).all_large_scale_params = zeros(length(eNodeBs(bb).neighbors_eNodeB)+1,7);
                    UE_positions(jj,:) = eNodeBs(bb).attached_UEs_vector(jj).pos;
                    eNodeBs(bb).attached_UEs_vector(jj).recalculate_3D_smallscale_fading = true; % In case the large scale parameters have changed, the fast scale fading is recalculated.
                    if ~isempty(eNodeBs(bb).neighbors_eNodeB)
                        for ii= 1:length(eNodeBs(bb).neighbors_eNodeB)
                            eNodeBs(bb).attached_UEs_vector(jj).recalculate_3D_smallscale_fading_i(ii) = true;
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
            % If only one UE is attached to the eNodeB than no cross-correlation is assumed
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
        if config.debug_level > 0
            fprintf('\t\tChannel generation: eNodeB %d, UE', bb);
        end
        for uu = 1:eNodeBs(bb).attached_UEs
            if config.debug_level > 0
                fprintf(' %d', eNodeBs(bb).attached_UEs_vector(uu).id);
            end
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
                
                % Channel matrix after FFT with size [nRx nTx 1 nTTI nRB]
                H_0_after_fft = eNodeBs(bb).attached_UEs_vector(uu).channels.get_RB_trace(eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_0);
                eNodeBs(bb).attached_UEs_vector(uu).H_0_final = H_0_after_fft;
            end
            
            % INTERFERING CHANNELS
            for kk = 1:length(eNodeBs(bb).neighbors_eNodeB) % indicates the interfering eNodeB index
                if eNodeBs(bb).attached_UEs_vector(uu).recalculate_3D_smallscale_fading_i(kk)
                    eNodeBs(bb).attached_UEs_vector(uu).sampled_channel_H_i = zeros(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).neighbors_eNodeB(kk).nTX, 1, config.simulation_time_blocks-networkClock.current_block+1, 200);                   
                    eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft = zeros(eNodeBs(bb).attached_UEs_vector(uu).nRX, eNodeBs(bb).neighbors_eNodeB(kk).nTX, 1, config.simulation_time_blocks-networkClock.current_block+1, params.Nsc*config.N_RB , length(eNodeBs(bb).neighbors_eNodeB));
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
                    eNodeBs(bb).attached_UEs_vector(uu).H_i_full_final(:,:,:,:,:,kk) = eNodeBs(bb).attached_UEs_vector(uu).H_i_after_fft(:,:,:,:,:,kk);
                end
            end
        end
        if config.debug_level > 0
            fprintf('\n');
        end
    end
end

% Fill up the cell of channel matrices that will be returned
channel_matrices = cell(params.nBS, params.nBS*params.nUE);
for uu = 1:params.nBS*params.nUE
    bb_interf = 1;
    for bb = 1:params.nBS
        if params.connection_table(bb, uu) == true
            % Desired channels
            channel_matrix = UEs(uu).H_0_final;
        else
            % Interferers channels
            channel_matrix = UEs(uu).H_i_full_final(:, :, :, :, :, bb_interf);
            bb_interf = bb_interf + 1;
        end
        % NOTE: Includes swap of TX <-> RX to get uplink channel
        channel_matrices{bb, uu} = permute(channel_matrix, [5, 3, 2, 1, 4]);
        switch params.ChanMod_config.filtering
            case 'BlockFading'
                channel_matrices{bb, uu} = repmat(channel_matrices{bb, uu},...
                    [1, params.Nsub, 1, 1, 1]);
            case 'FastFading'
                size_H = size(channel_matrices{bb, uu});
                channel_matrices{bb, uu} = reshape(channel_matrices{bb, uu},...
                    [size_H(1), params.Nsub, size_H(3), size_H(4), N_subframes]);
        end
    end
end

if config.debug_level > 0
    fprintf('\n');
end

end