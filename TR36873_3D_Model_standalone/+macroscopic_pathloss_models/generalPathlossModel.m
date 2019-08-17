classdef generalPathlossModel
    % Abstract class that represents a pathloss model. Also wraps the code
    % that calculates the macroscopic pathloss maps 
    % (c) Josep Colom Ikuno, INTHFT, 2008
    % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016
    
    properties
        % this model's name
        name
    end
    
    methods (Abstract)
        pathloss_in_db = pathloss(distance)
    end
    
    methods (Static)
        function calculate_pathloss_maps(config,eNodeBs,networkMacroscopicPathlossMap,varargin)
            % Calculates the pathloss map for each eNodeb sector
            % 
            if ~isempty(varargin)
                elevation_map     = varargin{1};
                elevation_map_set = true;
            else
                elevation_map     = 0;
                elevation_map_set = false;
            end
            
            % Calculates the pathloss maps for a given eNodeB set (cell set)
            if isprop(eNodeBs(1),'sectors')
                total_sectors = length([eNodeBs.sectors]);
            else
                total_sectors = length(eNodeBs);
            end
            
            data_res               = networkMacroscopicPathlossMap.data_res;
            roi_x                  = networkMacroscopicPathlossMap.roi_x;
            roi_y                  = networkMacroscopicPathlossMap.roi_y;
            roi_maximum_pixels     = LTE_common_pos_to_pixel( [roi_x(2) roi_y(2)], [roi_x(1) roi_y(1)], data_res);
            roi_height_pixels      = roi_maximum_pixels(2);
            roi_width_pixels       = roi_maximum_pixels(1);
            distance_matrix_2D     = zeros(roi_height_pixels,roi_width_pixels,length(eNodeBs));
            LOS_positions          = zeros(roi_height_pixels,roi_width_pixels,length(eNodeBs)); % LOS map.
            cell_pathloss_data_dB  = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            sector_antenna_gain_dB = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            sector_distances       = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            
            site_positions     = reshape([eNodeBs.pos],2,[])';
            site_positions_pix = LTE_common_pos_to_pixel( site_positions, [roi_x(1) roi_y(1)], data_res);
            
            % Generate distance and angle matrix
            position_grid_pixels      = zeros(roi_height_pixels*roi_width_pixels,2);
            position_grid_pixels(:,1) = reshape(repmat(1:roi_width_pixels,roi_height_pixels,1),1,roi_width_pixels*roi_height_pixels);
            position_grid_pixels(:,2) = repmat(1:roi_height_pixels,1,roi_width_pixels);
            position_grid_meters      = LTE_common_pixel_to_pos(...
                position_grid_pixels,...
                networkMacroscopicPathlossMap.coordinate_origin,...
                networkMacroscopicPathlossMap.data_res);
            
            
            %% Sector pathloss
            s_idx = 1;
            eNodeB_site_type = cell(1,total_sectors);
            all_sectors   = [eNodeBs.sectors];
            eNodeB_id_set = [all_sectors.eNodeB_id];
            
            for b_ = 1:length(eNodeBs)
                distances = sqrt(...
                    (position_grid_meters(:,1)-eNodeBs(b_).pos(1)).^2 + ...
                    (position_grid_meters(:,2)-eNodeBs(b_).pos(2)).^2);
                distance_matrix_2D(:,:,b_) = reshape(distances,roi_height_pixels,roi_width_pixels);
                
                current_site_sectors = eNodeBs(b_).sectors;
                
                LOS_positions(:,:,b_) = current_site_sectors(1).macroscopic_pathloss_model.generate_LOS_positions(distance_matrix_2D(:,:,b_), networkMacroscopicPathlossMap);
                
                for s_ = 1:length(current_site_sectors)
                    
                    % Distance matrix for each sector
                    sector_distances(:,:,s_idx) = distance_matrix_2D(:,:,b_);
                    
                    % Calculate macroscopic pathloss using the macroscopic pathloss model from each eNodeB
                    % - The output of the pathloss model TR 36.873 is in dB
                    % Calculate path loss and LOS only for first
                    % sector, otherwise copy.
                    cell_pathloss_data_dB(:,:,s_idx)  = current_site_sectors(s_).macroscopic_pathloss_model.pathloss(sector_distances(:,:,s_idx),config,networkMacroscopicPathlossMap, LOS_positions(:,:,b_));
                    
                    
                    % Horizontal angle grid: Convert the azimuth
                    % (0deg=North, 90deg=East, 180^=South, 270deg=West) degrees to cartesian
                    angle_grid = (180/pi)*(...
                        atan2(...
                        (position_grid_meters(:,2)-eNodeBs(b_).pos(2)),...
                        (position_grid_meters(:,1)-eNodeBs(b_).pos(1)))) - ...
                        utils.miscUtils.wrapTo359(current_site_sectors(s_).boresight);
                    
                    % Horizontal angle grid
                    horizontal_angle_grid   = reshape(angle_grid,roi_height_pixels,roi_width_pixels);
                    horizontal_angle_grid_s = utils.miscUtils.wrapTo359(horizontal_angle_grid);
                    
                    % For Tr 36.873 3D channel model antenna pattern
                    % Vertical angle grid
                    % (0deg--> zenith, 90deg--> perpendicular to the antenna) degrees to cartesian
                    distance_3D = sqrt(distance_matrix_2D(:,:,b_).^2 +(current_site_sectors(s_).tx_height - 1.5).^2);
                    theta_arrival  = acosd((current_site_sectors(s_).tx_height - config.rx_height)./distance_3D);
                    vertical_angle_grid_el       = 180 - theta_arrival;
                    
                    % Calculate antenna gain from 'TR36873_3DAntenna'
                    % Set phi to (-180,180)
                    horizontal_angle_grid_s = horizontal_angle_grid_s + 180;
                    horizontal_angle_grid_s = mod(horizontal_angle_grid_s,360);
                    horizontal_angle_grid_s = horizontal_angle_grid_s - 180;
%                     
%                        vertical_angle_grid_el = vertical_angle_grid_el + 180;
%                                 vertical_angle_grid_el = mod(vertical_angle_grid_el,360);
%                                 vertical_angle_grid_el = vertical_angle_grid_el - 180;
                    

                       sector_antenna_gain_dB(:,:,s_idx) = current_site_sectors(s_).antenna.calculate_array_single_column_field_pattern(...
                           config,...
                        vertical_angle_grid_el,...
                        horizontal_angle_grid_s,...
                        current_site_sectors(s_).pos);
                    
%                     sector_antenna_gain_dB(:,:,s_idx) = current_site_sectors(s_).antenna.gain(...
%                         vertical_angle_grid_el,...
%                         horizontal_angle_grid_s);
                    
                    % figure(); imagesc(sector_antenna_gain_dB(:,:,s_idx))
                    
                    % Mapping between s_idx and b_/s_ pair
                    current_eNodeB_id = current_site_sectors(s_).eNodeB_id;
                    networkMacroscopicPathlossMap.sector_idx_mapping(current_eNodeB_id,:) = [b_ s_];
                    networkMacroscopicPathlossMap.site_sector_mapping(b_,s_)              = current_eNodeB_id;
                    
                    % Site type, for which path loss map is generated
                    eNodeB_site_type{s_idx}  = eNodeBs(b_).site_type;
                    
                    s_idx = s_idx + 1;
                end
            end
            
            % Fill in pathloss data
            cell_pathloss_data_dB(isnan(cell_pathloss_data_dB) | (cell_pathloss_data_dB<0)) = 0;
            networkMacroscopicPathlossMap.pathloss(:,:,eNodeB_id_set)  = 10.^((cell_pathloss_data_dB - sector_antenna_gain_dB)/10);
            networkMacroscopicPathlossMap.distances(:,:,eNodeB_id_set) = sector_distances;
            
            if ~iscell(networkMacroscopicPathlossMap.site_type)
                networkMacroscopicPathlossMap.site_type                = cell(length(eNodeB_id_set),1);
                networkMacroscopicPathlossMap.site_type(eNodeB_id_set) = eNodeB_site_type;
            else
                networkMacroscopicPathlossMap.site_type(eNodeB_id_set) = eNodeB_site_type;
            end
            
            % In case of TR 36.873, also store LOS map.
            networkMacroscopicPathlossMap.LOS_map                                = LOS_positions;
            networkMacroscopicPathlossMap.excluding_region_near_BS               = sector_distances > config.min_UE_eNodeB_distance; 
            networkMacroscopicPathlossMap.sector_antenna_gain(:,:,eNodeB_id_set) = sector_antenna_gain_dB;
            
        end
        
        function UE_indoor_map = calculate_UE_indoor_map(config, networkMacroscopicPathlossMap)
            % Generates the indoor/outdoor map pixel-wise
            %             roi_height_pixels         ... size of the ROI in pixels in y axis
            %             roi_width_pixels          ... size of the ROI in pixels in x axis
            % input:      config.indoor_UE_fraction ... indoor/outdoor UE fraction
            %
            % output:     UE_indoor_map             ... matrix of boolean type with the size of ROI
            %
            % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016

            data_res               = networkMacroscopicPathlossMap.data_res;
            roi_x                  = networkMacroscopicPathlossMap.roi_x;
            roi_y                  = networkMacroscopicPathlossMap.roi_y;
            roi_maximum_pixels     = LTE_common_pos_to_pixel( [roi_x(2) roi_y(2)], [roi_x(1) roi_y(1)], data_res);
            roi_height_pixels      = roi_maximum_pixels(2);
            roi_width_pixels       = roi_maximum_pixels(1);
            UE_indoor_map          = binornd(1,config.indoor_UE_fraction,roi_height_pixels, roi_width_pixels);
        end
        
        function UE_height_map = calculate_UE_height_map(config,networkMacroscopicPathlossMap)
            % Generates the UE height map pixel-wise
            %             roi_height_pixels         ... size of the ROI in pixels in y axis
            %             roi_width_pixels          ... size of the ROI in pixels in x axis
            %             building_floor_height     ... height of the
            %             building given with the number of floors
            %             UE_indoor_standing_floor  ... standing floor of the UE
            % input:      config.min_floor_number   ... minimum number of floors per building
            %             config.max_floor_number   ... maximum number of floors per building
            %
            % output:     UE_indoor_map             ... height in [m] 
            %
            % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016
            data_res               = networkMacroscopicPathlossMap.data_res;
            roi_x                  = networkMacroscopicPathlossMap.roi_x;
            roi_y                  = networkMacroscopicPathlossMap.roi_y;
            roi_maximum_pixels     = LTE_common_pos_to_pixel( [roi_x(2) roi_y(2)], [roi_x(1) roi_y(1)], data_res);
            roi_height_pixels      = roi_maximum_pixels(2);
            roi_width_pixels       = roi_maximum_pixels(1);
            building_floor_height     = randi([config.min_floor_number,config.max_floor_number]);
            UE_indoor_standing_floor  = randi(building_floor_height,roi_height_pixels,roi_width_pixels);
            UE_height_map          = 3*(UE_indoor_standing_floor - 1) + 1.5;
            UE_height_map          = networkMacroscopicPathlossMap.UE_indoor_map.*UE_height_map + ~networkMacroscopicPathlossMap.UE_indoor_map.*1.5; 
        end
        
        function UE_indoor_dist_map = calculate_UE_indoor_distance(networkMacroscopicPathlossMap)
            % Generates the indoor distance in [m]
            % used in the O-to-I scenario (Outdoor-to-Indoor)
            % Determines the distance form the UE position to the outer
            % wall of the building where the Ue is located.
            %             roi_height_pixels         ... size of the ROI in pixels in y axis
            %             roi_width_pixels          ... size of the ROI in pixels in x axis
            %
            % output:     UE_indoor_dist_map        ... indoor distance 
            %
            % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016
            data_res               = networkMacroscopicPathlossMap.data_res;
            roi_x                  = networkMacroscopicPathlossMap.roi_x;
            roi_y                  = networkMacroscopicPathlossMap.roi_y;
            roi_maximum_pixels     = LTE_common_pos_to_pixel( [roi_x(2) roi_y(2)], [roi_x(1) roi_y(1)], data_res);
            roi_height_pixels      = roi_maximum_pixels(2);
            roi_width_pixels       = roi_maximum_pixels(1);
            indoor_distance = 25.*rand(roi_height_pixels,roi_width_pixels);
            UE_indoor_dist_map = networkMacroscopicPathlossMap.UE_indoor_map.*indoor_distance;
        end

        function macroscopic_pathloss_model = generateMacroscopicPathlossModel(LTE_config,macroscopic_pathloss_model_name,frequency,macroscopic_pathloss_model_settings)
            % Returns an appropriate pathloss model based on the provided information
            print_output = true;
            switch macroscopic_pathloss_model_name
                case 'TR 36.873'
                    macroscopic_pathloss_model = macroscopic_pathloss_models.TR36873PathlossModel(frequency,macroscopic_pathloss_model_settings.environment);
                    if print_output && LTE_config.debug_level>=1
                        fprintf('TR 36.873-recommended pathloss model, %s environment\n',macroscopic_pathloss_model_settings.environment);
                    end     
                otherwise
                    error('"%s" macroscopic pathloss model not supported',macroscopic_pathloss_model_name);
            end
        end
    end
end
