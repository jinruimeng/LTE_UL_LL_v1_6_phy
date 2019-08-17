function [filtered_rand_grid_OTOI, filtered_rand_grid_LOS, filtered_rand_grid_NLOS, min_roi_position, max_roi_position] = ...
    calculate_correlations(UE_positions, data_res, config)
% Correlation of Large Scale Parameters with the method proposed in 3GPP TR
% 36.873
% 1) Generate a grid of i.i.d. random variables for each large scale
% parameter (7 in total)
% 2) Define the matrix of dstances between the points in the Region of
% Interest (ROI), then construct the filter impulse response 
% 3) Filter the grid of i.i.d. random variables with the 2-dimensional FIR filter
%
% (c) Fjolla Ademaj, Martin Taranetz, ITC 2016

% Calculate ROI for random matrix generation
% lower left corner
min_position = min(UE_positions);
max_position = max(UE_positions);
extra_distance = 100;      % Reference: Winner channel model, extra_distance refers to 2D value in the TR.36.873, used to make the grid larger                    
min_roi_position = min_position + [-extra_distance -extra_distance];
max_roi_position = max_position + [extra_distance extra_distance];
roi_maximum_pixels     = LTE_common_pos_to_pixel( max_roi_position, min_roi_position, data_res);
% Generate seven grids with iid normal random variables ~N(0,1); 
% one for each large scale parameter; 
normal_rand_grid_OTOI = normrnd(0,1, [roi_maximum_pixels(1), roi_maximum_pixels(2), 6]);
normal_rand_grid_LOS  = normrnd(0,1, [roi_maximum_pixels(1), roi_maximum_pixels(2), 7]);
normal_rand_grid_NLOS = normrnd(0,1, [roi_maximum_pixels(1), roi_maximum_pixels(2), 6]);

% Calculate distance grid for filtering.
filter_maximum_pixels    = LTE_common_pos_to_pixel([2*extra_distance, 2*extra_distance], [0 0], data_res);
center_position          = LTE_common_pos_to_pixel([2*extra_distance, 2*extra_distance], [extra_distance extra_distance], data_res);
x_positions = repmat((1:filter_maximum_pixels(1)) - center_position(1), filter_maximum_pixels(1),1);
distance_matrix    = sqrt((data_res * x_positions).^2 + (data_res * x_positions').^2);

% LSP correlation distances in horizotal plane Table 7.3-6
corr_distances_OTOI_vector = [config.SF_lambda_OTOI, config.DS_lambda_OTOI, config.AS_D_lambda_OTOI, config.AS_A_lambda_OTOI, config.ZS_D_lambda_OTOI, config.ZS_A_lambda_OTOI]; 
corr_distances_LOS_vector = [config.SF_lambda_LOS,config.KF_lambda_LOS, config.DS_lambda_LOS, config.AS_D_lambda_LOS, config.AS_A_lambda_LOS, config.ZS_D_lambda_LOS, config.ZS_A_lambda_LOS];
corr_distances_NLOS_vector = [config.SF_lambda_NLOS, config.DS_lambda_NLOS, config.AS_D_lambda_NLOS, config.AS_A_lambda_NLOS, config.ZS_D_lambda_NLOS, config.ZS_A_lambda_NLOS];

corr_distances_OTOI = reshape(corr_distances_OTOI_vector, 1,1,6); 
corr_distances_LOS = reshape(corr_distances_LOS_vector, 1,1,7);
corr_distances_NLOS = reshape(corr_distances_NLOS_vector, 1,1,6);

% One filter for each large scale parameter.
distance_matrices_OTOI          = repmat(distance_matrix, [1,1,6]); 
distance_matrices_LOS           = repmat(distance_matrix, [1,1,7]);
distance_matrices_NLOS          = repmat(distance_matrix, [1,1,6]);

corr_distances_matrix_OTOI      = repmat(corr_distances_OTOI, [filter_maximum_pixels(1), filter_maximum_pixels(2), 1]); 
corr_distances_matrix_LOS       = repmat(corr_distances_LOS, [filter_maximum_pixels(1), filter_maximum_pixels(2), 1]);
corr_distances_matrix_NLOS      = repmat(corr_distances_NLOS, [filter_maximum_pixels(1), filter_maximum_pixels(2), 1]);

% Filter impulse response
H_OTOI            = exp(-distance_matrices_OTOI./corr_distances_matrix_OTOI); 
H_LOS             = exp(-distance_matrices_LOS./corr_distances_matrix_LOS);
H_NLOS            = exp(-distance_matrices_NLOS./corr_distances_matrix_NLOS);

% Filter Gaussian grid in 2D to get LSP auto-correlation

for ii=1:length(corr_distances_OTOI_vector)
    % Apply 2-dimensional FIR filter
    filtered_grid_OTOI = filter2(H_OTOI(:,:,ii), normal_rand_grid_OTOI(:,:,ii));  
    % Extract only the inner part not affected by zero paded edges
    filtered_grid_OTOI_cropped = filtered_grid_OTOI((size(H_OTOI(:,:,ii),1)+1)/2:end-(size(H_OTOI(:,:,ii),1)-1)/2, (size(H_OTOI(:,:,ii),1)+1)/2:end-(size(H_OTOI(:,:,ii),1)-1)/2); 
    % Normalize: Scale the filtered grid by std of the cropped and filtered grid to obtain the same std as before the fitlering
    filtered_grid_OTOI_scaled = filtered_grid_OTOI/std(filtered_grid_OTOI_cropped(:)); 
    filtered_rand_grid_OTOI(:,:,ii) = filtered_grid_OTOI_scaled;
end

for ii=1:length(corr_distances_LOS_vector)
     % Apply 2-dimensional FIR filter
    filtered_grid_LOS = filter2(H_LOS(:,:,ii),normal_rand_grid_LOS(:,:,ii));
    % Extract only the inner part not affected by zero paded edges
    filtered_grid_LOS_cropped = filtered_grid_LOS((size(H_LOS(:,:,ii),1)+1)/2:end-(size(H_LOS(:,:,ii),1)-1)/2, (size(H_LOS(:,:,ii),1)+1)/2:end-(size(H_LOS(:,:,ii),1)-1)/2);
    % Normalize: Scale the filtered grid by std of the cropped and filtered grid to obtain the same std as before the fitlering
    filtered_grid_LOS_scaled= filtered_grid_LOS/std(filtered_grid_LOS_cropped(:));
    filtered_rand_grid_LOS(:,:,ii) = filtered_grid_LOS_scaled;
end

for ii=1:length(corr_distances_NLOS_vector)
    filtered_grid_NLOS = filter2(H_NLOS(:,:,ii),normal_rand_grid_NLOS(:,:,ii));
    filtered_grid_NLOS_cropped = filtered_grid_NLOS((size(H_NLOS(:,:,ii),1)+1)/2:end-(size(H_NLOS(:,:,ii),1)-1)/2, (size(H_NLOS(:,:,ii),1)+1)/2:end-(size(H_NLOS(:,:,ii),1)-1)/2);
    filtered_grid_NLOS_scaled = filtered_grid_NLOS/std(filtered_grid_NLOS_cropped(:));
    filtered_rand_grid_NLOS(:,:,ii) = filtered_grid_NLOS_scaled;

end
