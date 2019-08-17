% Calculation of Large scale parameters according to Circular Angle
% Spread method proposed in 3GPP TR 25.996
% Load result files saved in .\resultsFolder. For calibration it is
% recommended to do at least 10 simulation runs

clear all; close all; clc;
% Load all results files from folder
foldername = './resultsFolder';
filelist  = dir(sprintf('%s/*.mat',foldername));


cluster_NLOS = 20; % for UMa -->20; for UMi --> 19
cluster_OTOI = 12; %12
cluster_LOS = 12; %12

zenith_dep_NLOS_all   = [];
zenith_arr_NLOS_all    = [];
cluster_power_NLOS_all = [];

zenith_dep_OTOI_all = [];
zenith_arr_OTOI_all= [];
cluster_power_OTOI_all= [];

zenith_dep_LOS_all= [];
zenith_arr_LOS_all= [];
cluster_power_LOS_all= [];

% extract results from each results file.
for ii = 1:length(filelist)
    filename = filelist(ii).name;
    load(fullfile(foldername,filename),'eNodeBs','UEs');
    fprintf(filelist(ii).name); fprintf('\n');
    % Filter UEs that were simulated.
    active_UE_idcs    = ~[UEs(:).deactivate_UE];
    
    extract_UE_data = reshape([eNodeBs(:,:).attached_UEs_vector],length(eNodeBs),[]);
    each_UE_data = reshape([extract_UE_data(:,:).is_indoor],length(eNodeBs(1,1).attached_UEs_vector),[]);
    OTOI_ind = logical([each_UE_data]);
    each_UE_data2 = reshape([extract_UE_data(:,:).is_LOS],length(eNodeBs(1,1).attached_UEs_vector),[]);
    LOS_ind  = logical(~OTOI_ind.*logical([each_UE_data2]));
    NLOS_ind = logical(~OTOI_ind.*~logical([each_UE_data2]));

    zenith_dep_NLOS = reshape([extract_UE_data(NLOS_ind).zenith_angles_of_departure_all],20,cluster_NLOS,[]); 
    zenith_arr_NLOS = reshape([extract_UE_data(NLOS_ind).zenith_angles_of_arrival_all],20,cluster_NLOS,[]); 
    cluster_power_NLOS = reshape([extract_UE_data(NLOS_ind).cluster_power_per_ray],20,cluster_NLOS,[]);
    
    zenith_dep_OTOI = reshape([extract_UE_data(OTOI_ind).zenith_angles_of_departure_all],20,cluster_OTOI,[]); 
    zenith_arr_OTOI = reshape([extract_UE_data(OTOI_ind).zenith_angles_of_arrival_all],20,cluster_OTOI,[]); 
    cluster_power_OTOI = reshape([extract_UE_data(OTOI_ind).cluster_power_per_ray],20,cluster_OTOI,[]);
    
    zenith_dep_LOS = reshape([extract_UE_data(LOS_ind).zenith_angles_of_departure_all],20,cluster_LOS,[]); 
    zenith_arr_LOS = reshape([extract_UE_data(LOS_ind).zenith_angles_of_arrival_all],20,cluster_LOS,[]); 
    cluster_power_LOS = reshape([extract_UE_data(LOS_ind).cluster_power_per_ray],20,cluster_LOS,[]);
    
    zenith_dep_NLOS_all =cat(3,zenith_dep_NLOS_all, zenith_dep_NLOS);
    zenith_arr_NLOS_all = cat(3,zenith_arr_NLOS_all, zenith_arr_NLOS);
    cluster_power_NLOS_all = cat(3,cluster_power_NLOS_all, cluster_power_NLOS);
    
    zenith_dep_OTOI_all = cat(3, zenith_dep_OTOI_all, zenith_dep_OTOI);
    zenith_arr_OTOI_all = cat(3, zenith_arr_OTOI_all, zenith_arr_OTOI);
    cluster_power_OTOI_all = cat(3, cluster_power_OTOI_all, cluster_power_OTOI);
    
    zenith_dep_LOS_all = cat(3, zenith_dep_LOS_all, zenith_dep_LOS);
    zenith_arr_LOS_all = cat(3, zenith_arr_LOS_all, zenith_arr_LOS);
    cluster_power_LOS_all = cat(3, cluster_power_LOS_all, cluster_power_LOS);
end

% Load the results extracted from 3GPP R1-143469
% ZSD and ZSA are saved for UMa and UMi scenarios
load('./Extracted_calibration_values_3GPP_R1_143469.mat')

%% Circular Angle Spread TR 25.996
% Apply delta shift
% delta = [10:10:100];

% Add delta shift to all angles
ones_matrix = ones(20,cluster_NLOS,size(zenith_dep_NLOS_all,3));
delta = 10:10:200;
for f = 1:length(delta)
    for k=1:size(zenith_dep_NLOS_all,3)
        for j = 1:size(zenith_dep_NLOS_all,2)
            for i = 1:20
                unit_matrix(i,j,k,f)=ones_matrix(i,j,k)*delta(f);
            end
        end
    end
end

for i_=1:size(unit_matrix,4)
%     azimuth_dep_NLOS_i(:,:,:,i_)=azimuth_dep_NLOS + unit_matrix(:,:,:,i_);
%     azimuth_arr_NLOS_i(:,:,:,i_)=azimuth_arr_NLOS + unit_matrix(:,:,:,i_);
    zenith_dep_NLOS_i(:,:,:,i_) = zenith_dep_NLOS_all + unit_matrix(:,:,:,i_);
    zenith_arr_NLOS_i(:,:,:,i_) = zenith_arr_NLOS_all + unit_matrix(:,:,:,i_);
end


% Define myy variable
for d = 1:size(zenith_dep_NLOS_i,4)
    for i = 1:size(zenith_dep_NLOS_i,3)
%         myy_azimuth_dep_NLOS(i,d) = sum(sum(azimuth_dep_NLOS_i(:,:,i,d).*cluster_power_NLOS_all(:,:,i)))./sum(sum(cluster_power_NLOS_all(:,:,i)));
%         myy_azimuth_arr_NLOS(i,d) = sum(sum(azimuth_arr_NLOS_i(:,:,i,d).*cluster_power_NLOS_all(:,:,i)))./sum(sum(cluster_power_NLOS_all(:,:,i)));
        myy_zenith_dep_NLOS(i,d) = sum(sum(zenith_dep_NLOS_i(:,:,i,d).*cluster_power_NLOS_all(:,:,i)))./sum(sum(cluster_power_NLOS_all(:,:,i)));
        myy_zenith_arr_NLOS(i,d) = sum(sum(zenith_arr_NLOS_i(:,:,i,d).*cluster_power_NLOS_all(:,:,i)))./sum(sum(cluster_power_NLOS_all(:,:,i)));
    end
end


for d = 1:size(zenith_dep_NLOS_i,4)
    for k = 1:size(zenith_dep_NLOS_i,3)
        for j= 1:size(zenith_dep_NLOS_i,2)
            for i = 1:size(zenith_dep_NLOS_i,1)
                if (zenith_dep_NLOS_i(i,j,k,d)-myy_zenith_dep_NLOS(k,d)) < -180
                    THETA_zenith_dep_NLOS(i,j,k,d) = 360 + (zenith_dep_NLOS_i(i,j,k,d)-myy_zenith_dep_NLOS(k,d));
                elseif abs(zenith_dep_NLOS_i(i,j,k,d)-myy_zenith_dep_NLOS(k,d)) <= 180 
                    THETA_zenith_dep_NLOS(i,j,k,d) = zenith_dep_NLOS_i(i,j,k,d)-myy_zenith_dep_NLOS(k,d);
                else (zenith_dep_NLOS_i(i,j,k,d)-myy_zenith_dep_NLOS(k,d)) > 180
                    THETA_zenith_dep_NLOS(i,j,k,d) = 360 - (zenith_dep_NLOS_i(i,j,k,d)-myy_zenith_dep_NLOS(k,d));
                end
                
                if (zenith_arr_NLOS_i(i,j,k,d)-myy_zenith_arr_NLOS(k,d)) < -180
                    THETA_zenith_arr_NLOS(i,j,k,d) = 360 + (zenith_arr_NLOS_i(i,j,k,d)-myy_zenith_arr_NLOS(k,d));
                elseif abs(zenith_arr_NLOS_i(i,j,k,d)-myy_zenith_arr_NLOS(k,d)) <= 180 
                    THETA_zenith_arr_NLOS(i,j,k,d) = zenith_arr_NLOS_i(i,j,k,d)-myy_zenith_arr_NLOS(k,d);
                else (zenith_arr_NLOS_i(i,j,k,d)-myy_zenith_arr_NLOS(k,d)) > 180
                    THETA_zenith_arr_NLOS(i,j,k,d) = 360 - (zenith_arr_NLOS_i(i,j,k,d)-myy_zenith_arr_NLOS(k,d));
                end
                
%                 sigma_zenith_dep_NLOS(k,d)= sqrt((sum(sum(THETA_zenith_dep_NLOS(:,:,k,d).^2.*cluster_power_NLOS_all(:,:,k))))./(sum(sum(cluster_power_NLOS_all(:,:,k)))));
%                 sigma_zenith_arr_NLOS(k,d)= sqrt((sum(sum(THETA_zenith_arr_NLOS(:,:,k,d).^2.*cluster_power_NLOS_all(:,:,k))))./(sum(sum(cluster_power_NLOS_all(:,:,k)))));
            end
        end
         sigma_zenith_dep_NLOS(k,d)= sqrt((sum(sum(THETA_zenith_dep_NLOS(:,:,k,d).^2.*cluster_power_NLOS_all(:,:,k))))./(sum(sum(cluster_power_NLOS_all(:,:,k)))));
         sigma_zenith_arr_NLOS(k,d)= sqrt((sum(sum(THETA_zenith_arr_NLOS(:,:,k,d).^2.*cluster_power_NLOS_all(:,:,k))))./(sum(sum(cluster_power_NLOS_all(:,:,k)))));
    end
end



% sigma_azimuth_dep_NLOS= min(sigma_azimuth_dep_NLOS,[],2);
% sigma_azimuth_arr_NLOS = min(sigma_azimuth_arr_NLOS,[],2);
sigma_zenith_dep_NLOS = min(sigma_zenith_dep_NLOS,[],2);
sigma_zenith_arr_NLOS = min(sigma_zenith_arr_NLOS,[],2);



% Add delta shift to all angles
ones_matrix_OTOI = ones(20,cluster_OTOI,size(zenith_dep_OTOI_all,3));
% delta = 10:10:200;
for f = 1:length(delta)
    for k=1:size(zenith_dep_OTOI_all,3)
        for j = 1:size(zenith_dep_OTOI_all,2)
            for i = 1:20
                unit_matrix_OTOI(i,j,k,f)=ones_matrix_OTOI(i,j,k)*delta(f);
            end
        end
    end
end

for i_=1:size(unit_matrix_OTOI,4)
%     azimuth_dep_OTOI_i(:,:,:,i_)=azimuth_dep_OTOI + unit_matrix_OTOI(:,:,:,i_);
%     azimuth_arr_OTOI_i(:,:,:,i_)=azimuth_arr_OTOI + unit_matrix_OTOI(:,:,:,i_);
    zenith_dep_OTOI_i(:,:,:,i_) = zenith_dep_OTOI_all + unit_matrix_OTOI(:,:,:,i_);
    zenith_arr_OTOI_i(:,:,:,i_) = zenith_arr_OTOI_all + unit_matrix_OTOI(:,:,:,i_);
end


% Define myy variable
for d = 1:size(zenith_dep_OTOI_i,4)
    for i = 1:size(zenith_dep_OTOI_i,3)
%         myy_azimuth_dep_OTOI(i,d) = sum(sum(azimuth_dep_OTOI_i(:,:,i,d).*cluster_power_OTOI_all(:,:,i)))./sum(sum(cluster_power_OTOI_all(:,:,i)));
%         myy_azimuth_arr_OTOI(i,d) = sum(sum(azimuth_arr_OTOI_i(:,:,i,d).*cluster_power_OTOI_all(:,:,i)))./sum(sum(cluster_power_OTOI_all(:,:,i)));
        myy_zenith_dep_OTOI(i,d) = sum(sum(zenith_dep_OTOI_i(:,:,i,d).*cluster_power_OTOI_all(:,:,i)))./sum(sum(cluster_power_OTOI_all(:,:,i)));
        myy_zenith_arr_OTOI(i,d) = sum(sum(zenith_arr_OTOI_i(:,:,i,d).*cluster_power_OTOI_all(:,:,i)))./sum(sum(cluster_power_OTOI_all(:,:,i)));
    end
end



for d = 1:size(zenith_dep_OTOI_i,4)
    for k = 1:size(zenith_dep_OTOI_i,3)
        for j= 1:size(zenith_dep_OTOI_i,2)
            for i = 1:size(zenith_dep_OTOI_i,1)
                if (zenith_dep_OTOI_i(i,j,k,d)-myy_zenith_dep_OTOI(k,d)) < -180
                    THETA_zenith_dep_OTOI(i,j,k,d) = 360 + (zenith_dep_OTOI_i(i,j,k,d)-myy_zenith_dep_OTOI(k,d));
                elseif abs(zenith_dep_OTOI_i(i,j,k,d)-myy_zenith_dep_OTOI(k,d)) <= 180
                    THETA_zenith_dep_OTOI(i,j,k,d) = zenith_dep_OTOI_i(i,j,k,d)-myy_zenith_dep_OTOI(k,d);
                else (zenith_dep_OTOI_i(i,j,k,d)-myy_zenith_dep_OTOI(k,d)) > 180
                    THETA_zenith_dep_OTOI(i,j,k,d) = 360 - (zenith_dep_OTOI_i(i,j,k,d)-myy_zenith_dep_OTOI(k,d));
                end
                
                if (zenith_arr_OTOI_i(i,j,k,d)-myy_zenith_arr_OTOI(k,d)) < -180
                    THETA_zenith_arr_OTOI(i,j,k,d) = 360 + (zenith_arr_OTOI_i(i,j,k,d)-myy_zenith_arr_OTOI(k,d));
                elseif abs(zenith_arr_OTOI_i(i,j,k,d)-myy_zenith_arr_OTOI(k,d)) <= 180
                    THETA_zenith_arr_OTOI(i,j,k,d) = zenith_arr_OTOI_i(i,j,k,d)-myy_zenith_arr_OTOI(k,d);
                else (zenith_arr_OTOI_i(i,j,k,d)-myy_zenith_arr_OTOI(k,d)) > 180
                    THETA_zenith_arr_OTOI(i,j,k,d) = 360 - (zenith_arr_OTOI_i(i,j,k,d)-myy_zenith_arr_OTOI(k,d));
                end
                
            end
        end
        sigma_zenith_dep_OTOI(k,d)= sqrt((sum(sum(THETA_zenith_dep_OTOI(:,:,k,d).^2.*cluster_power_OTOI_all(:,:,k))))./(sum(sum(cluster_power_OTOI_all(:,:,k)))));
        sigma_zenith_arr_OTOI(k,d)= sqrt((sum(sum(THETA_zenith_arr_OTOI(:,:,k,d).^2.*cluster_power_OTOI_all(:,:,k))))./(sum(sum(cluster_power_OTOI_all(:,:,k)))));
    end
end


% sigma_azimuth_dep_OTOI= min(sigma_azimuth_dep_OTOI,[],2);
% sigma_azimuth_arr_OTOI = min(sigma_azimuth_arr_OTOI,[],2);
sigma_zenith_dep_OTOI = min(sigma_zenith_dep_OTOI,[],2);
sigma_zenith_arr_OTOI = min(sigma_zenith_arr_OTOI,[],2);


% Add delta shift to all angles
ones_matrix_LOS = ones(20,cluster_LOS,size(zenith_dep_LOS_all,3));
for f = 1:length(delta)
    for k=1:size(zenith_dep_LOS_all,3)
        for j = 1:cluster_LOS
            for i = 1:20
                unit_matrix_LOS(i,j,k,f)=ones_matrix_LOS(i,j,k)*delta(f);
%                   unit_matrix(i,j,k,f)=unit_matrix(i,j,k)-f*10;
            end
        end
    end
end

for i_=1:size(unit_matrix_LOS,4)
% azimuth_dep_LOS_i(:,:,:,i_)=azimuth_dep_LOS +unit_matrix_LOS(:,:,:,i_);
% azimuth_arr_LOS_i(:,:,:,i_)=azimuth_arr_LOS +unit_matrix_LOS(:,:,:,i_);
zenith_dep_LOS_i(:,:,:,i_) = zenith_dep_LOS_all + unit_matrix_LOS(:,:,:,i_);
zenith_arr_LOS_i(:,:,:,i_) = zenith_arr_LOS_all + unit_matrix_LOS(:,:,:,i_);

end


% Define myy variable
for d = 1:size(zenith_dep_LOS_i,4)
    for i = 1:size(zenith_dep_LOS_i,3)
%         myy_azimuth_dep_LOS(i,d) = sum(sum(azimuth_dep_LOS_i(:,:,i,d).*cluster_power_LOS_all(:,:,i)))./sum(sum(cluster_power_LOS_all(:,:,i)));
%         myy_azimuth_arr_LOS(i,d) = sum(sum(azimuth_arr_LOS_i(:,:,i,d).*cluster_power_LOS_all(:,:,i)))./sum(sum(cluster_power_LOS_all(:,:,i)));
        myy_zenith_dep_LOS(i,d) = sum(sum(zenith_dep_LOS_i(:,:,i,d).*cluster_power_LOS_all(:,:,i)))./sum(sum(cluster_power_LOS_all(:,:,i)));
        myy_zenith_arr_LOS(i,d) = sum(sum(zenith_arr_LOS_i(:,:,i,d).*cluster_power_LOS_all(:,:,i)))./sum(sum(cluster_power_LOS_all(:,:,i)));
    end
end



for d = 1:size(zenith_dep_LOS_i,4)
    for k = 1:size(zenith_dep_LOS_i,3)
        for j= 1:size(zenith_dep_LOS_i,2)
            for i = 1:size(zenith_dep_LOS_i,1)
                if (zenith_dep_LOS_i(i,j,k,d)-myy_zenith_dep_LOS(k,d)) < -180
                    THETA_zenith_dep_LOS(i,j,k,d) = 360 + (zenith_dep_LOS_i(i,j,k,d)-myy_zenith_dep_LOS(k,d));
                elseif abs(zenith_dep_LOS_i(i,j,k,d)-myy_zenith_dep_LOS(k,d)) <= 180
                    THETA_zenith_dep_LOS(i,j,k,d) = zenith_dep_LOS_i(i,j,k,d)-myy_zenith_dep_LOS(k,d);
                else (zenith_dep_LOS_i(i,j,k,d)-myy_zenith_dep_LOS(k,d)) > 180
                    THETA_zenith_dep_LOS(i,j,k,d) = 360 - (zenith_dep_LOS_i(i,j,k,d)-myy_zenith_dep_LOS(k,d));
                end
                
                if (zenith_arr_LOS_i(i,j,k,d)-myy_zenith_arr_LOS(k,d)) < -180
                    THETA_zenith_arr_LOS(i,j,k,d) = 360 + (zenith_arr_LOS_i(i,j,k,d)-myy_zenith_arr_LOS(k,d));
                elseif abs(zenith_arr_LOS_i(i,j,k,d)-myy_zenith_arr_LOS(k,d)) <= 180
                    THETA_zenith_arr_LOS(i,j,k,d) = zenith_arr_LOS_i(i,j,k,d)-myy_zenith_arr_LOS(k,d);
                else (zenith_arr_LOS_i(i,j,k,d)-myy_zenith_arr_LOS(k,d)) > 180
                    THETA_zenith_arr_LOS(i,j,k,d) = 360 - (zenith_arr_LOS_i(i,j,k,d)-myy_zenith_arr_LOS(k,d));
                end
                
                
            end
        end
        sigma_zenith_dep_LOS(k,d)= sqrt((sum(sum(THETA_zenith_dep_LOS(:,:,k,d).^2.*cluster_power_LOS_all(:,:,k))))./(sum(sum(cluster_power_LOS_all(:,:,k)))));
        sigma_zenith_arr_LOS(k,d)= sqrt((sum(sum(THETA_zenith_arr_LOS(:,:,k,d).^2.*cluster_power_LOS_all(:,:,k))))./(sum(sum(cluster_power_LOS_all(:,:,k)))));
    end
end


% sigma_azimuth_dep_LOS = min(sigma_azimuth_dep_LOS,[],2);
% sigma_azimuth_arr_LOS = min(sigma_azimuth_arr_LOS,[],2);
sigma_zenith_dep_LOS = min(sigma_zenith_dep_LOS,[],2);
sigma_zenith_arr_LOS = min(sigma_zenith_arr_LOS,[],2);

% circ_azimuth_departure   = [sigma_azimuth_dep_NLOS',sigma_azimuth_dep_LOS', sigma_azimuth_dep_OTOI'];
% circ_azimuth_arrival     = [sigma_azimuth_arr_NLOS',sigma_azimuth_arr_LOS', sigma_azimuth_arr_OTOI'];
circ_zenith_departure    = [sigma_zenith_dep_NLOS',sigma_zenith_dep_LOS', sigma_zenith_dep_OTOI'];
circ_zenith_arrival      = [sigma_zenith_arr_NLOS',sigma_zenith_arr_LOS', sigma_zenith_arr_OTOI'];




%% PLOT Large Scale Parameters: UMa

figure(1); hold on; grid on;
%Plot 3GPP calibration results from
%'Extracted_calibration_values_3GPP_R1_143469.mat'
%%%%%%% ZSD
for i = 1:21
[f,x]=ecdf(ZSD_UMa(:,i)); 
 plot(x, f, 'Color', [0.5 0.5 0.5],'LineWidth',2); xlabel('Zenith spread of departure [°]'); ylabel('ECDF');xlim([0 60]);
end

figure(2); hold on; grid on; 
[f22,x22]=ecdf(circ_zenith_departure); [f23,x23]=ecdf(ZSD_UMa(:,22));
plot(x22,f22,'r',x23,f23,'g');xlabel('Zenith spread of departure [°]'); ylabel('ECDF'); xlim([0 60]);

copyobj(findobj(figure(2), 'Type', 'line'), gca(figure(1)));

%%%%%%% ZSA
figure(1); hold on; grid on;
for i = 1:21
[f,x]=ecdf(ZSA_UMa(:,i)); 
 plot(x, f, 'Color', [0.5 0.5 0.5],'LineWidth',2); xlabel('Zenith spread of arrival [°]'); ylabel('ECDF');xlim([0 60]);
end

figure(2); hold on; grid on; 
[f22,x22]=ecdf(circ_zenith_arrival); [f23,x23]=ecdf(ZSA_UMa(:,22)); 
plot(x22,f22,'r',x23,f23,'g');xlabel('Zenith spread of arrival [°]'); ylabel('ECDF'); xlim([0 60]);

copyobj(findobj(figure(2), 'Type', 'line'), gca(figure(1)));

%% PLOT Large Scale Parameters: UMi
%Plot 3GPP calibration results from
%'Extracted_calibration_values_3GPP_R1_143469.mat'
%%%%%%% ZSD
figure(1); hold on; grid on;
for i = 1:21
[f,x]=ecdf(ZSD_UMi(:,i)); 
 plot(x, f, 'Color', [0.5 0.5 0.5],'LineWidth',2); xlabel('Zenith spread of departure [°]'); ylabel('ECDF');xlim([0 60]);
end

figure(2); hold on; grid on; 
[f22,x22]=ecdf(circ_zenith_departure); [f23,x23]=ecdf(ZSD_UMi(:,22));
plot(x22,f22,'b',x23,f23,'g');xlabel('Zenith spread of depatrure [°]'); ylabel('ECDF'); xlim([0 60]);

copyobj(findobj(figure(2), 'Type', 'line'), gca(figure(1)));

%%%%%%% ZSA
figure(1); hold on; grid on;
for i = 1:21
[f,x]=ecdf(ZSA_UMi(:,i)); 
 plot(x, f, 'Color', [0.5 0.5 0.5],'LineWidth',2); xlabel('Zenith spread of arrival [°]'); ylabel('ECDF');xlim([0 60]);
end

figure(2); hold on; grid on; 
[f22,x22]=ecdf(circ_zenith_arrival); [f23,x23]=ecdf(ZSA_UMi(:,22)); 
plot(x22,f22,'b',x23,f23,'g');xlabel('Zenith spread of arrival [°]'); ylabel('ECDF'); xlim([0 60]);

copyobj(findobj(figure(2), 'Type', 'line'), gca(figure(1)));






