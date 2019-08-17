function config = load_UMa_uplink_parameters(config)
% 3D-UMa input parameters from TR 36.873
% config.C_constant                 - Scaling factors for AOA, AOD generation Table 7.3-2
% config.C_constant_elevation       - Scaling factors for ZOA, ZOD generation Table 7.3.4
% config.PerClusterRays             - Number of rays per cluster
% config.NumClusters_LOS            - Number of clusters
% config.r_DS_LOS                   - Delays spread proportionality factor
% config.PerClusterAS_D_LOS         - Per cluster rms azimuth spread of departure angles [deg]
% config.PerClusterAS_A_LOS         - Per cluster rms azimuth spread of arrival angles [deg]
% config.PerClusterZS_A_LOS         - Per cluster rms zenith spread of arrival angles
% config.LNS_ksi_LOS                - Per cluster shadowing std
% config.asD_ds_LOS                 - Cross correlation coefficient for departure AS vs delay spread
% config.asA_ds_LOS                 - Cross correlation coefficient for arrival AS vs delay spread
%
%
% (c) Fjolla Ademaj, Martin Taranetz, ITC 2016

config.C_constant                   = [4,5,8,10,11,12,14,15,16,19,20;...
                        0.779,0.860,1.018,1.090,1.123,1.146,1.190,1.211,1.226,1.273,1.289];
config.C_constant_elevation         = [12,19,20;1.104,1.184,1.178];

config.PerClusterRays               = 20;  % Same for all scenarios and LOS/NLOS cases
%All input parameters as in table 7.3-6
% 3DUMa, LOS   (3D Urban Macro)
% Fixed scenario specific parameters
config.NumClusters_LOS = 12;         % Number of clusters    [Table 7.3-6]
config.r_DS_LOS   = 2.5;             % Delays spread proportionality factor
config.PerClusterAS_D_LOS = 11;      % Per cluster rms azimuth spread of departure angles [deg]
config.PerClusterAS_A_LOS = 5;       % Per cluster rms azimuth spread of arrival angles [deg]
% NOTE: PerClusterZS_A_LOS not used, so no swap needed
config.PerClusterZS_A_LOS = 7;       % Per cluster rms zenith spread of arrival angles
config.LNS_ksi_LOS = 3;              % LNS ksi [dB], per cluster shadowing std
% Cross correlation coefficients for azimuth case
config.asD_ds_LOS = 0.8;             % departure AS vs delay spread
config.asA_ds_LOS = 0.4;             % arrival AS vs delay spread
config.asA_sf_LOS = -0.5;            % arrival AS vs shadowing std
config.asD_sf_LOS = -0.5;            % departure AS vs shadowing std
config.ds_sf_LOS  = -0.4;            % delay spread vs shadowing std
config.asD_asA_LOS = 0;              % departure AS vs arrival AS
config.asD_kf_LOS = -0.2;            % departure AS vs k-factor
config.asA_kf_LOS = 0;               % arrival AS vs k-factor
config.ds_kf_LOS = -0.4;             % delay spread vs k-factor
config.sf_kf_LOS = 0;                % shadowing std vs k-factor
% Cross correlation coefficients for elevation (zenith) case
config.zsD_sf_LOS = -0.8;            % departure ZS vs shadowing std
config.zsA_sf_LOS = 0;               % arrival ZS vs shadowing std
config.zsD_kf_LOS = 0;               % departure ZS vs k-factor
config.zsA_kf_LOS = 0;               % arrival ZS vs k-factor
config.zsD_ds_LOS = 0;               % departure ZS vs delay spread
config.zsA_ds_LOS = -0.2;            % arrival ZS vs delay spread
config.zsD_asD_LOS = 0.4;            % departure ZS vs departure AS
config.zsA_asD_LOS = -0.3;           % arrival ZS vs departure AS
config.zsD_asA_LOS = 0;              % departure ZS vs arrival AS
config.zsA_asA_LOS = 0.5;            % arrival ZS vs arrival AS
config.zsD_zsA_LOS = 0;              % departure ZS vs arrival ZS
% Polarisation parameters
config.xpr_mu_LOS    = 8;            % XPR mean [dB]
config.xpr_sigma_LOS = 4;            % XPR std  [dB]
% Dispersion parameters [1, Table 4.5]. Log-normal distributions
config.DS_mu_LOS      = -7.03;       % delay spread, mean [log10(s)]
config.DS_sigma_LOS   = 0.66;        % delay spread, std [log10(s)]
config.AS_D_mu_LOS    = 1.81;        % azimuth departure angle spread, mean [log10(deg)]
config.AS_D_sigma_LOS = 0.2;         % azimuth departure angle spread, std [log10(deg)]
config.AS_A_mu_LOS    = 1.15;        % azimuth arrival angle spread, mean [log10(deg)]
config.AS_A_sigma_LOS = 0.28;        % azimuth arrival angle spread, std [log10(deg)]
% NOTE: ZS_A_mu_LOS and ZS_A_sigma_LOS not swapped here but in
% network_elements.eNodeB.sigmas_LOS()
config.ZS_A_mu_LOS    = 0.95;        % zenith arrival angle spread, mean [log10(deg)]
config.ZS_A_sigma_LOS = 0.16;        % zenith arrival angle spread, std [log10(deg)]
config.SF_sigma_LOS   = 4;           % shadowing std [dB] (zero mean)
config.KF_mu_LOS      = 9;           % K-factor mean [dB]
config.KF_sigma_LOS   = 3.5;         % K-factor std [dB]
% "Correlation distances" in horizontal plane [m]
config.DS_lambda_LOS   = 30;         % [m], delay spread
config.AS_D_lambda_LOS = 15;         % [m], departure azimuth spread
config.AS_A_lambda_LOS = 18;         % [m], arrival azimuth spread
config.SF_lambda_LOS   = 37;         % [m], shadowing
config.KF_lambda_LOS   = 12;         % [m], k-factor
config.ZS_A_lambda_LOS = 15;         % [m], arrival zenith spread
config.ZS_D_lambda_LOS = 15;         % [m], departure zenith spread
% Cluster angle spread value
config.cluster_ASD_LOS = 11;
config.cluster_ASA_LOS = 5;
% NOTE: cluster_ZSA_LOS not swapped here but in
% channel_gain_wrappers.TR36873_Fading_3D_Channel.zenith_angle_of_arrival_LOS()
% and zenith_angle_of_departure_LOS()
config.cluster_ZSA_LOS = 7;

% 3DUMa, NLOS  (3D Urban Macro)
% Fixed scenario specific parameters
config.NumClusters_NLOS = 20;         % Number of clusters    [Table 7.3-6]
config.r_DS_NLOS   = 2.3;             % Delays spread proportionality factor
config.PerClusterAS_D_NLOS = 15;      % Per cluster rms azimuth spread of departure angles [deg]
config.PerClusterAS_A_NLOS = 2;       % Per cluster rms azimuth spread of arrival angles [deg]
% NOTE: PerClusterZS_A_NLOS not used, so no swap needed
config.PerClusterZS_A_NLOS = 7;       % Per cluster rms zenith spread of arrival angles
config.LNS_ksi_NLOS = 3;              % LNS ksi [dB], per cluster shadowing std
% Cross correlation coefficients for azimuth case
config.asD_ds_NLOS = 0.6;             % departure AS vs delay spread
config.asA_ds_NLOS = 0.4;             % arrival AS vs delay spread
config.asA_sf_NLOS = -0.6;            % arrival AS vs shadowing std
config.asD_sf_NLOS = 0;               % departure AS vs shadowing std
config.ds_sf_NLOS  = -0.4;            % delay spread vs shadowing std
config.asD_asA_NLOS = 0.4;            % departure AS vs arrival AS
% Cross correlation coefficients for elevation (zenith) case
config.zsD_sf_NLOS = -0.4;            % departure ZS vs shadowing std
config.zsA_sf_NLOS = 0;               % arrival ZS vs shadowing std
config.zsD_ds_NLOS = 0;               % departure ZS vs delay spread
config.zsA_ds_NLOS = -0.5;            % arrival ZS vs delay spread
config.zsD_asD_NLOS = 0;              % departure ZS vs departure AS
config.zsA_asD_NLOS = 0;              % arrival ZS vs departure AS
config.zsD_asA_NLOS = -0.1;           % departure ZS vs arrival AS
config.zsA_asA_NLOS = 0.5;            % arrival ZS vs arrival AS
config.zsD_zsA_NLOS = 0;              % departure ZS vs arrival ZS
% Polarisation parameters
config.xpr_mu_NLOS    = 7;            % XPR mean [dB]
config.xpr_sigma_NLOS = 3;            % XPR std  [dB]
% Dispersion parameters [1, Table 4.5].Log-normal distributions
config.DS_mu_NLOS      = -6.44;       % delay spread, mean [log10(s)]
config.DS_sigma_NLOS   = 0.39;        % delay spread, std [log10(s)]
config.AS_D_mu_NLOS    = 1.87;        % azimuth departure angle spread, mean [log10(deg)]
config.AS_D_sigma_NLOS = 0.11;        % azimuth departure angle spread, std [log10(deg)]
config.AS_A_mu_NLOS    = 1.41;        % azimuth arrival angle spread, mean [log10(deg)]
config.AS_A_sigma_NLOS = 0.28;        % azimuth arrival angle spread, std [log10(deg)]
% NOTE: ZS_A_mu_NLOS and ZS_A_sigma_NLOS not swapped here but in
% network_elements.eNodeB.sigmas_NLOS()
config.ZS_A_mu_NLOS    = 1.26;        % zenith arrival angle spread, mean [log10(deg)]
config.ZS_A_sigma_NLOS = 0.16;        % zenith arrival angle spread, std [log10(deg)]
config.SF_sigma_NLOS   = 6;           % shadowing std [dB] (zero mean)
% "Decorrelation distances" in horizontal plane [m]
config.DS_lambda_NLOS   = 40;         % [m], delay spread
config.AS_D_lambda_NLOS = 50;         % [m], departure azimuth spread
config.AS_A_lambda_NLOS = 50;         % [m], arrival azimuth spread
config.SF_lambda_NLOS   = 50;         % [m], shadowing
config.ZS_A_lambda_NLOS = 50;         % [m], arrival zenith spread
config.ZS_D_lambda_NLOS = 50;         % [m], departure zenith spread
% CLuster angle spread value
config.cluster_ASD_NLOS = 15;
config.cluster_ASA_NLOS = 2;
% NOTE: cluster_ZSA_NLOS not swapped here but in
% channel_gain_wrappers.TR36873_Fading_3D_Channel.zenith_angle_of_arrival_NLOS()
% and zenith_angle_of_departure_NLOS()
config.cluster_ZSA_NLOS = 7;

%%Added from tri_sector_3D line 271 to 319
%Change 'LTE_config' for config

 % 3DUMa, OTOI  (3D Urban Macro) Outdoor-TO-Indoor
% Fixed scenario specific parameters
config.NumClusters_OTOI = 12;         % Number of clusters    [Table 7.3-6]
config.r_DS_OTOI   = 2.2;             % Delays spread proportionality factor
config.PerClusterAS_D_OTOI = 8;       % Per cluster rms azimuth spread of departure angles [deg]
config.PerClusterAS_A_OTOI = 5;       % Per cluster rms azimuth spread of arrival angles [deg]
% NOTE: PerClusterZS_A_OTOI not used, so no swap needed
config.PerClusterZS_A_OTOI = 3;       % Per cluster rms zenith spread of arrival angles
config.LNS_ksi_OTOI = 4;              % LNS ksi [dB], per cluster shadowing std
% Cross correlation coefficients for azimuth case
config.asD_ds_OTOI = 0.4;             % departure AS vs delay spread
config.asA_ds_OTOI = 0.4;             % arrival AS vs delay spread
config.asA_sf_OTOI = 0.2;             % arrival AS vs shadowing std
config.asD_sf_OTOI = 0;               % departure AS vs shadowing std
config.ds_sf_OTOI  = -0.5;            % delay spread vs shadowing std
config.asD_asA_OTOI = 0;              % departure AS vs arrival AS
% Cross correlation coefficients for elevation (zenith) case
config.zsD_sf_OTOI = 0;               % departure ZS vs shadowing std
config.zsA_sf_OTOI = 0;               % arrival ZS vs shadowing std
config.zsD_ds_OTOI = -0.2;            % departure ZS vs delay spread
config.zsA_ds_OTOI = -0.6;            % arrival ZS vs delay spread
config.zsD_asD_OTOI = 0.5;            % departure ZS vs departure AS
config.zsA_asD_OTOI = 0;              % arrival ZS vs departure AS
config.zsD_asA_OTOI = 0;              % departure ZS vs arrival AS
config.zsA_asA_OTOI = -0.2;           % arrival ZS vs arrival AS
config.zsD_zsA_OTOI = 0.5;            % departure ZS vs arrival ZS
% Polarisation parameters
config.xpr_mu_OTOI    = 9;            % XPR mean [dB]
config.xpr_sigma_OTOI = 5;            % XPR std  [dB]
% Dispersion parameters [1, Table 4.5].Log-normal distributions
config.DS_mu_OTOI      = -6.62;       % delay spread, mean [log10(s)]
config.DS_sigma_OTOI   = 0.32;        % delay spread, std [log10(s)]
config.AS_D_mu_OTOI    = 1.76;        % azimuth departure angle spread, mean [log10(deg)]
config.AS_D_sigma_OTOI = 0.16;        % azimuth departure angle spread, std [log10(deg)]
config.AS_A_mu_OTOI    = 1.25;        % azimuth arrival angle spread, mean [log10(deg)]
config.AS_A_sigma_OTOI = 0.42;        % azimuth arrival angle spread, std [log10(deg)]
% NOTE: ZS_A_mu_OTOI and ZS_A_sigma_OTOI not swapped here but in
% network_elements.eNodeB.sigmas_OTOI()
config.ZS_A_mu_OTOI    = 1.01;        % zenith arrival angle spread, mean [log10(deg)]
config.ZS_A_sigma_OTOI = 0.43;        % zenith arrival angle spread, std [log10(deg)]
config.SF_sigma_OTOI   = 7;           % shadowing std [dB] (zero mean)
% "Decorrelation distances" in horizontal plane [m]
config.DS_lambda_OTOI   = 10;         % [m], delay spread
config.AS_D_lambda_OTOI = 17;         % [m], departure azimuth spread
config.AS_A_lambda_OTOI = 11;         % [m], arrival azimuth spread
config.SF_lambda_OTOI   = 7;          % [m], shadowing
config.ZS_A_lambda_OTOI = 25;         % [m], arrival zenith spread
config.ZS_D_lambda_OTOI = 25;         % [m], departure zenith spread
% CLuster angle spread value
config.cluster_ASD_OTOI = 8;
config.cluster_ASA_OTOI = 5;
% NOTE: cluster_ZSA_OTOI not swapped here but in
% channel_gain_wrappers.TR36873_Fading_3D_Channel.zenith_angle_of_arrival_OTOI()
% and zenith_angle_of_departure_OTOI()
config.cluster_ZSA_OTOI = 3; 
end