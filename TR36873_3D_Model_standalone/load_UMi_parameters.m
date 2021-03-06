function config = load_UMi_parameters(config)
% 3D-UMi input paramaters from TR 36.873
% (c) Fjolla Ademaj, Martin Taranetz, ITC 2016

config.C_constant                   = [4,5,8,10,11,12,14,15,16,19,20;...
                        0.779,0.860,1.018,1.090,1.123,1.146,1.190,1.211,1.226,1.273,1.289];
config.C_constant_elevation         = [12,19,20;1.104,1.184,1.178];

config.PerClusterRays               = 20;  % Same for all scenarios and LOS/NLOS cases
% 3DUMi, LOS   (3D Urban Micro)
% Fixed scenario specific parameters
config.NumClusters_LOS = 12;         % Number of clusters    [Table 7.3-6]
config.r_DS_LOS   = 3.2;             % Delays spread proportionality factor
config.PerClusterAS_D_LOS = 3;       % Per cluster rms azimuth spread of departure angles [deg]
config.PerClusterAS_A_LOS = 17;      % Per cluster rms azimuth spread of arrival angles [deg]
config.PerClusterZS_A_LOS = 7;       % Per cluster rms zenith spread of arrival angles
config.LNS_ksi_LOS = 3;                 % LNS ksi [dB], per cluster shadowing std
% Cross correlation coefficients for azimuth case
config.asD_ds_LOS = 0.5;             % departure AS vs delay spread
config.asA_ds_LOS = 0.8;             % arrival AS vs delay spread
config.asA_sf_LOS = -0.4;            % arrival AS vs shadowing std
config.asD_sf_LOS = -0.5;            % departure AS vs shadowing std
config.ds_sf_LOS  = -0.4;            % delay spread vs shadowing std
config.asD_asA_LOS = 0.4;            % departure AS vs arrival AS
config.asD_kf_LOS = -0.2;            % departure AS vs k-factor
config.asA_kf_LOS = -0.3;            % arrival AS vs k-factor
config.ds_kf_LOS = -0.7;             % delay spread vs k-factor
config.sf_kf_LOS = 0.5;              % shadowing std vs k-factor
% Cross correlation coefficients for elevation (zenith) case
config.zsD_sf_LOS = 0;               % departure ZS vs shadowing std
config.zsA_sf_LOS = 0;               % arrival ZS vs shadowing std
config.zsD_kf_LOS = 0;               % departure ZS vs k-factor
config.zsA_kf_LOS = 0;               % arrival ZS vs k-factor
config.zsD_ds_LOS  = 0;              % departure ZS vs delay spread
config.zsA_ds_LOS = 0.2;             % arrival ZS vs delay spread
config.zsD_asD_LOS = 0.5;            % departure ZS vs departure AS
config.zsA_asD_LOS = 0.3;            % arrival ZS vs departure AS
config.zsD_asA_LOS = 0;              % departure ZS vs arrival AS
config.zsA_asA_LOS = 0;              % arrival ZS vs arrival AS
config.zsD_zsA_LOS = 0;              % departure ZS vs arrival ZS
% Polarisation parameters
config.xpr_mu_LOS    = 9;            % XPR mean [dB]
config.xpr_sigma_LOS = 3;            % XPR std  [dB]
% Dispersion parameters [1, Table 4.5]. Log-normal distributions
config.DS_mu_LOS      = -7.19;       % delay spread, mean [log10(s)]
config.DS_sigma_LOS   = 0.40;        % delay spread, std [log10(s)]
config.AS_D_mu_LOS    = 1.20;        % azimuth departure angle spread, mean [log10(deg)]
config.AS_D_sigma_LOS = 0.43;        % azimuth departure angle spread, std [log10(deg)]
config.AS_A_mu_LOS    = 1.75;        % azimuth arrival angle spread, mean [log10(deg)]
config.AS_A_sigma_LOS = 0.19;        % azimuth arrival angle spread, std [log10(deg)]
config.ZS_A_mu_LOS    = 0.60;        % zenith arrival angle spread, mean [log10(deg)]
config.ZS_A_sigma_LOS = 0.16;        % zenith arrival angle spread, std [log10(deg)]
config.SF_sigma_LOS   = 3;           % shadowing std [dB] (zero mean)
config.KF_mu_LOS = 9;                % K-factor mean [dB]
config.KF_sigma_LOS = 5;             % K-factor std [dB]
% "Decorrelation distances" in horizontal plane [m]
config.DS_lambda_LOS   = 7;          % [m], delay spread
config.AS_D_lambda_LOS = 8;          % [m], departure azimuth spread
config.AS_A_lambda_LOS = 8;          % [m], arrival azimuth spread
config.SF_lambda_LOS   = 10;         % [m], shadowing
config.KF_lambda_LOS   = 15;         % [m], k-factor
config.ZS_A_lambda_LOS = 12;         % [m], arrival zenith spread
config.ZS_D_lambda_LOS = 12;         % [m], departure zenith spread
 % CLuster angle spread value
config.cluster_ASD_LOS = 3;
config.cluster_ASA_LOS = 17;
config.cluster_ZSA_LOS = 7;     

% 3DUMi, NLOS     (3D Urban Micro)
% Fixed scenario specific parameters
config.NumClusters_NLOS = 19;         % Number of clusters    [Table 7.3-6]
config.r_DS_NLOS   = 3;               % Delays spread proportionality factor
config.PerClusterAS_D_NLOS = 10;      % Per cluster rms azimuth spread of departure angles [deg]
config.PerClusterAS_A_NLOS = 22;      % Per cluster rms azimuth spread of arrival angles [deg]
config.PerClusterZS_A_NLOS = 7;       % Per cluster rms zenith spread of arrival angles
config.LNS_ksi_NLOS = 3;              % LNS ksi [dB], per cluster shadowing std
% Cross correlation coefficients for azimuth case
config.asD_ds_NLOS = 0;               % departure AS vs delay spread
config.asA_ds_NLOS = 0.4;             % arrival AS vs delay spread
config.asA_sf_NLOS = -0.4;            % arrival AS vs shadowing std
config.asD_sf_NLOS = 0;               % departure AS vs shadowing std
config.ds_sf_NLOS  = -0.7;            % delay spread vs shadowing std
config.asD_asA_NLOS = 0;              % departure AS vs arrival AS
% Cross correlation coefficients for elevation (zenith) case
config.zsD_sf_NLOS = 0;               % departure ZS vs shadowing std
config.zsA_sf_NLOS = 0;               % arrival ZS vs shadowing std
config.zsD_ds_NLOS  = -0.5;           % departure ZS vs delay spread
config.zsA_ds_NLOS = 0;               % arrival ZS vs delay spread
config.zsD_asD_NLOS = 0.5;            % departure ZS vs departure AS
config.zsA_asD_NLOS = 0.5;            % arrival ZS vs departure AS
config.zsD_asA_NLOS = 0;              % departure ZS vs arrival AS
config.zsA_asA_NLOS = 0.2;            % arrival ZS vs arrival AS
config.zsD_zsA_NLOS = 0;              % departure ZS vs arrival ZS
% Polarisation parameters
config.xpr_mu_NLOS    = 8;            % XPR mean [dB]
config.xpr_sigma_NLOS = 3;            % XPR std  [dB]
% Dispersion parameters [1, Table 4.5]. Log-normal distributions
config.DS_mu_NLOS      = -6.89;       % delay spread, mean [log10(s)]
config.DS_sigma_NLOS   = 0.54;        % delay spread, std [log10(s)]
config.AS_D_mu_NLOS    = 1.41;        % azimuth departure angle spread, mean [log10(deg)]
config.AS_D_sigma_NLOS = 0.17;        % azimuth departure angle spread, std [log10(deg)]
config.AS_A_mu_NLOS    = 1.84;        % azimuth arrival angle spread, mean [log10(deg)]
config.AS_A_sigma_NLOS = 0.15;        % azimuth arrival angle spread, std [log10(deg)]
config.ZS_A_mu_NLOS    = 0.88;        % zenith arrival angle spread, mean [log10(deg)]
config.ZS_A_sigma_NLOS = 0.16;        % zenith arrival angle spread, std [log10(deg)]
config.SF_sigma_NLOS   = 4;           % shadowing std [dB] (zero mean)
% "Decorrelation distances" in horizontal plane [m]
config.DS_lambda_NLOS   = 10;         % [m], delay spread
config.AS_D_lambda_NLOS = 10;         % [m], departure azimuth spread
config.AS_A_lambda_NLOS = 9;          % [m], arrival azimuth spread
config.SF_lambda_NLOS   = 13;         % [m], shadowing
config.ZS_A_lambda_NLOS = 10;         % [m], arrival zenith spread
config.ZS_D_lambda_NLOS = 10;         % [m], departure zenith spread
% CLuster angle spread value
config.cluster_ASD_NLOS = 10;
config.cluster_ASA_NLOS = 22;
config.cluster_ZSA_NLOS = 7; 

%Added from tri_sector_3D line 436 to 484
%Change 'LTE_config' for config

 % 3DUMi, OTOI  (3D Urban Micro) Outdoor-TO-Indoor
% Fixed scenario specific parameters
config.NumClusters_OTOI = 12;         % Number of clusters    [Table 7.3-6]
config.r_DS_OTOI   = 2.2;             % Delays spread proportionality factor
config.PerClusterAS_D_OTOI = 5;       % Per cluster rms azimuth spread of departure angles [deg]
config.PerClusterAS_A_OTOI = 8;      % Per cluster rms azimuth spread of arrival angles [deg]
config.PerClusterZS_A_OTOI = 3;       % Per cluster rms zenith spread of arrival angles
config.LNS_ksi_OTOI = 4;                 % LNS ksi [dB], per cluster shadowing std
% Cross correlation coefficients for azimuth case
config.asD_ds_OTOI = 0.4;             % departure AS vs delay spread
config.asA_ds_OTOI = 0.4;             % arrival AS vs delay spread
config.asA_sf_OTOI = 0;            % arrival AS vs shadowing std
config.asD_sf_OTOI = 0.2;            % departure AS vs shadowing std
config.ds_sf_OTOI  = -0.5;            % delay spread vs shadowing std
config.asD_asA_OTOI = 0;            % departure AS vs arrival AS
% Cross correlation coefficients for elevation (zenith) case
config.zsD_sf_OTOI = 0;               % departure ZS vs shadowing std
config.zsA_sf_OTOI = 0;            % arrival ZS vs shadowing std
config.zsD_ds_OTOI  = -0.6;           % departure ZS vs delay spread
config.zsA_ds_OTOI = -0.2;               % arrival ZS vs delay spread
config.zsD_asD_OTOI = -0.2;            % departure ZS vs departure AS
config.zsA_asD_OTOI = 0;           % arrival ZS vs departure AS
config.zsD_asA_OTOI = 0;              % departure ZS vs arrival AS
config.zsA_asA_OTOI = 0.5;              % arrival ZS vs arrival AS
config.zsD_zsA_OTOI = 0.5;              % departure ZS vs arrival ZS
% Polarisation parameters
config.xpr_mu_OTOI    = 9;            % XPR mean [dB]
config.xpr_sigma_OTOI = 5;            % XPR std  [dB]
% Dispersion parameters [1, Table 4.5].Log-normal distributions
config.DS_mu_OTOI      = -6.62;       % delay spread, mean [log10(s)]
config.DS_sigma_OTOI   = 0.32;        % delay spread, std [log10(s)]
config.AS_D_mu_OTOI    = 1.25;        % azimuth departure angle spread, mean [log10(deg)]
config.AS_D_sigma_OTOI = 0.42;        % azimuth departure angle spread, std [log10(deg)]
config.AS_A_mu_OTOI    = 1.76;        % azimuth arrival angle spread, mean [log10(deg)]
config.AS_A_sigma_OTOI = 0.16;        % azimuth arrival angle spread, std [log10(deg)]
config.ZS_A_mu_OTOI    = 1.01;        % zenith arrival angle spread, mean [log10(deg)]
config.ZS_A_sigma_OTOI = 0.43;        % zenith arrival angle spread, std [log10(deg)]
config.SF_sigma_OTOI   = 7;           % shadowing std [dB] (zero mean)
% "Decorrelation distances" in horizontal plane [m]
config.DS_lambda_OTOI   = 10;         % [m], delay spread
config.AS_D_lambda_OTOI = 11;         % [m], departure azimuth spread
config.AS_A_lambda_OTOI = 17;         % [m], arrival azimuth spread
config.SF_lambda_OTOI   = 7;         % [m], shadowing
config.ZS_A_lambda_OTOI = 25;         % [m], arrival zenith spread
config.ZS_D_lambda_OTOI = 25;         % [m], departure zenith spread
% CLuster angle spread value
config.cluster_ASD_OTOI = 5;
config.cluster_ASA_OTOI = 8;
config.cluster_ZSA_OTOI = 3; 
end