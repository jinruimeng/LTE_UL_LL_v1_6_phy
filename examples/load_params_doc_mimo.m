% Configuration file for uplink
% based on ver. 1.6
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% Jan Prokopec,
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% Contains code the following sources:
%   - Coded Modulation Library (CML) by Iterative Solutions


global LTE_params;
LTE_params.cqi_i = cqi_i;
LTE_params.SNR_vec = SNR_vec;

%% General
LTE_params.Simulation_type = 'uplink';      % 'uplink'
LTE_params.simulation_method = 'parallel';    % 'parallel' or 'normal' to parallelize the SNR loop in LTE_sim_main_parallel.m
LTE_params.nUE = 1;                         % number of user equipments to simulate
LTE_params.nBS = 1;                         % number of base stations to simulate (hard-coded to 1)
LTE_params.Bandwidth = 1.4e6;               % in Hz, allowed values: 1.4 MHz, 3 MHz, 5 MHz, 10 MHz, 15 MHz, 20MHz => number of resource blocks 6, 15, 25, 50, 75, 100

%% PRIMITIVE PARAMETERS
%DEFINED IN STANDARD 3GPP TS 36.104 V8.1.0 (2008-03), page 10
%Table 5.1-1 Transmission bandwidth configuration
LTE_params.HARQ_processes = 8;                          % Number of HARQ processes
LTE_params.max_HARQ_retransmissions = 0;                % max num of HARQ retransmissions, NOT including the original tx. 0, 1, 2 or 3
LTE_params.SubcarrierSpacing = 15e3;                    % in Hz, 15 kHz only for uplink
LTE_params.CyclicPrefix = 'normal';                     % 'normal' or 'extended'
LTE_params.specialSubframeConf = 0;
LTE_params.introduce_frequency_offset =  false;         % 'false', not implemented yet    
LTE_params.Nslot = 2;                                   % number of slots in one subframe
LTE_params.Nsfr = 10;                                   % number of subframes in one frame
LTE_params.downlink_delay = 0;                          % Delay the downlink channel will introduce (in TTIs), only 0 implemented for now
LTE_params.trace_subcarrier_SNR = false;                % 'false', not implemented
LTE_params.N_seed_reset = 1;                            % resets the random number generator seeds to a new value after LTE_params.N_seed_reset subframes
LTE_params.carrier_freq_UP = 1.9e9;                     % UPLINK carrier frequency [Hz]
LTE_params.speed_of_light = 299792458;                  % [m/s]

LTE_params.DFT_spreading_off = false;                   % turn off DFT spreading to have OFDM performance


%% Plots
LTE_params.show_plots = false;
LTE_params.plot_confidence_intervals = true;

%% Define User parameters (identical settings).
LTE_params.UE_config.mode = 4;                      % DEFINED IN STANDARD 3GPP TS 36.213 V11.4.0 Section 8.0, page 105
                                                    % 1: Single-antenna port
                                                    % 4: Closed-loop spatial multiplexing
LTE_params.UE_config.nTX = 2;                       % number of UE transmit antennas {1,2,4}
LTE_params.UE_config.user_speed = 0;          % [m/s]
LTE_params.UE_config.N_soft = 1e6*LTE_params.HARQ_processes;    % Defines the total number of soft channel bits available for HARQ processing (TS 36.306 4.2.1.3). Set to a high enough value.


%% Define BS parameters
LTE_params.BS_config.channel_estimation_method = 'PERFECT';         % 'PERFECT','LS_AV','LS_SAV', 'LS_QS', 'DFT'
LTE_params.BS_config.channel_interpolation_method = 'linear';       % 'linear','cubic','spline'
LTE_params.BS_config.channel_interpolation_past_points = 0;         % number of previous channel estimates (slots) used for interpolation {0,1,2}
LTE_params.BS_config.channel_est_time_smoothing = 0;                % time smoothing contraint, currently not in use
LTE_params.BS_config.channel_est_frequency_smoothing = 8;           % frequency smoothing constraint for LS_QS channel estimation, 8
LTE_params.BS_config.nRX = 4;                                       % number of BS receive antennas
LTE_params.BS_config.receiver = 'ZF';                               % 'ZF', 'MMSE'
LTE_params.BS_config.turbo_iterations = 8;                          % Number of iterations of the turbo decoder
LTE_params.BS_config.autocorrelation_matrix_type = 'estimated';         % 'ideal', 'estimated' type of autocorrelation matrix
LTE_params.BS_config.realization_num = 0;                           % Number of realizations of channel, used for averaging fo channel autocorrelation matrix
LTE_params.BS_config.realization_num_total = 20;                    % First xy number of channel realizations are used just for estimation of autocorrelation matrix

%% RS higher layer configuration
LTE_params.cyclicShift = 0;                 % {0 to 7}
LTE_params.CSFiDCIFormat = 000;             % Cyclic Shift Field in DCI format 0
LTE_params.activate_DMRS_with_OCC = true;   % activate OCC, 36.211 section 5.5.2.1.1
LTE_params.deltass = 0;                     % {0..29}
LTE_params.groupHopping = true;             % enable sequence group hopping
LTE_params.sequenceHopping = false;         % enable sequence hopping
LTE_params.nID_PUSCH = NaN;                 % can be provided, otherwise NAN        
LTE_params.nID_PUCCH = NaN;                 % can be provided, otherwise NAN 
LTE_params.srsSubframeConfig = 0;           % {0..15} 15..SRS off
LTE_params.Csrs = 1;                        % srs-BandwidthConfig {0..7}
LTE_params.srsMaxUpPTS = false;             % {true, false}
LTE_params.frameStructure = 1;              % '1', FDD only
LTE_params.Nsp = 2;                         % number of uplink to downlink switching points {1,2}

%% higher layer SRS configuration
LTE_params.UE_config.Isrs = 0;                  % srs configuration index {0..1023}
LTE_params.UE_config.SRS_triggerType = 0;       % {0,1}
LTE_params.UE_config.SRS_duration = 'periodic'; % 'periodic' 
LTE_params.UE_config.nCS_SRS = 4;               % cyclicShift {0..7}
LTE_params.UE_config.Bsrs = 2;                  % Bandwidth {0..3}
LTE_params.UE_config.b_hop = 1;                 % HoppingBandwidth {0..3}
LTE_params.UE_config.nRRC = 15;                 % freqDomainPosition {0..31}
LTE_params.UE_config.kTC = 1;                   % transmission comb {0,1}

%% Define ChanMod parameters - now it is only possible to have same channel parameters for BS and UE
LTE_params.ChanMod_config.filtering = 'BlockFading';                             % 'BlockFading', 'FastFading'
LTE_params.ChanMod_config.type = 'TU';                                          % 'PedA', 'PedB', 'PedBcorr', 'AWGN', 'flat Rayleigh','VehA','VehB','TU','RA','HT','winner_II', 'ePDP'
LTE_params.ChanMod_config.interpolation_method = 'shift_to_nearest_neighbor';   % 'shift_to_nearest_neighbor' 'sinc_interpolation'
LTE_params.ChanMod_config.time_correlation = 'independent';                      % 'independent', 'correlated'
LTE_params.ChanMod_config.corr_coefRX = 0;                                      % correlation between receive antennas
LTE_params.ChanMod_config.corr_coefTX = 0;                                      % correlation between transmit antennas
LTE_params.ChanMod_config.sin_num = 10;                                         % Number of sin realizations
%% Parameter of exponential PDP
LTE_params.ChanMod_config.tau_rms = 1.5e-6;
%% Channel matrix source
LTE_params.channel_matrix_source     = 'generated';     % 'generated' to generate every time, 'trace' to load it from a trace
LTE_params.channel_matrix_tracefile  = 'auto';          % filename where the trace is stored, only if 'trace' mode is used
LTE_params.store_channel_trace       = false;           % if mode is 'generated', the channel trace will be saved at the end of the simulation
%% winner model settings
LTE_params.ChanMod_config.winner_settings.Scenario                 = 11;        % 1=A1, 2=A2, 3=B1, 4=B2, 5=B3, 6=B4, 7=B5a, 8=B5c, 9=B5f, 10=C1,
                                                                                % 11=C2, 12=C3, 13=C4, 14=D1 and 15=D2a
LTE_params.ChanMod_config.winner_settings.PropagCondition          = 'NLOS';    % [LOS,{NLOS}]
LTE_params.ChanMod_config.winner_settings.SampleDensity            = 2;         % number of time samples per half wavelength [ {2} ]
LTE_params.ChanMod_config.winner_settings.UniformTimeSampling      = 'yes';     % Use same time sampling grid for all links [ yes | {no} ]
LTE_params.ChanMod_config.winner_settings.FixedPdpUsed             = 'no';      % nonrandom path delays and powers [ yes | {no}]
LTE_params.ChanMod_config.winner_settings.FixedAnglesUsed          = 'no';      % nonrandom AoD/AoAs [ yes | {no} ]
LTE_params.ChanMod_config.winner_settings.PolarisedArrays          = 'yes';     % usage of dual polarised arrays [ {yes} | no ]
LTE_params.ChanMod_config.winner_settings.TimeEvolution            = 'no';      % usage of time evolution  [ yes | {no} ]
LTE_params.ChanMod_config.winner_settings.PathLossModelUsed        = 'no';      % usage of path loss model [ yes | {no} ]
LTE_params.ChanMod_config.winner_settings.ShadowingModelUsed       = 'no';      % usage of shadow fading model [ yes | {no} ]
LTE_params.ChanMod_config.winner_settings.PathLossModel            = 'pathloss';% path loss model function name [ {pathloss} ]
LTE_params.ChanMod_config.winner_settings.PathLossOption           = 'CR_light';% ['{CR_light}' or 'CR_heavy' or 'RR_light' or 'RR_heavy', CR = Corridor-Room, RR = Room-Room nlos}
LTE_params.ChanMod_config.winner_settings.RandomSeed               = 27;        % sets random seed [ {[empty]} ]
LTE_params.ChanMod_config.winner_settings.UseManualPropCondition   = 'yes';     % whether to use manual propagation condition (los/nlos) setting or not. If not, the propagation condition is drawn from probabilities.  [ {yes} | no]

%% Scheduler settings
LTE_params.scheduler.type = 'round robin';
% Available options are:
%   - 'round robin'                             : Will serve equally all of the available users
%   - 'best cqi'                                : Will serve only users that maximize the CQI for specific RB
%   - 'best cqi traffic'                        : Best CQI adapted to traffic models
%   - 'proportional fair'
%   - 'proportional fair traffic'               : PF adapted to traffic models
%   - 'proportional fair RT-Priority traffic'   : Real Time users are prioritized over Non_Real_time users like http and ftp. The spare RB  are assignes with the normal PF metric.
%   - 'proportional fair criteria traffic'
%   - 'proportional fair exp/pf traffic'
%   - 'max throughput'
%   - 'max min'
%   - 'resource fair'
%   - 'Kwan'
%   - 'opt throughput'                          : Only implemented for nUE=2
    
LTE_params.scheduler.av_window = 10;            % Scheduler averaging window size in subframes (not supported by all schedulers)
LTE_params.scheduler.fairness = 0.8;            % desired fairness for the variable fairness scheduler (only supported by the var fair scheduler)

%% Traffic Models
LTE_params.trafficmodel.usetraffic_model = false;   % for simulation (true) or fullbuffer assumption (false) (note: this is just used if the scheduler supports it!)

%% Initial state of random number generator
LTE_params.simulate_with_all_zero_sequences = false;% true if you want that the transmitted data is an all-zero sequence (useful for interleaver testing)
LTE_params.random_data_seeding = true;              % Whether the seed for the random number generator that generates the transmitted databits is set
LTE_params.data_seed = 10;                          % Only used if the upper variable is set to 'true'
LTE_params.random_channel_param_seeding = true;     % Whether the seed for the random number generator that generates the channel parameters is set
LTE_params.channel_param_seed = 175; % 3            % Only used if the upper variable is set to 'true'

%% Configuration of the random noise generation
LTE_params.random_noise_seeding = false;            % Whether the seed for the random number generator that generates the noise is set
LTE_params.noise_seed = 10;                         % Only used if the upper variable is set to 'true'

%% CQI mapping
LTE_params.CQI_mapping.coeffs = [0.5223 4.6176];    % these values hold for BLER 0.1 at first transmission
LTE_params.CQI_mapping.table = [-500;-6.4;-4.8;-2.95;-1.1;0.9;2.70;4.85;6.65;8.5;10.45;12.5;14.25;15.9;17.95;19.4;21];  % UL BLER 0.1 at first transmission

%% Feedback Parameters
LTE_params.UE_config.PMI_fb = true;                             % PMI feedback activated (just used  with LTE_params.UE_config.mode = 4)
LTE_params.UE_config.RI_fb = true;                              % RI feedback activated (just used with LTE_params.UE_config.mode = 4)
LTE_params.UE_config.CQI_fb = true;                             % CQI feedback activated (CQI from cqi vector in batch file is used when this is turned off)
LTE_params.UE_config.PMI  = 0;                                  % used when PMI_fb = false
LTE_params.UE_config.RI  = 2;                                   % used when RI_fb = false
LTE_params.UE_config.predict = false;                           % channel prediction activated (used in the feedback calculation)
LTE_params.UE_config.ignore_channel_estimation = false;          % ignores the channel estimation MSE for the feedback calculation if set
LTE_params.UE_config.ignore_ISI_ICI = false;                     % ignores the introduced ISI and ICI (dependend on the channel model)
LTE_params.UE_config.channel_averaging = true;                   % use channel averaging for feedback calculation
LTE_params.UE_config.CQI_fb_granularity = 'coarse';             % CQI feedback granularity, 'coarse','fine' or number in multiples of RBs

LTE_params.UE_config.SINR_averaging.averager = 'HARM_MEAN';     % SINR averaging method used for the feedback calculation 'EESM','MIESM','HARM_MEAN'
LTE_params.UE_config.SINR_averaging.EESMbetas = [5,5.01,5.01,0.84,1.67,1.61,1.64,3.87,5.06,6.4,12.59,17.59,23.33,29.45,33.05,35.41];    % weigthed EESM beta values
LTE_params.UE_config.SINR_averaging.MIESMbetas = [4,3.07,4.41,0.6,1.16,1.06,1.06,0.87,1.01,1.04,1.03,1.11,1.01,1.07,1,1.05];            % weighted MIESM beta values
LTE_params.UE_config.SINR_averaging.MCSs = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];


