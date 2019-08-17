% Generate AWGN BLER curves
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

clear
clear global
close all
clc

cd ..

% true if you are simulating the LTE CQIs. otherwise, a simulation with the
% MCSs defined in R1-071967 is employed
LTE_CQIS_simulation = true;

% Just plot, no simulations (for the case you already have the results)
skip_simulations = true;

% the BLER cures can also be output as CSV (comma-separated values) files
output_CSVs      = false;

%% DEBUG level
global DEBUG_LEVEL;
DEBUG_LEVEL = 1; 

%% SNR setting
if LTE_CQIS_simulation
    SNR_begin =  [-11 -8 -6 -3.5 -1.5 1   3   5  7  9 11 12.5 14.5 16 18];
    SNR_end   =  [-4 -3 -1  1    2.5 4.5 6.5 8 10 12 14 16   17.5 20 22];   % uplink
    %SNR_end   = [-5 -4 -2  0    1.5 3.5 5.5 7  9 11 13 15   16.5 19 21];   % downlink
else
    %SNR_begin = [-7 -6 -5 -4 -3 -2 -1 0 2 3 4 4 5 6  6   7  7  8 10 11 12 12 13 14 15 16 17];  % downlink
    SNR_begin =  [-8 -6 -5 -4 -3 -2 -1 0 2 3 4 4 5 6  6   7  7  8 10 11 12 12 13 14 15 16 17];   % uplink
    SNR_end   =  [-3 -2 -1  0  1  2  3 4 6 7 8 8 9 10 10 11 11 12 14 15 16 16 17 18 19 20 21];
end
SNR_end   = SNR_end + 0.5;
SNR_stepsize = 0.25;

if LTE_CQIS_simulation
    cqi_list = 1:15;
    results_folder = './examples/AWGN_BLER_curves';
else
    cqi_list = 101:127;
    results_folder = './examples/R1-071967';
end
BW_vect = 1.4e6;



%% Actual simulations
if ~skip_simulations
    for set_BW = BW_vect
        for cqi_idx = 1:length(cqi_list)
            cqi_i       = cqi_list(cqi_idx);
            N_subframes = 5000;
            SNR_vec     = SNR_begin(cqi_idx):SNR_stepsize:SNR_end(cqi_idx);
            
            LTE_UL_load_parameters;
            %% parameters for SUSISO 
            LTE_params.simulation_method = 'parallel';  % 'parallel' or 'normal' to parallelize the SNR loop in LTE_sim_main_parallel.m
            LTE_params.nUE = 1;                         % number of user equipments to simulate
            LTE_params.nBS = 1;                         % number of base stations to simulate (hard-coded to 1)
            LTE_params.Bandwidth = set_BW;              % in Hz, allowed values: 1.4 MHz, 3 MHz, 5 MHz, 10 MHz, 15 MHz, 20MHz => number of resource blocks 6, 15, 25, 50, 75, 100
            LTE_params.max_HARQ_retransmissions = 0;    % max num of HARQ retransmissions, NOT including the original tx. 0, 1, 2 or 3
            LTE_params.CyclicPrefix = 'normal';         % 'normal' or 'extended'
            LTE_params.downlink_delay = 0;              % Delay the downlink channel will introduce (in TTIs), only 0 implemented for now
            LTE_params.show_plots = false;
            LTE_params.UE_config.mode = 1;              % 1: Single-antenna port
            LTE_params.UE_config.nTX = 1;               % number of UE transmit antennas {1,2,4}
            LTE_params.UE_config.user_speed = 0/3.6;    % [m/s]

            LTE_params.BS_config.channel_estimation_method = 'PERFECT';         %'PERFECT','LS','MMSE'
            LTE_params.BS_config.channel_interpolation_method = 'linear';       %'linear','cubic','spline','sinc_freq','sinc_time','T-F'
            LTE_params.BS_config.nRX = 1;                                       % number of BS receive antennas
            LTE_params.BS_config.receiver = 'ZF';                               % 'SSD','ZF'
            LTE_params.ChanMod_config.filtering = 'BlockFading';                            % 'BlockFading','FastFading'
            LTE_params.ChanMod_config.type = 'AWGN';                                       
            LTE_params.ChanMod_config.time_correlation = 'independent';                      % 'independent', 'correlated'
            LTE_params.ChanMod_config.interpolation_method = 'shift_to_nearest_neighbor';   % 'shift_to_nearest_neighbor' 'sinc_interpolation'

            LTE_params.scheduler.type = 'round robin';
                        
            LTE_params.UE_config.PMI_fb = false;                             % PMI feedback activated (just used  with LTE_params.UE_config.mode = 4)
            LTE_params.UE_config.RI_fb = false;                              % RI feedback activated (just used with LTE_params.UE_config.mode = 4)
            LTE_params.UE_config.CQI_fb = false;                             % CQI feedback activated (CQI from cqi vector in batch file is used when this is turned off)
            LTE_params.UE_config.PMI  = 0;                                  % used when PMI_fb = false
            LTE_params.UE_config.RI  = 1;                                   % used when RI_fb = false
            LTE_params.UE_config.predict = false;                           % channel prediction activated (used in the feedback calculation)
            LTE_params.UE_config.ignore_channel_estimation = true;          % ignores the channel estimation MSE for the feedback calculation if set
            LTE_params.UE_config.ignore_ISI_ICI = true;                     % ignores the introduced ISI and ICI (dependend on the channel model)
            LTE_params.UE_config.channel_averaging = true;                   % use channel averaging for feedback calculation
            LTE_params.UE_config.CQI_fb_granularity = 'coarse';             % CQI feedback granularity, 'coarse','fine' or number in multiples of RBs
            
            % Load dependent Parameters and generate Elements
            LTE_UL_load_parameters_dependent;
            LTE_UL_load_parameters_generate_elements;
            LTE_UL_check_parameters;
            
            LTE_UL_sim_main;
            
            % Code to generate the output filename
            output_filename = sprintf('cqi_%.0f_%.1fMHz',cqi_i,set_BW/1e6);
            filename_suffix = [];
            
            save(fullfile(results_folder,[output_filename filename_suffix '.mat']));
        end
    end
    
    %% Compact results (too big to put them into the release, if not)
    for set_BW = BW_vect
        for cqi = cqi_list
            filename_suffix         = [];
            output_filename         = sprintf('cqi_%.0f_%.1fMHz',cqi,set_BW/1e6);
            output_filename_compact = sprintf('cqi_%.0f_%.1fMHz_compact',cqi,set_BW/1e6);
            
            full_filename         = fullfile(results_folder,[output_filename filename_suffix '.mat']);
            full_filename_compact = fullfile(results_folder,[output_filename_compact filename_suffix '.mat']);
            
            utils.resultsFileReader.delete_results_file_info_except_throughput_and_BLER(full_filename,full_filename_compact);
        end
    end
end

%% Plot files
%LTE_load_parameters;
%[LTE_params,~,~] = LTE_load_parameters_dependent(LTE_params,1);
cqi_i = 1;
SNR_vec = 1;
N_subframes = 1;
LTE_UL_load_parameters;
LTE_UL_load_parameters_dependent;

for set_BW = BW_vect
    for cqi = cqi_list
        % Code to generate the output filename
        output_filename = sprintf('cqi_%.0f_%.1fMHz',cqi,set_BW/1e6);
        filename_suffix = [];
        
        full_filename = fullfile(results_folder,[output_filename filename_suffix '.mat']);
        
        % Load data
        result_files{cqi}    = utils.resultsFileReader(full_filename);
        sim_summary(cqi)     = result_files{cqi}.get_config_params_summary;
        SNR_vect{cqi}        = result_files{cqi}.get_SNR_vector;
        cell_throughput{cqi} = result_files{cqi}.get_throughput_over_SNR;
        cell_BLER{cqi}       = result_files{cqi}.get_UE_BLER_over_SNR;
    end
    
    all_sims_lengths = unique([sim_summary.simulation_length]);
    if length(all_sims_lengths)==1
        all_sims_lengths_string = sprintf('%g',all_sims_lengths);
    else
        all_sims_lengths_string = sprintf('%g-%g',min(all_sims_lengths),max(all_sims_lengths));
    end
    
    % Plot BLER data
    figure;
    axes('YScale','log');
    hold all;
    for cqi = cqi_list
        modulation  = LTE_params.CQI_params(cqi).modulation;
        coding_rate = LTE_params.CQI_params(cqi).coding_rate_x_1024/1024;
        if LTE_CQIS_simulation
            CQI_name = sprintf('CQI %.0f',cqi);
            CSV_filename = fullfile(results_folder,sprintf('BLER_CQI%0.f.csv',cqi));
        else
            CQI_name = sprintf('MCS %.0f: %s, %3.2f',cqi-100,modulation,coding_rate);
            CSV_filename = fullfile(results_folder,sprintf('BLER_MCS%0.f.csv',cqi-100));
        end
        plot(SNR_vect{cqi},cell_BLER{cqi},'DisplayName',CQI_name);
        if output_CSVs
            csvwrite(CSV_filename,[SNR_vect{cqi}(:),cell_BLER{cqi}(:)].')
        end
    end
    ylim([1e-3 1]);
    grid on;
    legend('show','Location','SouthEastOutside');
    title(sprintf('AWGN BLER: %.1f MHz, %s subframes',set_BW/1e6,all_sims_lengths_string));
    hgsave(fullfile(results_folder,sprintf('plot_BLER_%0.1fMHz.fig',set_BW/1e6)));
    
    % Plot throughput data
    figure;
    hold all;
    for cqi = cqi_list
        modulation  = LTE_params.CQI_params(cqi).modulation;
        coding_rate = LTE_params.CQI_params(cqi).coding_rate_x_1024/1024;
        if LTE_CQIS_simulation
            CQI_name = sprintf('CQI %.0f',cqi);
            CSV_filename = fullfile(results_folder,sprintf('throughput_CQI%0.f.csv',cqi));
        else
            CQI_name = sprintf('MCS %.0f: %s, %3.2f',cqi-100,modulation,coding_rate);
            CSV_filename = fullfile(results_folder,sprintf('throughput_MCS%0.f.csv',cqi-100));
        end
        plot(SNR_vect{cqi},cell_throughput{cqi},'DisplayName',CQI_name);
        if output_CSVs
            csvwrite(CSV_filename,[SNR_vect{cqi}(:),cell_BLER{cqi}(:)].')
        end
    end
    grid on;
    legend('show','Location','SouthEastOutside');
    title(sprintf('AWGN throughput: %.1f MHz, %s subframes',set_BW/1e6,all_sims_lengths_string));
    hgsave(fullfile(results_folder,sprintf('plot_throughput_%0.1fMHz.fig',set_BW/1e6)));
end

