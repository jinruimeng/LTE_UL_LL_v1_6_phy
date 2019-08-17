% Basic batch simulation script
% uplink, based on ver. 1.6
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
% www.urel.feec.vutbr.cz


close all; 
clc;
clear classes;
%clearvars;
%clear all
%warning off

%% DEBUG level
global DEBUG_LEVEL 
DEBUG_LEVEL = 1;    % 0: no text output at all
                    % 1: shows current simulation status
                    % 2: shows additional variables and figures

global LTE_params;

%% Actual simulations
cqi_vec = [1];      % only used for static scheduling
SNR_vec= [1000];
sim_result_temp.UE_specific.throughput_coded =[];
user_speed=1;

for cqi_i = cqi_vec
  
    N_subframes =2;%number of subframe
    
    % load Simulation Parameters
    LTE_UL_load_parameters;  
    % test data intialization
   
    
    % Load dependent Parameters and generate Elements
    LTE_UL_load_parameters_dependent;
    LTE_UL_load_parameters_generate_elements;
    LTE_UL_check_parameters;
    % actual Simulation
    tic
    LTE_UL_sim_main;
    toc    

end



save(fullfile('results', [output_filename_UL '.mat']));