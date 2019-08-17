% generate UE, BS ...
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%

cqi_i = LTE_params.cqi_i;

%% Generate UEs
UE = network_elements.UE(LTE_params.UE_config,LTE_params.Ntot,LTE_params.Nsub,LTE_params.BS_config.nRX,LTE_params.ChanMod_config.filtering,LTE_params.frameStructure,10);
for uu=1: (LTE_params.nUE * LTE_params.nBS) % first assumption: all user equipments have the same capabilities
    UE(uu) = network_elements.UE(LTE_params.UE_config,LTE_params.Ntot,LTE_params.Nsub,LTE_params.BS_config.nRX,LTE_params.ChanMod_config.filtering,LTE_params.frameStructure,10);
end


%% Generate eNodeBs
LTE_params.BS_config.AtPort = 0:LTE_params.BS_config.nRX-1;
LTE_params.BS_config.nAtPort = length(LTE_params.BS_config.AtPort);
LTE_params.BS_config.maxStreams = 2;
LTE_params.BS_config.HARQ_processes = LTE_params.HARQ_processes;

BS = network_elements.eNodeB(...
    1,...
    LTE_params.BS_config.nRX,...
    LTE_params.nUE,...
    LTE_params.BS_config.AtPort,...
    LTE_params.BS_config.nAtPort,...
    LTE_params.BS_config.maxStreams,...
    LTE_params.BS_config.HARQ_processes,...
    LTE_params.max_HARQ_retransmissions,...
    LTE_params.BS_config);
for bb=1:LTE_params.nBS % first assumption: all base stations have the same capabilities
    BS(bb) = network_elements.eNodeB(...
        bb,...
        LTE_params.BS_config.nRX,...
        LTE_params.nUE,...
        LTE_params.BS_config.AtPort,...
        LTE_params.BS_config.nAtPort,...
        LTE_params.BS_config.maxStreams,...
        LTE_params.BS_config.HARQ_processes,...
        LTE_params.max_HARQ_retransmissions,...
        LTE_params.BS_config);
end

%% update number of channel realizations
for bb=1:LTE_params.nBS
    BS(bb).realization_num_total = BS(bb).realization_num_total * BS(bb).nRX * ChanMod.nTX; %channel of every subframe consitcs of nRx*nTx different channel
end

%% Load or allocate channel autocorrelation matrix
for bb=1:LTE_params.nBS
     if strcmp(BS(bb).channel_estimation_method,'MMSE') || strcmp(BS(bb).channel_estimation_method,'MMSE_2D')
            switch ChanMod.interpolation_method
                case 'shift_to_nearest_neighbor'
                    switch LTE_params.ChanMod_config.type
                        case {'AWGN','flat Rayleigh'}
                              H = 1;
                        case {'winner_II', 'TR 36.873'}
                            error('MMSE channel estimation currently not supported for the 3D channel model or the winner II channel model.');
                        otherwise
                            NTap = size(ChanMod.PDP_dB,2);% Number of Taps
                            tap_delays = round(ChanMod.PDP_dB(2,:)*LTE_params.Fs);
                            H = zeros(tap_delays(end)+1,1);

                            for tap_i = 1:tap_delays(end)+1
                                H(tap_i) = sqrt(sum(10.^(ChanMod.PDP_dB(1,tap_delays == tap_i-1)./10)));
                            end
                            H = H./ChanMod.normH;    
                    end
                case 'sinc_interpolation'
                    switch LTE_params.ChanMod_config.type
                        case {'winner_II', 'TR 36.873'}
                            error('MMSE channel estimation currently not supported for the 3D channel model or the winner II channel model.');
                        otherwise
                            H = sum(ChanMod.interpolator.precomputed_sincs)./ChanMod.interpolator.norm_h;
                    end
                otherwise
                    error('channel interpolation method unkonwn');
            end
            DFT_matrix = dftmtx(LTE_params.Nfft);
            time_correlation_matrix = zeros(LTE_params.Nfft);
            time_correlation_matrix(1:length(H),1:length(H)) = diag(H.^2);
            % calculate perfect autocorrelation function R_ff = DFT_mat * R_tt * DFT_mat'
            autocorrelation_matrix_help = DFT_matrix * time_correlation_matrix * DFT_matrix';          
            % pick the correct values, remove DC ...
            autocorrelation_matrix = autocorrelation_matrix_help([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 1:LTE_params.Ntot/2],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 1:LTE_params.Ntot/2]);
            for uu=1:LTE_params.nUE*LTE_params.nBS
                UE(uu).channel_autocorrelation_matrix = autocorrelation_matrix;
            end
        if strcmp(BS(bb).autocorrelation_matrix_type,'estimated')
            for uu=1:LTE_params.nUE
                Data=load('cov_matrix.mat');
                UE(uu).channel_autocorrelation_matrix = Data.real_cov;
            end
        end
    end
end

%% Channel parameters dependent - now only the same channel parameters for each user and BS are allowed
% load Correlation Matrices
if strcmp(ChanMod.type,'PedA') || strcmp(ChanMod.type,'PedB') || strcmp(ChanMod.type,'VehA') || strcmp(ChanMod.type,'VehB') ||strcmp(ChanMod.type,'TU') || strcmp(ChanMod.type,'RA') || strcmp(ChanMod.type,'HT') || strcmp(ChanMod.type,'ePDP')
    ChanMod.corrRX = ones(size(ChanMod.PDP_dB,2),BS(1).nRX,BS(1).nRX);
    ChanMod.corrTX = ones(size(ChanMod.PDP_dB,2),UE(1).nTX,UE(1).nTX);
    for kk = 1:size(ChanMod.PDP_dB,2)
        ChanMod.corrRX(kk,:,:) = eye(BS(1).nRX);
        ChanMod.corrTX(kk,:,:) = eye(UE(1).nTX);
    end
elseif strcmp(ChanMod.type,'PedBcorr') || strcmp(ChanMod.type,'EVehA') || strcmp(ChanMod.type,'ETU')
    ChanMod.corrRX = ones(size(ChanMod.PDP_dB,2),BS(1).nRX,BS(1).nRX);
    ChanMod.corrTX = ones(size(ChanMod.PDP_dB,2),UE(1).nTX,UE(1).nTX);
    for kk = 1:size(ChanMod.PDP_dB,2)
        ChanMod.corrRX(kk,:,:) = eye(BS(1).nRX) + ChanMod.corr_coefRX*ones(BS(1).nRX) - ChanMod.corr_coefRX*eye(BS(1).nRX);
        ChanMod.corrTX(kk,:,:) = eye(UE(1).nTX) + ChanMod.corr_coefTX*ones(UE(1).nTX) - ChanMod.corr_coefTX*eye(UE(1).nTX);
    end
end

%% generate psi, theta, phi for Rosa-Zheng Channel Model
switch ChanMod.time_correlation 
    case 'correlated'
    switch ChanMod.type
        case {'PedA', 'PedB', 'PedBcorr','VehA','VehB','TU','RA','HT','EPedA','EVehA','ETU','ePDP'}
            for uu=1:(LTE_params.nUE*LTE_params.nBS)
                
                bb = utils.findBS( LTE_params.connection_table, uu );

                number_of_taps = ChanMod.interpolator.num_faders;

                UE(uu).channel_coef_rosa_zheng.theta = (rand(LTE_params.channel_param_RandStream,BS(bb).nRX,UE(uu).nTX,number_of_taps)*2 -1) * pi;
                UE(uu).channel_coef_rosa_zheng.phi   = (rand(LTE_params.channel_param_RandStream,BS(bb).nRX,UE(uu).nTX,number_of_taps,ChanMod.sin_num)*2 -1) * pi; % The Original paper states \phi as not changing for each sinusoid, but that is not the case (see T. Zemen and C.F. Mecklenbr�auker, �Time-Variant Channel Estimation Using Discrete Prolate Spheroidal Sequences,� IEEE Transactions on Signal Processing, vol. 53, no. 9, pp. 3597�3607, Sept. 2005.)
                UE(uu).channel_coef_rosa_zheng.psi   = (rand(LTE_params.channel_param_RandStream,BS(bb).nRX,UE(uu).nTX,number_of_taps,ChanMod.sin_num)*2 -1) * pi;

            end
        otherwise
    end
    case 'independent'
end

%% Scheduler initial initialization

LTE_params.scheduler.PMI = zeros(LTE_params.Nrb,2);

% if LTE_params.downlink_delay==0
%     LTE_params.scheduler.zero_delay = true;
% else
%     LTE_params.scheduler.zero_delay = false;
% end

LTE_params.scheduler.UE_specific = BS(1).UE_specific; 

% Copy CQI mapping configuration to the scheduler to allow it to call the CQI mapping funciton with correct parameters
LTE_params.scheduler.CQI_mapping_params = LTE_params.CQI_mapping;

LTE_params.MI_data = network_elements.MIdata_load;   %added for LTE_feedback_getBICM.m

%% Further scheduler initialization: Create scheduler. Assume the eNodeB where the scheduler is to be put is the first one
LTE_params.scheduler.max_HARQ_retx = LTE_params.max_HARQ_retransmissions;

switch LTE_params.UE_config.SINR_averaging.averager
    case 'EESM'
        averager        = network_elements.eesmAverager(LTE_params.UE_config.SINR_averaging.EESMbetas,LTE_params.UE_config.SINR_averaging.MCSs);
    case 'MIESM'
        averager        = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,LTE_params.UE_config.SINR_averaging.MCSs);
    case 'HARM_MEAN'
        averager        = network_elements.harmMeanAverager;     
    otherwise
        error('SINR averager not supported');
end



for bb = 1:LTE_params.nBS
    
    bb_users = LTE_params.connection_table(bb,:) == 1; 
    
    switch LTE_params.scheduler.type
        case 'fixed' 
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.fixedScheduler(LTE_params,UE(bb_users),averager,mapping_data);
        case 'fixed MU MIMO'
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.fixed_MU_MIMO(LTE_params,UE(bb_users),averager,mapping_data);
        case 'max MU MIMO'
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.optMaxMUMIMO(LTE_params,UE(bb_users),averager,mapping_data);
        case 'random MU MIMO'
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.randomMUMIMO(LTE_params,UE(bb_users),averager,mapping_data);
        case 'greedy MU MIMO'
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.greedyMUMIMO(LTE_params,UE(bb_users),averager,mapping_data);
        case 'greedy frobenius MU MIMO'
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.frobMUMIMO(LTE_params,UE(bb_users),averager,mapping_data);
        case 'round robin' 
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.roundRobinScheduler(LTE_params,UE(bb_users),averager,mapping_data);
        case 'opt max throughtput'
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.optMaxThroughputScheduler(LTE_params,UE(bb_users),averager,mapping_data);
        case 'approx max throughtput'
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.approxMaxThroughputScheduler(LTE_params,UE(bb_users),averager,mapping_data);
        case 'best CQI'
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.bestCQIScheduler(LTE_params,UE(bb_users),averager,mapping_data);
        otherwise
            error('scheduler not supported');
    end
end

