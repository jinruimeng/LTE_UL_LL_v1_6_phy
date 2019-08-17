% Check parameter consistency
% Author: Josep Colom, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

SNR_vec = LTE_params.SNR_vec;

if N_subframes < 3 && LTE_params.confidence_interval_probability > 0
    warning('At least 3 subframes are required for confidence intervals.');
end


% bandwidth  1.4 MHz, 3 MHz, 5 MHz, 10 MHz, 15 MHz, 20MHz
if~((LTE_params.Bandwidth==1.4e6) || (LTE_params.Bandwidth==3e6) || (LTE_params.Bandwidth==5e6) || (LTE_params.Bandwidth==10e6) || (LTE_params.Bandwidth==15e6) || (LTE_params.Bandwidth==20e6))
	error('selected bandwidth is not standardized');
end

% plotting
% if length(SNR_vec)==1 && LTE_params.show_plots
%     error('plots are not usefull for only one SNR point');
% end

% subcarrier spacing
if(LTE_params.SubcarrierSpacing ~= 15e3)
    error('in the LTE-A uplink, only 15kHz subcarrier spacing is allowed');
end

% CP length
if(strcmp(LTE_params.CyclicPrefix,'extended') && LTE_params.SubcarrierSpacing == 7.5e3)
    error('For this combination of subcarrier spacing and cyclic prefix (MBSFN transmissions) reference and synchronization channels not implemented');
end

% channel estimation
if (LTE_params.BS_config.channel_interpolation_past_points > 0) && strcmp(LTE_params.ChanMod_config.time_correlation,'independent')
    error('interpolation exploiting previous channel estimates only useful for correlated channel realizations');
end

if (LTE_params.downlink_delay > 0) && (LTE_params.nUE > 1)
    error('Simulations with downlink feedback channel delay only supported with a single user.');
end

if strcmp(LTE_params.UE_config.MCS_and_scheduling_CSI, 'estimated') && (LTE_params.downlink_delay > 0) && strcmp(LTE_params.ChanMod_config.filtering, 'BlockFading')
    error('For LTE_params.UE_config.MCS_and_scheduling_CSI = ''estimated'' and LTE_params.downlink_delay > 0 only LTE_params.ChanMod_config.filtering = ''FastFading'' is supported to perform channel prediction.');
end


if strcmp(LTE_params.BS_config.channel_estimation_method,'PERFECT') && strcmp(LTE_params.UE_config.MCS_and_scheduling_CSI, 'estimated')
    error('Estimated ''MCS and scheduling CSI'' only available with non perfect channel estimation method.');
end

if strcmp(LTE_params.BS_config.channel_estimation_method,'PERFECT') && LTE_params.BS_config.channel_prediction
    error('Channel prediction is only supported in combination with channel estimation.');
end

if (LTE_params.BS_config.channel_interpolation_past_points > 0) && (LTE_params.nUE > 1)
    error('Exploiting previous points for channel interpolation only supported with a single user.');
end

if (LTE_params.BS_config.channel_interpolation_past_points > 0) && strcmp(LTE_params.ChanMod_config.filtering, 'BlockFading')
    error('On a block-fading channel previous points for channel interpolation are not considered as no interpolation is performed.');
end

if LTE_params.BS_config.channel_prediction && (LTE_params.UE_config.mode == 5)
    error('channel prediction not supported for MU-MIMO');
end

if LTE_params.BS_config.channel_prediction && (LTE_params.nUE > 1)
    error('Channel prediction only supported with a single user.');
end

if LTE_params.BS_config.channel_prediction && strcmp(LTE_params.ChanMod_config.time_correlation, 'independent')
    error('Channel prediction only useful for a time correlated channel.');
end

if LTE_params.BS_config.channel_prediction && strcmp(LTE_params.ChanMod_config.filtering, 'BlockFading')
    warning('Channel prediction on a BlockFading channel will probably result in a high prediction error.');
end

if strcmp(LTE_params.BS_config.channel_estimation_method,'MMSE_2D') && strcmp(LTE_params.ChanMod_config.filtering, 'BlockFading')
    error('The 2D MMSE channel estimator is only defined for FastFading');
end

if (strcmp(LTE_params.BS_config.channel_estimation_method,'MMSE_2D') && ~strcmp(LTE_params.BS_config.channel_interpolation_method, 'MMSE_2D')) || ...
   (~strcmp(LTE_params.BS_config.channel_estimation_method,'MMSE_2D') && strcmp(LTE_params.BS_config.channel_interpolation_method, 'MMSE_2D'))
    error('To use 2D MMSE estimation please set estimation and interpolation method to MMSE_2D.');
end

% HARQ
if LTE_params.max_HARQ_retransmissions > 3 || LTE_params.max_HARQ_retransmissions < 0
    error('Maximum HARQ retransmissions cannot be higher than 3 or negative');
end
if LTE_params.HARQ_processes > 8
    error('The standard does not allow more than 8 HARQ processes. Deleting this error message may result in errors in rate matching process.');
end

% check number of transmit antennas
for bb=1:LTE_params.nBS
    if (BS(bb).nRX == 3 || BS(bb).nRX < 1)
        error('number of antennas not supported');
    end
    if (BS(bb).nAtPort ~= BS(bb).nRX)
        error('number of antenna ports not consistent with number of antennas');
    end
end

% channel model
if (LTE_params.UE_config.user_speed > eps) && (strcmp(LTE_params.ChanMod_config.time_correlation,'independent'))...
        && strcmp(LTE_params.ChanMod_config.filtering,'BlockFading') && ~strcmp(LTE_params.ChanMod_config.type, 'TR 36.873')
    warning('UE speed has no effect for independent block fading channel realizations (LTE_params.ChanMod_config.time_correlation=''independent'')');
end
if (LTE_params.UE_config.user_speed < eps) && (strcmp(LTE_params.ChanMod_config.filtering,'FastFading'))
    error('please use block fading model for zero speed simulations');
end

if strcmp(LTE_params.ChanMod_config.time_correlation, 'correlated') && strcmp(LTE_params.ChanMod_config.type, 'TR 36.873')
    warning('Time correlation setting has no effect when 3D channel model is used.');
end

% check feedback parameters
if (LTE_params.UE_config.RI_fb) && (~LTE_params.UE_config.PMI_fb) && (LTE_params.UE_config.PMI~=0)
    error('for simulation with RI_fb but without PMI_fb, the PMI hast to be zero');
end
if (LTE_params.UE_config.mode == 4) && ( (~LTE_params.UE_config.RI_fb) && (~LTE_params.UE_config.PMI_fb) )
    warning('make sure that LTE_params.UE_config.PMI is a valid codebook entry');
end

if strcmp(LTE_params.UE_config.MCS_and_scheduling_CSI,'estimated') && (LTE_params.nUE>1)
    error('the estimated channel can only be used for link adaption with a single user for now');
end

if strcmp(LTE_params.UE_config.MCS_and_scheduling_CSI,'perfect') && (LTE_params.BS_config.channel_prediction)
    warning('impact of downlink feedback channel delay will not be visible with perfect CSI feedback information');
end

if strcmp(LTE_params.UE_config.MCS_and_scheduling_CSI,'estimated') && strcmp(LTE_params.BS_config.channel_estimation_method,'PERFECT')
    warning('please note: in case of LTE_params.BS_config.channel_estimation_method=''PERFECT'', the CSI for link adaptation is perfect as well');
end
if (~LTE_params.UE_config.RI_fb) && (LTE_params.UE_config.RI > LTE_params.UE_config.nTX)
    error('rank indicator (RI) must not be greater then number of transmit antennas');
end

% scheduling
if ~(strcmp(LTE_params.scheduler.type, 'round robin') || strcmp(LTE_params.scheduler.type, 'fixed')) && LTE_params.UE_config.mode == 4
    error('currently CLSM (mode 4) is not supported by the selected scheduler');
end
if ~(strcmp(LTE_params.scheduler.type, 'round robin') || strcmp(LTE_params.scheduler.type, 'fixed')) && LTE_params.UE_config.CQI_fb==false
    warning('the selected scheduling method has no effect without link adaptation');
end

if (LTE_params.UE_config.mode ~= 5) && (strcmp(LTE_params.scheduler.type, 'fixed MU MIMO') || strcmp(LTE_params.scheduler.type, 'max MU MIMO') || ...
       strcmp(LTE_params.scheduler.type, 'random MU MIMO') ||  strcmp(LTE_params.scheduler.type, 'greedy MU MIMO') || strcmp(LTE_params.scheduler.type, 'greedy frobenius MU MIMO'))
    error('MU MIMO schedulers only work together with transmission mode 5');
end

% maybe introduce a MU MIMO property?
if (LTE_params.UE_config.mode == 5) && ~(strcmp(LTE_params.scheduler.type, 'fixed MU MIMO') || strcmp(LTE_params.scheduler.type, 'max MU MIMO') ...
        || strcmp(LTE_params.scheduler.type, 'random MU MIMO') || strcmp(LTE_params.scheduler.type, 'greedy MU MIMO') || strcmp(LTE_params.scheduler.type, 'greedy frobenius MU MIMO'))
   error('currently only the fixed -, max-, random-, greedy- and greedy frobenius MU MIMO scheduler are supporting transmission mode 5');
end

% no MUMIMO with IAMMSE for now
if (LTE_params.UE_config.mode == 5) &&  strcmp(LTE_params.BS_config.receiver, 'IAMMSE')
   error('currently IAMMSE is not supported with MU MIMO');
end

% IAMMSE
% only perfect channel knowledge
if ( ~strcmp(LTE_params.BS_config.channel_estimation_method, 'PERFECT') && strcmp(LTE_params.BS_config.receiver, 'IAMMSE') )
    error('only PERFECT channel knowledge is supported while using IAMMSE receiver');
end

% only mode 1
if ( strcmp(LTE_params.BS_config.receiver, 'IAMMSE') && LTE_params.UE_config.mode ~= 1 )
    error('only transmission mode 1 is supported with IAMMSE');
end


% transmission mode
switch LTE_params.UE_config.mode
    case 1  % MODE 1 - single transmit antenna
        if LTE_params.UE_config.nTX > 1
            error('Mode 1 only working with 1 UE transmit antenna');
        end
    case 4  % MODE 4 - CLSM
        if LTE_params.UE_config.nTX == 1
            error('Closed loop spatial multiplexing not supported for one transmit antenna');
        end
    case 5
        if LTE_params.UE_config.nTX > 1
            error('Mode 5 only working with 1 UE transmit antenna');
        end
    otherwise
        error('transmission mode invalid');
end 
