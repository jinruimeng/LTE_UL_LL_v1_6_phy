function [BS_output, UE_output,H_test_results] = LTE_UL_sim_main_subframe(LTE_params, ChanMod, SNR, UE, UE_output, BS, BS_output, subframe_i, out, downlinkChannel, cqi_i, channel)

global sampled_covariance
global samples
global DEBUG_LEVEL

nUE = LTE_params.nUE;
nBS = LTE_params.nBS;

ChanMod_output = cell(nBS, nUE*nBS);

%% channel matrix, for scheduler purposes, zero delay from RX
switch ChanMod.type
    case 'winner_II' % check this ( mbs )
        for bb = 1:LTE_params.nBS 
            for uu = 1:LTE_params.nUE *LTE_params.nBS % parallel channels for multi-user
                % NOTE: every user should have a different MIMO channel matrix
                switch ChanMod.filtering
                    case 'BlockFading'
                        ChanMod_output{bb,uu} = LTE_UL_channel_matrix(LTE_params, ChanMod, UE_output, UE, uu, LTE_params.TxSymbols, subframe_i,channel{uu}(:,:,:,subframe_i),out);
                    case 'FastFading'
                        switch ChanMod.time_correlation
                            case 'correlated'
                                [channel, delays, out] = LTE_winner_channel_model(LTE_params.TxSymbols,LTE_params.Arrays,out);
                            case 'independent'
                                [channel, delays, out] = LTE_winner_channel_model(LTE_params.TxSymbols,LTE_params.Arrays);
                        end      
                        ChanMod_output{bb,uu} = LTE_UL_channel_matrix(LTE_params, ChanMod, UE_output, UE, uu, LTE_params.TxSymbols, subframe_i,channel{uu},out);
                    otherwise
                        error('channel filtering type invalid');
                end
            end
        end
    case 'TR 36.873'    
        for bb = 1:LTE_params.nBS 
            for uu = 1:LTE_params.nUE*LTE_params.nBS % parallel channels for multi-user
                ChanMod_output{bb, uu}.genie.H_fft = channel{uu}(:, :, :, :, subframe_i);
            end
        end
    otherwise
        for bb = 1:LTE_params.nBS
            for uu = 1:LTE_params.nUE*LTE_params.nBS  
                ChanMod_output{bb,uu} = LTE_UL_channel_matrix(LTE_params, ChanMod, UE_output, UE, uu, LTE_params.TxSymbols, subframe_i);
            end
        end
end
% Receive feedback from the previous subframe
UE_output = downlinkChannel.receive_feedback;

% ACK of the previous frame. If this is the first frame, set the
% ACK to correct so that the HARQ handling generates new data
if subframe_i==1 || (subframe_i <= LTE_params.downlink_delay && ~strcmp(LTE_params.ChanMod_config.time_correlation,'independent'))% (number of max HARQ processes that will be used)
    for uu = 1:LTE_params.nUE*LTE_params.nBS
        UE_output(uu).ACK    = true(1,2);
        UE_output(uu).rv_idx = zeros(1,2);
    end
end

if DEBUG_LEVEL > 3
    HH=ChanMod_output{1}.genie.H_fft;
    HH=squeeze(HH(:,1,1,1));
    HH=HH(:);
    
    samples = samples + 1;
    sampled_covariance = sampled_covariance+HH*HH';
    %     sampled_covariance = sampled_covariance+hh*hh';
    
    imagesc(abs(sampled_covariance/samples));
    colorbar;
end

%% Transmitter
% changed 18.11.2011, for each user is signal generated seperately,
% channels are independent, perfect synchronization due to
% alignment after this loop
subframe_corr = mod(subframe_i-1,LTE_params.Nsfr)+1;
srs_subframe = max(mod(subframe_corr-1,LTE_params.SRS_period) == LTE_params.SRS_offset);
if LTE_params.SRS_period == 0 || srs_subframe == 0;
    SyncUsedElements = zeros(LTE_params.Nrb,2);
    CHusedElements = 12 * ones(LTE_params.Nrb*2,1); % one symbol in each slot used for DMRS
    srs_subframe = 0;
else
    SyncUsedElements = zeros(LTE_params.Nrb,2);
    SyncUsedElements(:,2) = 12;                     % SRS in last OFDM symbol in subframe
    CHusedElements = 12 * ones(LTE_params.Nrb*2,1); % one symbol in each slot used for DMRS
end

%% Scheduling and link adaptation
% for each base station bb
for bb = 1:LTE_params.nBS
    % prepare channel for link adaptation
    ind = find(LTE_params.connection_table(bb,:)); 
    if(subframe_i <= LTE_params.downlink_delay) || (LTE_params.downlink_delay==0)       % as long as there is no feedback yet
        tmp = ChanMod_output(bb,logical(LTE_params.connection_table(bb,:).'));
        for ul=1:length(ind)
            channel{ul} = tmp{ul}.genie.H_fft;                                          % take perfect channel
            
            % channel prediction error calculation
            uG = utils.localToGlobalUser(LTE_params.connection_table, bb, ul);
            UE_output(uG).channel_prediction_error = 0;
        end
    else
        for ul=1:length(ind)
            uG = utils.localToGlobalUser(LTE_params.connection_table, bb, ul);
            
            if strcmp(LTE_params.UE_config.MCS_and_scheduling_CSI,'perfect')
                channel{ul} = UE_output(uG).UE_genie.perfect_channel;
            else
                channel{ul} = UE_output(uG).UE_genie.channel_estimate_complete;    % take estimated channel (already delayed)
            end
            
            % channel prediction error calculation
            if LTE_params.BS_config.channel_prediction
                UE_output(uG).channel_prediction_error = mean(mean(mean(mean( (abs(UE_output(uG).UE_genie.perfect_channel - UE_output(uG).UE_genie.channel_estimate_complete).^2) ))));
            end
        end
    end
    
    % prepare channel estimation MSE for link adaptation
    for ul = 1:length(ind)
        uG = utils.localToGlobalUser(LTE_params.connection_table, bb, ul);
        if strcmp(LTE_params.BS_config.channel_estimation_method,'PERFECT') || (subframe_i <= LTE_params.downlink_delay) || (subframe_i==1)
            channel_MSE(uG) = 0;                                                        % as long as there is no estimate, assume zero MSE
        else
            if ~isempty(UE_output(uG).channel_estimation_error)
                UE(uG).MSE_buffer = [ UE(uG).MSE_buffer(2:end), UE_output(uG).channel_estimation_error ];
            end
            if(UE(uG).MSE_buffer>0)                                                     % if there is at least one nonzero MSE
                channel_MSE(uG) = mean(UE(uG).MSE_buffer(UE(uG).MSE_buffer~=0));        % average MSE from previous subframes
            else
                channel_MSE(uG) = 0;                                                    % set MSE to zero if buffer is all zero
            end
        end
    end
       
    % actual scheduling and link adaptation
    UE_MCS_and_scheduling_info(bb,:) = BS(bb).scheduler.scheduler_users(subframe_corr,SyncUsedElements,BS(bb),CHusedElements,cqi_i, channel, 10.^(-SNR/10), channel_MSE);
end

%% assign a cyclic shift field in DCI format to each UE
% get number of scheduled users per base station
% n_schedUE = zeros(nBS,1);
% for uu=1:nUE*nBS
%     bb = utils.findBS(LTE_params.connection_table, uu);
%     local_uu = utils.globalToLocalUser(LTE_params.connection_table, bb, uu);
%     
%     if UE_MCS_and_scheduling_info(bb,local_uu).assigned_RBs
%         n_schedUE(bb) = n_schedUE(bb) + 1;
%     end
% end

% assign cyclic shift field to users

scheduled_users_per_BS = cell(nBS,1);
for bb = 1:nBS
    scheduled_users_global = [];
    for ul = 1:nUE
        uG = utils.localToGlobalUser(LTE_params.connection_table, bb, ul);
        if UE_MCS_and_scheduling_info(bb,ul).assigned_RBs > 1
            scheduled_users_global = [scheduled_users_global, uG]; %#ok
        end
    end
    scheduled_users_per_BS{bb} = scheduled_users_global;
end

for bb = 1:nBS
    if length(scheduled_users_per_BS{bb}) <= 4
        for ul_scheduled = 1:length(scheduled_users_per_BS{bb})
            uG = scheduled_users_per_BS{bb}(ul_scheduled);
            UE(uG).CSFiDCIFormat = LTE_params.CSF_mapping{1}(ul_scheduled);
        end
    else
        for ul_scheduled = 1:length(scheduled_users_per_BS{bb})
            uG = scheduled_users_per_BS{bb}(ul_scheduled);
            UE(uG).CSFiDCIFormat = LTE_params.CSF_mapping{2}(ul_scheduled);
        end
    end
        
end

% for uu=1:nUE*nBS
%     bb = utils.findBS(LTE_params.connection_table, uu);
%     local_uu = utils.globalToLocalUser(LTE_params.connection_table, bb, uu);
%     
%     if(n_schedUE(bb) <= 4)
%         if UE_MCS_and_scheduling_info(bb,local_uu).assigned_RBs
%             UE(uu).CSFiDCIFormat = LTE_params.CSF_mapping{1}(local_uu);
%         end
%     else
%         if UE_MCS_and_scheduling_info(bb,local_uu).assigned_RBs
%             UE(uu).CSFiDCIFormat = LTE_params.CSF_mapping{2}(local_uu);
%         end
%     end
% end


%% Transmission for each user
for uu=1:nUE*nBS
    bb = utils.findBS(LTE_params.connection_table, uu);
    local_uu = utils.globalToLocalUser(LTE_params.connection_table, bb, uu);
    
    % RS allocation
    LTE_params.DM_SRS_allocation = LTE_UL_allocate_RS(LTE_params,UE_MCS_and_scheduling_info(bb,local_uu), UE(uu), subframe_i, subframe_corr, srs_subframe);
    
    % update the UE signaling
    BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL = UE_MCS_and_scheduling_info(bb,local_uu);
    
    % actual transmission for every scheduled UE
    if (utils.isScheduled(UE_MCS_and_scheduling_info, bb, uu, LTE_params.connection_table))
        % generate transmit signal
        LTE_UL_TX(LTE_params, UE(1,uu), BS, subframe_i, UE_output(1,uu), BS_output, uu, bb, srs_subframe); % ref and sync is repeating on slot basis
    else
%         UE(uu).traffic_model.decrease_data_buffer(0,0,0,UE_output(1,uu));
    end
end


%% channel model filtering
[ChanMod_output, BS_input] = LTE_UL_channel_model(LTE_params, ChanMod, ChanMod_output, UE_output, SNR(1), UE_MCS_and_scheduling_info);

y_alice=BS_input{1,1}.input;
y_bob=BS_input{2,1}.input;
%% Receiver
%%
for bb = 1:LTE_params.nBS
    if (LTE_params.UE_config.mode == 5) % MUMIMO
        LTE_UL_RX_MUMIMO(LTE_params, ChanMod_output, ChanMod, SNR, subframe_i, BS(bb), UE, UE_output, BS_output, y_alice, UE_MCS_and_scheduling_info, bb, length(scheduled_users_per_BS{bb}));
    else
        if strcmp(LTE_params.BS_config.receiver, 'IAMMSE')
            LTE_UL_RX_ICI(LTE_params, ChanMod_output, ChanMod, SNR, subframe_i, BS, UE, UE_output, BS_output, y_alice, UE_MCS_and_scheduling_info, bb);
        else
            LTE_UL_RX_alice(LTE_params, ChanMod_output, ChanMod, SNR, subframe_i, BS(bb), UE, UE_output, BS_output, y_alice, UE_MCS_and_scheduling_info, bb);
            H_test_results{1,1}= LTE_UL_RX_alice(LTE_params, ChanMod_output, ChanMod, SNR, subframe_i, BS(bb), UE, UE_output, BS_output, y_alice, UE_MCS_and_scheduling_info, bb);
            
        end
    end
end

%%
%% Receiver
%%
for bb = 1:LTE_params.nBS
    if (LTE_params.UE_config.mode == 5) % MUMIMO
        LTE_UL_RX_MUMIMO(LTE_params, ChanMod_output, ChanMod, SNR, subframe_i, BS(bb), UE, UE_output, BS_output, BS_input, UE_MCS_and_scheduling_info, bb, length(scheduled_users_per_BS{bb}));
    else
        if strcmp(LTE_params.BS_config.receiver, 'IAMMSE')
            LTE_UL_RX_ICI(LTE_params, ChanMod_output, ChanMod, SNR, subframe_i, BS, UE, UE_output, BS_output, BS_input, UE_MCS_and_scheduling_info, bb);
        else
           LTE_UL_RX_bob(LTE_params, ChanMod_output, ChanMod, SNR, subframe_i, BS(bb), UE, UE_output, BS_output, y_bob, UE_MCS_and_scheduling_info, bb);
            
            H_test_results{1,2}=  LTE_UL_RX_bob(LTE_params, ChanMod_output, ChanMod, SNR, subframe_i, BS(bb), UE, UE_output, BS_output, y_bob, UE_MCS_and_scheduling_info, bb);
        end
    end
end
%% Feedback through downlink channel
% store feedback for next subframe 
downlinkChannel.insert_feedback(UE_output);

end

