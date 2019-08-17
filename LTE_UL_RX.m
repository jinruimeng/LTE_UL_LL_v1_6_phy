function H_test_results=LTE_UL_RX(LTE_params, ChanMod_output, ChanMod, SNR, subframe_i, BS, UE, UE_output, BS_output, y_alice, UE_MCS_and_scheduling_info, bb)
% LTE uplink receiver for a specific user.
% [UE_output, ACK, UE] = LTE_RX(chan_output, SNR, AtPort, subframe_i, BS, UE, BS_output)
% based on LTE_RX, by Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% Authors Jan Propopec, Jiri Blumenstein, prokopec@feec.vutbr.cz
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
% www.urel.feec.vutbr.cz
%
% input :   chan_output         ... [1 x 1]struct - channel output:
%                                   [LTE_params.TxSymbols x 1]double receive signal y_rx
%           SNR                 ... [1 x 1]double current SNR value
%           AtPort              ... [1 x 1]double antenna port number
%           subframe_i          ... [1 x 1]double number of the subframe transmitted
%           BS                  ... [1 x nBS]struct - Base stations parameters from LTE_load_parameters.m
%           UE                  ... [1 x nUE]struct - User equipments capabilities from LTE_load_parameters.m
%           BS_output           ... [1 x nBS]struct - Base Station output
%                                   [LTE_params.TxSymbols x 1]double transmit signal y_tx,
%                                   [1 x nUE]struct - scheduler output and coding parameters, e.g:
%                                   [1 x # data bits]logical genie.data_bits
%                                   [1 x # sent bits(coded)]logical
%                                   genie.sent_bits
%                                   [1 x nUE]struct - user specific parameters and variables
% output:   UE_output           ... [1 x nUE]struct UE output for rx_coded_bits bits before decoder and after demapper
%                                                             and logical decoded bits
%           ACK                 ... [N_subframes x length(SNR_vec)]logical positive or negative acknowledgement
%           UE                  ... [1 x 1]struct - Updated UE struct


global DEBUG_LEVEL
global global_ISI

all_UE_numbers = 1:(LTE_params.nUE*LTE_params.nBS);
global_UE_numbers = all_UE_numbers(logical(LTE_params.connection_table(bb,:)));

for ul = 1:LTE_params.nUE

    uu = global_UE_numbers(ul); %% uu is the global number (1... nUE*nBS), while ul is the local number (1..nUE)
    UE_output(1,uu).UE_scheduled = 0;
    
    if (utils.isScheduled(UE_MCS_and_scheduling_info, bb, uu, LTE_params.connection_table)) 
        
        UE_signaling_UL = BS_output.UE_signaling_UL(bb,ul);
        genie_data = UE_output(1,uu).UE_genie;

        sigma_n2 = 10^(-SNR/10);
        % Perform reception only if there is data to receive
        if(UE_signaling_UL.MCS_and_scheduling_UL.assigned_RBs)

            % Get necessary data and convert in in the correct format
            nLayers    = UE_signaling_UL.MCS_and_scheduling_UL.nLayers;
            nCodewords = UE_signaling_UL.MCS_and_scheduling_UL.nCodewords;
            %UE_mapping = BS_output.UE_signaling_UL(1,uu).MCS_and_scheduling_UL.UE_mapping;

            %% Calculate the correct number of subframe from 1 - 10
            subframe_corr = mod(subframe_i,10);
            if(subframe_corr == 0)
                subframe_corr = 10;
            end

            % must override downlink mapping
            CHmapping = UE_output(1,uu).UE_genie.CH_mapping;

            y_rx_sync = y_alice{bb}.input;

            %% Remove CP, FFT, remove zeros
            y_rx_assembled = zeros(LTE_params.Nsc*LTE_params.Nrb,LTE_params.Nsub,BS.nRX);
            for rr = 1:BS.nRX
                y_rx_resolved = cell(4,1);

                y_rx_resolved{1} = y_rx_sync(1:LTE_params.NfftCP{1},rr);
                y_rx_resolved{2} = reshape( y_rx_sync(LTE_params.NfftCP{1}+1:LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),rr), LTE_params.NfftCP{2}, LTE_params.Ns-1);
                y_rx_resolved{3} = y_rx_sync(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),rr);
                y_rx_resolved{4} = reshape( y_rx_sync(2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:end,rr), LTE_params.NfftCP{2}, LTE_params.Ns-1);

                y_rx_assembled_ifft = [y_rx_resolved{1}(LTE_params.Index_RxCyclicPrefix{1},:)...
                                       y_rx_resolved{2}(LTE_params.Index_RxCyclicPrefix{2},:)...
                                       y_rx_resolved{3}(LTE_params.Index_RxCyclicPrefix{1},:)...
                                       y_rx_resolved{4}(LTE_params.Index_RxCyclicPrefix{2},:)];

                % FFT
                y_rx_assembled_padded = fft(y_rx_assembled_ifft)/sqrt(LTE_params.Nfft);

                % shift from baseband
                y_rx_assembled_shifted = circshift(y_rx_assembled_padded,(LTE_params.Ntot)/2);

                % remove zero padding
                y_rx_assembled(:,:,rr) = y_rx_assembled_shifted(1:LTE_params.Ntot,:);         
            end    

            %% Disassemble reference symbols
            RefMapping = genie_data.DM_mapping; % user allocation of RE
            RefSym = genie_data.rs;
            rx_ref_symbols = y_rx_assembled(repmat(RefMapping,[1,1,LTE_params.BS_config.nRX]));
            rx_ref_symbols = reshape(rx_ref_symbols,size(RefSym,1),[],LTE_params.BS_config.nRX);

            %% Disassemble symbols
            CHmapping = repmat(CHmapping,[1,1,ChanMod.nRX]);
            rx_user_symbols = y_rx_assembled(CHmapping);
            rx_user_symbols = reshape(rx_user_symbols,size(RefSym,1),[],ChanMod.nRX);

            %% Channel and noise estimation
            W = UE_output(1,uu).LayerMapping;                  % precoding matrix W
            CHmapping = repmat(CHmapping,[1,1,1,ChanMod.nTX]);
            SCmapping = repmat(CHmapping(:,1,1,1),[1 size(CHmapping,2) ChanMod.nRX ChanMod.nTX]);

            perfectChannel = ChanMod_output{bb,uu}.genie.H_fft(SCmapping);
            perfectChannel = reshape(perfectChannel, [size(RefSym,1) size(SCmapping,2) ChanMod.nRX ChanMod.nTX]);

            % multiply perfect channel with precoding matrix to get the effective
            % channel
            if max(size(W)) > 1                 % always fullfilled in transmission mode 4
                sizes = size(perfectChannel);
                sizes(end) = size(W,2);         % size(W,2) is the number or layers

                perfectChannel = reshape(perfectChannel,[],size(W,1));
                perfectChannel = perfectChannel*W;
                perfectChannel = reshape(perfectChannel,sizes);
            else
                perfectChannel = perfectChannel*W;
            end

            [H_est,H_pred,H_test] = LTE_UL_channel_estimator(LTE_params,ChanMod,BS,rx_ref_symbols,RefSym,RefMapping,perfectChannel,sigma_n2,genie_data,UE_signaling_UL.MCS_and_scheduling_UL.nLayers,subframe_i,SNR,UE(uu).channel_autocorrelation_matrix,UE_output(1,uu).LayerMapping,UE(uu));  

            % save channel for downlink delay
            UE_output(1,uu).UE_genie.channel_estimate_complete  = H_pred;
            UE_output(1,uu).UE_genie.perfect_channel            = perfectChannel;
            UE_output(1,uu).channel_estimation_error            = mean(mean( sum(sum((abs(perfectChannel - H_est).^2),1),2)/(size(perfectChannel,1)*size(perfectChannel,2)) ));
          

            %% Perform detection
            [LLR_SD,M,H_back] = LTE_detecting(UE_signaling_UL.MCS_and_scheduling_UL,BS.nAtPort,rx_user_symbols,BS.nRX,H_est,ChanMod.filtering,LTE_params,BS.receiver,sigma_n2,UE_output,uu, bb);


            %% Undo layer mapping
            % So we have codewords again
            switch nLayers
                case 1
                    LLR_SS{1} = reshape(LLR_SD(M(1):-1:1,:),1,[]).';
                case 2
                    LLR_SS{1} = reshape(LLR_SD(M(1):-1:1,:),1,[]).';
                    LLR_SS{2} = reshape(LLR_SD(M(1)+M(2):-1:M(1)+1,:),1,[]).';
                case 3
                    LLR_SS{1} = reshape(LLR_SD(M(1):-1:1,:),1,[]).';
                    LLR_SS{2} = reshape(LLR_SD([M(1)+M(2):-1:M(1)+1 M(1)+M(2)+M(3):-1:M(1)+M(2)+1],:),1,[]).';
                case 4
                    LLR_SS{1} = reshape(LLR_SD([M(1):-1:+1 M(1)+M(2):-1:M(1)+1],:),1,[]).';
                    LLR_SS{2} = reshape(LLR_SD([M(3)+M(2)+M(1):-1:M(2)+M(1)+1 M(4)+M(3)+M(2)+M(1):-1:M(3)+M(2)+M(1)+1],:),1,[]).';
            end

            rx_scrambled_bits = cell(size(LLR_SS));
            LLR_SS_descrambled = cell(size(LLR_SS));

            %% Decoding
            for i = 1:UE_signaling_UL.MCS_and_scheduling_UL.nCodewords
                %%%%%%%%%%%%% Soft Sphere Decoder Solution
                rx_scrambled_bits{i} = (1+sign(LLR_SS{i}.'))/2; % NOTE: I don't get what this line is for. Is to get the uncoded BER performance? => yup, exactly, and uncoded throughput and so on

                %% Decrambling of the bits, as of TS 36.211 V11.4.0 (2013-09) Section 5.3.1
                % for rx_coded_bits there is mode 'scramble', because descrambler operates on a bit basis that is equal to scrambling operation
                UE_output(1,uu).rx_coded_bits{i} = LTE_common_scrambling(rx_scrambled_bits{i},BS.NIDcell,uu,subframe_corr,i,'scramble');    
                % for rx_data_bits the LLRs have to be descrambled first
                LLR_SS_descrambled{i} = LTE_common_scrambling(LLR_SS{i},BS.NIDcell,uu,subframe_corr,i,'descramble');
                %% PUSCH: Decoding of the bits
                [UE_output(1,uu).rx_data_bits{i}, ~, UE_output(1,uu).ACK(i)] = LTE_UL_rx_ULSCH_decode(LTE_params,LLR_SS_descrambled{i},UE_signaling_UL,BS,UE(uu),BS.UE_specific(1,ul),genie_data,i);
            end

        %     % Feedback performed only when data is received
            feedback_rv_idx = [UE_signaling_UL.turbo_rate_matcher_UL.rv_idx];
            UE_output(1,uu).rv_idx = feedback_rv_idx;
        %     
        %     % Precoding feedback calculation
        %     if LTE_params.downlink_delay ~= 0
        %         error('not yet supported') % spratsch; feedback only implemented for 0 delay
        %         [RI_tmp,PMI_tmp,CQI_tmp]=LTE_UL_feedback(BS.nAtPort,sigma_n2,LTE_params,H_est,UE(uu),uu,LTE_params.UE_config.mode);
        %         UE_output(1,uu).RI = RI_tmp;
        %         UE_output(1,uu).PMI = PMI_tmp-1;
        %         UE_output(1,uu).CQI = CQI_tmp;
        %     end
        %     
        %     %% check the deadline of the data
        %     UE(uu).traffic_model.check_deadline(BS.clock.current_TTI);
        else
        %     % Feedback for the cases where no data was transmitted
            UE_output(1,uu).rv_idx = zeros(1,2);
            UE_output(1,uu).ACK(:) = 0;
        %     H_est = ChanMod_output.genie.H_fft;
        %     % Precoding feedback calculation. Use QPSK modulation for calculations when no data was assigned (assume the channel was really bad)
        %     if LTE_params.downlink_delay ~= 0
        %         error('not yet supported')
        %         [RI_tmp,PMI_tmp,CQI_tmp]=LTE_UL_feedback(BS.nAtPort,sigma_n2,LTE_params,H_est,UE(uu),uu,LTE_params.UE_config.mode);
        %         UE_output(1,uu).RI = RI_tmp;
        %         UE_output(1,uu).PMI = PMI_tmp-1;
        %         UE_output(1,uu).CQI = CQI_tmp;
        %     end
        %     
        %     %% check the deadline of the data
        %     UE(uu).traffic_model.check_deadline(BS.clock.current_TTI);
        end


        % Finish feedback construction
        % For now use perfect channel information. Can be changed to an SINR computed from channel estimation
        % SINR = BS_output.cell_genie.SINR(uu,:);
        % SINRs_to_CQI = SINR(1:12:end);
        % UE_output(uu).CQI_feedback = floor(LTE_common_CQI_mapping(LTE_params.CQI_mapping,SINRs_to_CQI)); % only one CQI value over bandwidth, because there is only noise added
        % UE_output(uu).CQI = UE_output(uu).CQI_feedback;
        % UE_output(uu).PMI = 1;
        % HARQ process id of this current ACK
        UE_output(1,uu).HARQ_process = UE_signaling_UL.MCS_and_scheduling_UL.HARQ_process_id;
        UE_output(1,uu).UE_scheduled = (UE_signaling_UL.MCS_and_scheduling_UL.assigned_RBs>0);
    end
     H_test_results=H_test;
end
 
