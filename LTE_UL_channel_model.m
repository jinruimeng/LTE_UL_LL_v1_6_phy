function [ChanMod_output, BS_input] = LTE_UL_channel_model(LTE_params,ChanMod,ChanMod_output, UE_output, SNR, UE_MCS_and_scheduling_info)
% LTE channel model - to filter the output of the transmitter.
% [chan_output] = LTE_channel_model(UE_output, SNR)
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input :   UE_output           ... [1x1]struct - UE output
%                                   [LTE_params.TxSymbols x 1]double transmit signal y_tx,
%                                   [1x1]struct - scheduler output and coding parameters, e.g:
%                                   [1 x # data bits]logical genie.data_bits
%                                   [1 x # sent bits(coded)]logical genie.sent_bits
%           SNR                 ... [1]double current SNR value
% output:   chan_output         ... [nBS x (nUE * nBS)] struct - channel output:
%                                   [LTE_params.TxSymbols x 1]double receive signal y_rx

%% Channel Model
% y_tx_filtered = zeros(size(UE_output(1).y_tx,1) + channel_size - 1,ChanMod.nRX);

for bb = 1:LTE_params.nBS
    for uu = 1:LTE_params.nUE*LTE_params.nBS
        if(utils.isScheduled(UE_MCS_and_scheduling_info, 1:LTE_params.nBS, uu, LTE_params.connection_table))
            switch ChanMod.type
                case {'AWGN','flat Rayleigh'}
                    y_tx_filtered = UE_output(uu).y_tx*transpose(ChanMod_output{bb,uu}.H);
                case {'Rayleigh', 'TR 36.873'}
                    H_full_fft = permute(ChanMod_output{bb,uu}.genie.H_fft,[3,4,1,2]);
                    y_tx_fft = UE_output(uu).UE_genie.y_tx_assembled;
                    c1 = size(ChanMod_output{bb,uu}.genie.H_fft,1);
                    c2 = size(ChanMod_output{bb,uu}.genie.H_fft,2);
                    y_filter_fft = zeros(c1,c2,ChanMod(1).nRX);
                    for cc1 = 1:c1
                        for cc2 = 1:c2
                            y_filter_fft(cc1,cc2,:) =  H_full_fft(:,:,cc1,cc2)*y_tx_fft(:,cc1,cc2);
                        end
                    end
                    Nsub = LTE_params.Nsub;
                    Ns = LTE_params.Ns;
                    for nn = 1:ChanMod(1).nRX
                        y_filter_padded = [y_filter_fft(:,:,nn); zeros(LTE_params.Nfft - LTE_params.Ntot,Nsub)];
                        y_filter_shifted = circshift(y_filter_padded,-LTE_params.Ntot/2);
                        y_filter_ifft = sqrt(LTE_params.Nfft)*ifft(y_filter_shifted);
                        y_tx_filtered(:,nn) = [    y_filter_ifft(LTE_params.Index_TxCyclicPrefix{1},1);...
                            reshape(y_filter_ifft(LTE_params.Index_TxCyclicPrefix{2},2:Ns),[],1);...
                            y_filter_ifft(LTE_params.Index_TxCyclicPrefix{1},Ns+1);...
                            reshape(y_filter_ifft(LTE_params.Index_TxCyclicPrefix{2},Ns+2:end),[],1) ];
                    end
                    %y_tx_filtered = UE_output(uu).y_tx*transpose(ChanMod_output{bb,uu}.H); % spratsch: does this work for multiple BS ?
                    %       y_tx_filtered = y_tx_filtered.';
                case {'PedA','PedB','PedBcorr','VehA','VehB','TU','RA','HT','EPedA','EVehA','ETU','winner_II','ePDP'}
                    switch ChanMod.filtering
                        case 'BlockFading'
                            channel_size = size(ChanMod_output{bb,uu}.H,3);
                        case 'FastFading'
                            channel_size = size(ChanMod_output{bb,uu}.H,4);
                    end
                    y_tx_filtered = zeros(size(UE_output(uu).y_tx,1) + channel_size - 1,ChanMod.nRX);
                    for rr = 1:ChanMod.nRX
                        for tt = 1:ChanMod.nTX
                            switch ChanMod.filtering
                                case 'BlockFading'
                                    y_tx_filtered(:,rr) = y_tx_filtered(:,rr) + conv(squeeze(ChanMod_output{bb,uu}.H(rr,tt,:)), UE_output(uu).y_tx(:,tt));
                                case 'FastFading'
                                    taps_num = size(squeeze(ChanMod_output{bb,uu}.H(rr,tt,:,:)),2);    %number of channel taps
                                    y_tx_filtered_matrix = [repmat(UE_output(uu).y_tx(:,tt),1,taps_num).*squeeze(ChanMod_output{bb,uu}.H(rr,tt,:,:))].'; %multiplication of the samples with their channel
                                    
                                    if(strcmp(LTE_params.CyclicPrefix,'normal'))
                                        mapping_longer = logical(toeplitz([1;zeros(LTE_params.NfftCP{1}-1,1)],[ones(taps_num,1);zeros(LTE_params.NfftCP{1}-1,1)])); %the structur for time variant convolution for first and 7th symbols
                                        mapping_shorter = mapping_longer(1:LTE_params.NfftCP{2},1:(LTE_params.NfftCP{2} + taps_num - 1));   %mapping for other symbols, same as previous, but we just removed last column and row
                                        y_tx_filtered_part = zeros(length(UE_output(uu).y_tx)+taps_num-1,1);
                                        start = 1;
                                        % time-variant convolution
                                        %loop over ofdm symbols
                                        for symbol_i = 1:LTE_params.Nsub
                                            if(symbol_i == 1 || symbol_i == 7) %the calculation of the output for longer symbols
                                                y_tx_filtered_part_shifted_longer = zeros(LTE_params.NfftCP{1} + taps_num - 1,LTE_params.NfftCP{1});    %reseting of the convolution structure
                                                stop = start + LTE_params.NfftCP{1} + taps_num - 1 - 1;
                                                y_tx_filtered_part_shifted_longer(mapping_longer') = reshape(y_tx_filtered_matrix(:,start:stop - taps_num + 1),[],1);
                                                y_tx_filtered_part(start:stop,1) = y_tx_filtered_part(start:stop,1) + sum(y_tx_filtered_part_shifted_longer,2);
                                                
                                                start = start + LTE_params.NfftCP{1};
                                            else %the calculation of the output for shorter symbols
                                                y_tx_filtered_part_shifted_shorter = zeros(LTE_params.NfftCP{2} + taps_num - 1,LTE_params.NfftCP{2});    %reseting of the convolution structure
                                                stop = start + LTE_params.NfftCP{2} + taps_num - 1 - 1;
                                                y_tx_filtered_part_shifted_shorter(mapping_shorter') = reshape(y_tx_filtered_matrix(:,start:stop - taps_num + 1),[],1);
                                                y_tx_filtered_part(start:stop,1) = y_tx_filtered_part(start:stop,1) + sum(y_tx_filtered_part_shifted_shorter,2);
                                                start = start + LTE_params.NfftCP{2};
                                            end
                                            %
                                        end
                                    else
                                        mapping = logical(toeplitz([1;zeros(LTE_params.NfftCP{1}-1,1)],[ones(taps_num,1);zeros(LTE_params.NfftCP{1}-1,1)]));
                                        y_tx_filtered_part = zeros(length(UE_output(uu).y_tx)+taps_num-1,1);
                                        start = 1;
                                        % time-variant convolution
                                        %loop over ofdm symbols
                                        for symbol_i = 1:LTE_params.Nsub
                                            y_tx_filtered_part_shifted_longer = zeros(LTE_params.NfftCP{1} + taps_num - 1,LTE_params.NfftCP{1});    %reseting of the convolution structure
                                            stop = start + LTE_params.NfftCP{1} + taps_num - 1 - 1;
                                            y_tx_filtered_part_shifted_longer(mapping') = reshape(y_tx_filtered_matrix(:,start:stop - taps_num + 1),[],1);
                                            y_tx_filtered_part(start:stop,1) = y_tx_filtered_part(start:stop,1) + sum(y_tx_filtered_part_shifted_longer,2);
                                            
                                            start = start + LTE_params.NfftCP{1};
                                            %
                                        end
                                    end
                                    %
                                    y_tx_filtered(:,rr) = y_tx_filtered(:,rr) + y_tx_filtered_part;
                            end
                        end
                    end
            end
            %     %% Add Noise
%                  noise_RandStream = LTE_params.noise_RandStream;
%                  n = (randn(noise_RandStream,size(y_tx_filtered)) + 1i*randn(noise_RandStream,size(y_tx_filtered)))*10^(-SNR/20)/sqrt(2);
%                  v = n;
%             
%                  y_tx_filtered = y_tx_filtered + v;
            
            ChanMod_output{bb,uu}.y_rx = y_tx_filtered(1:length(UE_output(uu).y_tx),:);
            
            %% Add genie information
            %     UE_output(uu).UE_genie.v = v;     %noise at RX-antennas (preFFT), the comment must be fixed
            %     UE_output(uu).UE_genie.n = n;     %noise after FFT at receiver (postFFT) the comment must be fixed
        end
    end
end % end of bb loop

%%

[BS_input,n_rstream, signal_size] = LTE_UL_generate_BS_input(LTE_params, ChanMod_output, SNR, UE_MCS_and_scheduling_info);

%%ADD NOISE
BS_input{1,1}.input = BS_input{1,1}.input+( randn(n_rstream, signal_size) + 1i*randn(n_rstream, signal_size) ) * 10^(-SNR/20)/sqrt(2);
BS_input{2,1}.input = BS_input{1,1}.input+( randn(n_rstream, signal_size) + 1i*randn(n_rstream, signal_size) ) * 10^(-SNR/20)/sqrt(2);
end % end of function

