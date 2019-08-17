function [y_rx_sync, freq_offset_est] = LTE_rx_freq_sync(y_rx, UE, freq_offset_est, RefSym, RefMapping, PrimSync, PrimMapping, SecSync, SecMapping)
% LTE carrier frequency offset compensation
% [y_rx_sync, freq_offset_est] = LTE_rx_freq_sync(y_rx, UE, RefSym, RefMapping)
% Author: Qi Wang, qwang@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input :   y_rx                        ... [ XX x nRX ] received sequence
%           UE                          ... UE parameter
%           freq_offset_est             ... struct
%           RefSym                      ... Reference Symbol
%           RefMapping                  ... Reference Symbol Mapping
%           PrimSync                    ... Primary Synchronization signals
%           PrimMapping                 ... Primary Synchronization signals mapping
%           SecSync                     ... Secondary Synchronization signals
%           SecMapping                  ... Secondary Synchronization signals mapping
% output:   y_rx_sync                   ... [ XX x nRX ] sequence with carrier frequency offset compensated
%           freq_offset_est             ... struct
%
% date of creation: 2009/04/24
% last changes: 2009/06/23 Wang

global LTE_params;

if UE.perfect_freq_sync
    y_rx_sync = y_rx.*exp(-1i*2*pi*UE.carrier_freq_offset*repmat((0:1:(size(y_rx,1)-1))',1,UE.nRX)/LTE_params.Nfft);
    freq_offset_est.error = 0;
else
    nTX = size(RefSym,3);
    if (exist('PrimSync','var'))&&(exist('PrimMapping','var'))
        %% estimate the fractional carrier frequency offset
        if (length(LTE_params.Ng) == 2)
            Index_extractCP_temp{1} = 1:LTE_params.Ng(1);
            Index_extractCP_temp{2} = reshape(LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-1)),LTE_params.NfftCP{2},(LTE_params.Ns-1));
            Index_extractCP_temp{3} = LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+(1:LTE_params.Ng(1));
            Index_extractCP_temp{4} = reshape(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-1)),LTE_params.NfftCP{2},(LTE_params.Ns-1));
            Index_extractCP = [Index_extractCP_temp{1} reshape(Index_extractCP_temp{2}(1:LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3} reshape(Index_extractCP_temp{4}(1:LTE_params.Ng(2),:),1,[])];
        else
            Index_extractCP_temp = reshape((1:LTE_params.NfftCP*LTE_params.Nsub),LTE_params.NfftCP,LTE_params.Nsub);
            Index_extractCP = reshape(Index_extractCP_temp(1:LTE_params.Ng,:),1,[]);
        end
        rx_CP1 = y_rx(Index_extractCP,:);
        rx_CP2 = y_rx(Index_extractCP+LTE_params.Nfft,:);
        expp = sum(sum(rx_CP1 .* conj(rx_CP2)));
        freq_offset_est.frac = -angle(expp)/2/pi;
        y_rx_sync1 = y_rx.*exp(-1i*2*pi*freq_offset_est.frac*repmat((0:1:(size(y_rx,1)-1))',1,UE.nRX)/LTE_params.Nfft);
        %% estimate the integer carrier frequency offset
        for nn = 1:UE.nRX
            y_rx_resolved = cell(4,1);
            if(length(LTE_params.Ng)==2)
                y_rx_resolved{1} = y_rx_sync1(1:LTE_params.NfftCP{1},nn);
                y_rx_resolved{2} = reshape(y_rx_sync1(LTE_params.NfftCP{1}+1:LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),nn),LTE_params.NfftCP{2},LTE_params.Ns-1);
                y_rx_resolved{3} = y_rx_sync1(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),nn);
                y_rx_resolved{4} = reshape(y_rx_sync1(2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:end,nn),LTE_params.NfftCP{2},LTE_params.Ns-1);
                y_rx_assembled_ifft = [y_rx_resolved{1}(LTE_params.Index_RxCyclicPrefix{1},:) y_rx_resolved{2}(LTE_params.Index_RxCyclicPrefix{2},:)...
                    y_rx_resolved{3}(LTE_params.Index_RxCyclicPrefix{1},:) y_rx_resolved{4}(LTE_params.Index_RxCyclicPrefix{2},:)];
            else
                y_rx = reshape(y_rx_sync1(:,nn),LTE_params.NfftCP,LTE_params.Nsub);
                y_rx_assembled_ifft = y_rx(LTE_params.Index_RxCyclicPrefix,:);
            end
            y_rx_assembled_shifted = fft(1/sqrt(LTE_params.Nfft)/(sqrt(LTE_params.Nfft/LTE_params.Ntot))*y_rx_assembled_ifft);
            y_rx_assembled_padded = circshift(y_rx_assembled_shifted,LTE_params.Ntot/2);
            % remove zero DC carrier
            y_rx_assembled_temp = y_rx_assembled_padded([1:LTE_params.Ntot/2,LTE_params.Ntot/2+2:LTE_params.Ntot+1],:);
            PrimSync_rx(:,nn) = y_rx_assembled_temp(PrimMapping);
            SecSync_rx(:,nn) = y_rx_assembled_temp(SecMapping);
            y_rx_assembled(:,:,nn) = y_rx_assembled_temp;
        end

        % ML estimator
        for n_shift = -31:1:31
            for rr = 1:UE.nRX
                obs(rr,n_shift+32) = sum(circshift(PrimSync_rx(:,rr).*conj(SecSync_rx(:,rr)),n_shift).*conj(PrimSync).*SecSync);
            end
        end
        Tp = sum(obs,1).*exp(1i*2*pi*(-31:31)*(LTE_params.NfftCP{2}/LTE_params.Nfft-1));
        [peak, argument] = max(Tp);
        freq_offset_est.int = 32 - argument;

        y_rx_sync = y_rx_sync1.*exp(-1i*2*pi*freq_offset_est.int*repmat((0:1:(size(y_rx,1)-1))',1,UE.nRX)/LTE_params.Nfft);
    else
        y_rx_sync = y_rx.*exp(-1i*2*pi*(freq_offset_est.frac+freq_offset_est.int)*repmat((0:1:(size(y_rx,1)-1))',1,UE.nRX)/LTE_params.Nfft);
    end
    %% estimate the residual carrier frequency offset
    % to frequency domain
    for nn = 1:UE.nRX
        y_rx_resolved = cell(4,1);
        if(length(LTE_params.Ng)==2)
            y_rx_resolved{1} = y_rx_sync(1:LTE_params.NfftCP{1},nn);
            y_rx_resolved{2} = reshape(y_rx_sync(LTE_params.NfftCP{1}+1:LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),nn),LTE_params.NfftCP{2},LTE_params.Ns-1);
            y_rx_resolved{3} = y_rx_sync(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),nn);
            y_rx_resolved{4} = reshape(y_rx_sync(2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:end,nn),LTE_params.NfftCP{2},LTE_params.Ns-1);
            y_rx_assembled_ifft = [y_rx_resolved{1}(LTE_params.Index_RxCyclicPrefix{1},:) y_rx_resolved{2}(LTE_params.Index_RxCyclicPrefix{2},:)...
                y_rx_resolved{3}(LTE_params.Index_RxCyclicPrefix{1},:) y_rx_resolved{4}(LTE_params.Index_RxCyclicPrefix{2},:)];
        else
            y_rx_block = reshape(y_rx_sync(:,nn),LTE_params.NfftCP,LTE_params.Nsub);
            y_rx_assembled_ifft = y_rx_block(LTE_params.Index_RxCyclicPrefix,:);
        end
        y_rx_assembled_shifted = fft(1/sqrt(LTE_params.Nfft)/(sqrt(LTE_params.Nfft/LTE_params.Ntot))*y_rx_assembled_ifft);
        y_rx_assembled_padded = circshift(y_rx_assembled_shifted,LTE_params.Ntot/2);
        % remove zero DC carrier
        y_rx_assembled(:,:,nn) = y_rx_assembled_padded([1:LTE_params.Ntot/2,LTE_params.Ntot/2+2:LTE_params.Ntot+1],:);
    end
        % at most reference symbols from 2 TX are used, even if there can be 4 TX
    if nTX~= 4
        nRefUse = nTX;
    else
        nRefUse = 2;
    end
    RefSym = RefSym(:,:,1:nRefUse,:); % extract the reference symbols that are used for RFO estimation
    RefSym_rx = nan([size(RefSym,1), size(RefSym,2), nRefUse, UE.nRX]);   %allocate memory for reference signal for every channel, the number of transmitt antennas is included in RefSym
    for tt = 1:nRefUse
        for rr = 1:UE.nRX
            RefSym_rx_help = y_rx_assembled(:,:,rr);   %use received symbols from one recieve antenna
            RefSym_rx_help = RefSym_rx_help(RefMapping(:,:,tt));  %extract the signal on pilots positons
            RefSym_rx(:,:,tt,rr) = reshape(RefSym_rx_help,size(RefSym(:,:,tt)));  %finally place the reference signal on allocated position
        end
    end
    switch LTE_params.UE_config.rfo_correct_method
        case 'none'
            freq_offset_est.res = 0;
        case 'subframe'
            for rr = 1:UE.nRX
                freq_offset_est.phi(:,:,:,rr) = RefSym_rx(:,[1,2],:,rr).*conj(RefSym_rx(:,[3,4],:,rr)).*conj(RefSym(:,[1,2],:)).*RefSym(:,[3,4],:);
            end
            U = sum(sum(sum(sum(freq_offset_est.phi))));
            freq_offset_est.res = -angle(U)*LTE_params.Nfft/(length(y_rx_sync)*pi);
        case 'IIR'
            if (isfield(freq_offset_est,'Upre'))
                res_pre = freq_offset_est.sum-freq_offset_est.frac-freq_offset_est.int; % previousCFO - currentFFO -currentIFO 
                freq_offset_est.Upre = exp(-j*2*pi*res_pre*length(y_rx_sync)/2/LTE_params.Nfft/LTE_params.Ns);
            else 
                freq_offset_est.Upre = 0;
            end
            for rr = 1:UE.nRX
                freq_offset_est.phi(:,:,:,rr) = RefSym_rx(:,[1,2],:,rr).*conj(RefSym_rx(:,[3,4],:,rr)).*conj(RefSym(:,[1,2],:)).*RefSym(:,[3,4],:);
            end
            Utemp = sum(sum(sum(sum(freq_offset_est.phi))));
            U = 1/5*Utemp/abs(Utemp)+4/5*freq_offset_est.Upre;
            freq_offset_est.res = -1/2/pi*angle(U)*LTE_params.Nfft*LTE_params.Ns/length(y_rx_sync)/2;
            freq_offset_est.Upre = U/abs(U);
        case '5subframe'
            for rr = 1:UE.nRX
                freq_offset_est.phi(:,:,:,rr) = RefSym_rx(:,[1,2],:,rr).*conj(RefSym_rx(:,[3,4],:,rr)).*conj(RefSym(:,[1,2],:)).*RefSym(:,[3,4],:);
            end
            U = sum(sum(sum(sum(freq_offset_est.phi))));
            if (exist('PrimSync','var'))
                freq_offset_est.phisum = U;
            else
                freq_offset_est.phisum = freq_offset_est.phisum + U;
            end
            freq_offset_est.res = -angle(freq_offset_est.phisum)*LTE_params.Nfft/(length(y_rx_sync)*pi);
            
        otherwise
            error('unknown RFO estimation method.')
    end
    freq_offset_est.sum = freq_offset_est.frac + freq_offset_est.int + freq_offset_est.res;
    freq_offset_est.error = UE.carrier_freq_offset-freq_offset_est.sum;
    y_rx_sync = y_rx_sync.*exp(-1i*2*pi*freq_offset_est.res*repmat((0:1:(size(y_rx,1)-1))',1,UE.nRX)/LTE_params.Nfft);
end