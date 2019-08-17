function [H_cal] = LTE_UL_estimate_ISI_power(ChanMod_config,Tg,SamplingTime,Nfft,Ntot,channel,uu)
%% Werner Henkel and Georg Tauboeck; "The Cyclic Prefix of OFDM/DMT - An Analysis"
% calculating ISI and ICI power caused by insufficient guard interval
% length in terms of a PSD which can then be used in the feedback function


% get impulse responses for all possible paths

% if strcmp(ChanMod_config.filtering,'BlockFading')
%     channel_h = ChanMod_output{uu}.H;
% else% FastFading
%     temp= ChanMod_output{uu}.H;
%     channel_h=mean(temp,3);      %take the mean of all LTV channels
%     channel_h=reshape(channel_h,size(channel_h,1),size(channel_h,2),[]);
% end


channel_temp = circshift(squeeze(mean(channel{uu},2)),[0,-Ntot/2]);
channel_h = ifft(channel_temp,Ntot)./Ntot;

% channel sizes
% N_RX = size(channel_h,1);
% N_TX = size(channel_h,2);
% N_h  = size(channel_h,3);
% N_cp = round(max(Tg./SamplingTime));  % taking the smaller value, which is correct for 6 out of 7 symbols in normal CP

% channel sizes
N_RX = size(channel_h,2);
N_TX = size(channel_h,3);
N_h  = 30;
N_cp = round(max(Tg./SamplingTime));  % taking the smaller value, which is correct for 6 out of 7 symbols in normal CP


%=============================================================
% DFT of the tail of the channel impulse response for every channel entry
H_cal = zeros(N_h-N_cp, Nfft, N_RX, N_TX);
for rr=1:N_RX
    for tt=1:N_TX
        for m=1:(N_h-N_cp)
            H_cal(m,:,rr,tt) = fft(transpose(squeeze(channel_h( (N_cp+m):end, rr, tt))), Nfft);
        end
    end
end
%=============================================================

% PSD according eq.14 from paper above
% ISI_power_PSD = 2.*diag( squeeze(H_cal(:,:,1,1))'*squeeze(H_cal(:,:,1,1)) ) / Nfft;


end

