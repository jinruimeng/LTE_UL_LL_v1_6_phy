function MI_data = MIdata_load


% 4QAM data
load('./+network_elements/MI_data/BICMC_4.mat','snr','I_mp');
MI_data(1).SNR = 10.^(linspace(snr(1),snr(end),1000)/10);
MI_data(1).BICM = interp1(snr,I_mp,10*log10(MI_data(1).SNR));
MI_data(1).max_SNR = max(MI_data(1).SNR);
MI_data(1).min_SNR = min(MI_data(1).SNR);
[~,MI_data(1).sat_SNR] = max(MI_data(1).BICM);

% 16QAM data
load('./+network_elements/MI_data/BICMC_16.mat','snr','I_mp');
MI_data(2).SNR = 10.^(linspace(snr(1),snr(end),1000)/10);
MI_data(2).BICM = interp1(snr,I_mp,10*log10(MI_data(2).SNR));
MI_data(2).max_SNR = max(MI_data(2).SNR);
MI_data(2).min_SNR = min(MI_data(2).SNR);
I = (diff(MI_data(2).BICM) < 0.000001);
[~,MI_data(2).sat_SNR] = max(I);
% [~,MI_data(2).sat_SNR] = max(MI_data(2).BICM);

% 64QAM data
load('./+network_elements/MI_data/BICMC_64.mat','snr','I_mp');
MI_data(3).SNR = 10.^(linspace(snr(1),snr(end),1000)/10);
MI_data(3).BICM = interp1(snr,I_mp,10*log10(MI_data(3).SNR));
MI_data(3).max_SNR = max(MI_data(3).SNR);
MI_data(3).min_SNR = min(MI_data(3).SNR);
[~,MI_data(3).sat_SNR] = max(MI_data(3).BICM);

% lines through BICM at different CQI values
% load('./+network_elements/MI_data/BICM_k_d.mat','k','d');
load('./+network_elements/MI_data/BICM_k_d_MSE.mat','k','d');
% load('./+network_elements/MI_data/Efficiency_k_d_MSE.mat','k','d');
MI_data(1).k = k;
MI_data(1).d = d;
