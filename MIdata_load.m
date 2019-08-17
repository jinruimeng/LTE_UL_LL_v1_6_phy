function MI_data = MIdata_load

plot_MI = false;

% 4QAM data
load('./+network_elements/MI_data/BICMC_4.mat','snr','I_mp');
MI_data(1).SNR = 10.^(snr/10);
MI_data(1).SNR_dB = snr;
MI_data(1).BICM = I_mp;
MI_data(1).max_SNR = max(MI_data(1).SNR);
MI_data(1).min_SNR = min(MI_data(1).SNR);
MI_data(1).max_SNR_dB = max(MI_data(1).SNR_dB);
MI_data(1).min_SNR_dB = min(MI_data(1).SNR_dB);
[C,MI_data(1).sat_SNR] = max(MI_data(1).BICM);

% 16QAM data
load('./+network_elements/MI_data/BICMC_16.mat','snr','I_mp');
MI_data(2).SNR = 10.^(snr/10);
MI_data(2).SNR_dB = snr;
MI_data(2).BICM = I_mp;
MI_data(2).max_SNR = max(MI_data(2).SNR);
MI_data(2).min_SNR = min(MI_data(2).SNR);
MI_data(2).max_SNR_dB = max(MI_data(2).SNR_dB);
MI_data(2).min_SNR_dB = min(MI_data(2).SNR_dB);
[C,MI_data(2).sat_SNR] = max(MI_data(2).BICM);

% 64QAM data
load('./+network_elements/MI_data/BICMC_64.mat','snr','I_mp');
MI_data(3).SNR = 10.^(snr/10);
MI_data(3).SNR_dB = snr;
MI_data(3).BICM = I_mp;
MI_data(3).max_SNR = max(MI_data(3).SNR);
MI_data(3).min_SNR = min(MI_data(3).SNR);
MI_data(3).max_SNR_dB = max(MI_data(3).SNR_dB);
MI_data(3).min_SNR_dB = min(MI_data(3).SNR_dB);
[C,MI_data(3).sat_SNR] = max(MI_data(3).BICM);

if plot_MI
    figure;
    hold all;
    for i_=1:length(MI_data)
        plot(MI_data(i_).SNR_dB,MI_data(i_).BICM);
    end
    xlabel('SNR [dB]');
    ylabel('BICM Capacity');
    title('BICM capacity');
    grid on;
end