function [ ICI_power ] = LTE_UL_estimate_ICI_power(LTE_params)
% ICI power according to;
% Yang-Seok Choi and Peter J. Voltz and Frank A. Cassara; "On Channel
% Estimation and Detection for Multicarrier Signals in Fast and Selective
% Rayleigh fading Channels". IEEE Trans. on Communications, Vol. 49, No. 8,
% August 2001


fd = LTE_params.UE_config.user_speed/LTE_params.speed_of_light*LTE_params.carrier_freq_UP; % Doppler frequency 

T_bar = 1e-3/LTE_params.Nsub; % OFDM symbol duration (approximately...)
bessel_store = zeros(LTE_params.Ntot,1);
for kk = 1:LTE_params.Nfft
    bessel_store(kk) = besselj(0,2*pi*fd*kk*T_bar/LTE_params.Nfft);
end

ICI_power = zeros(1, LTE_params.Ntot);


for ii = 1:LTE_params.Ntot 
    sum_val = 0;
    for pp = 1:LTE_params.Ntot
        if pp ~= ii
            sum_val = sum_val + LTE_params.Nfft;
            temp_sum = 0;
            for kk = 1:LTE_params.Nfft                    
               temp_sum = temp_sum + (LTE_params.Nfft-kk)*cos(2*pi*kk*(pp-ii)/LTE_params.Nfft)*bessel_store(kk);
            end
            sum_val = sum_val + 2*temp_sum;
        end
    end
    sum_val = 1/LTE_params.Nfft^2*sum_val;
    ICI_power(ii) = sum_val;
end

end