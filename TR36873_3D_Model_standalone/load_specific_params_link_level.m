function config = load_specific_params_link_level(config)
% Load specific channel model parameters from Table 7.3-6 and system
% parameters for the FFT of the channel impulse response
%
% (c) Fjolla Ademaj, Martin Taranetz, ITC 2016


switch config.channel_model.type
    case '3D_UMi_fading'
        config = load_UMi_parameters(config);
    case '3D_UMa_fading'
        config = load_UMa_parameters(config);
    otherwise
        error('Channel model type not valid')
end


switch config.bandwidth
    case 1.4e6
        config.N_RB = 6;
        config.fft_points = 128;
        if strcmp(config.CP_length,'normal')
            config.CP_length_samples = 9;
        else
            config.CP_length_samples = 32;
        end
    case 3e6
        config.N_RB = 15;
        config.fft_points = 256;
        if strcmp(config.CP_length,'normal')
            config.CP_length_samples = 18;
        else
            config.CP_length_samples = 64;
        end
    case 5e6
        config.N_RB = 25;
        config.fft_points = 512;
        if strcmp(config.CP_length,'normal')
            config.CP_length_samples = 36;
        else
            config.CP_length_samples = 128;
        end
    case 10e6
        config.N_RB = 50;
        config.fft_points = 1024;
        if strcmp(config.CP_length,'normal')
            config.CP_length_samples = 72;
        else
            config.CP_length_samples = 256;
        end
    case 15e6
        config.N_RB = 75;
        config.fft_points = 1536;
        if strcmp(config.CP_length,'normal')
            config.CP_length_samples = 108;
        else
            config.CP_length_samples = 384;
        end
    case 20e6
        config.N_RB = 100;
        config.fft_points = 2048;
        if strcmp(config.CP_length,'normal')
            config.CP_length_samples = 144;
        else
            config.CP_length_samples = 512;
        end
    otherwise
        error('Bandwidth not supported');
end
config.Ntot = config.N_RB*12;
config.fs = 15e3*config.fft_points;
config.RB_bandwidth    = 180e3;
% config.sample_points = (0:config.fft_points-1)/config.fs;
config.sample_points = (0:config.CP_length_samples-1)/config.fs;

config.add_femtocells = false;
config.N_sym = 14;
end