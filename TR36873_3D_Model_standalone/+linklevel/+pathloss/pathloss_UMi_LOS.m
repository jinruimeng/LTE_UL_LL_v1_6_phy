function pathloss_dB = pathloss_UMi_LOS(f_c, h_BS, h_UT, d_2D, varargin)
%pathloss_UMi_LOS - Calculates pathloss in dB for UMi LOS
%TR 36.873 Table 7.2-1, Urban Micro Line Of Sight
%
% Syntax:  
%  pathloss_dB = pathloss_UMi_LOS(f_c, h_BS, h_UT, d_2D)
%  pathloss_dB = pathloss_UMi_LOS(f_c, h_BS, h_UT, d_2D, rand_seed)
%
% Inputs:
%  f_c - Carrier frequency [Hz]
%  h_BS - Height of base station [m]
%  h_UT - Height of user [m]
%  d_2D - 2D distance between base station and user [m]
%  rand_seed - Seed for the random number generator
%
% Outputs:
%  pathloss_dB - Pathloss [dB]
%
% Example: 
%  pathloss = pathloss_UMi_LOS(2.1e9, 25, 1.5, 100)

% Author: Markus Gasser, markus.gasser@nt.tuwien.ac.at
% http://www.nt.tuwien.ac.at
% (c) 2016 ITC, TU Wien


% Check conditions
if h_BS ~= 10
    warning('BS height not within bounds');
end
if h_UT < 1.5 || h_UT > 22.5
    warning('UE height not within bounds');
end

% Set RNG
if numel(varargin) > 0
    rand_state = rng;
    rng(varargin{1});
end

% Calculate pathloss

h_BS_eff = h_BS - 1;
h_UT_eff = h_UT - 1;
d_BP_eff = 4 * h_BS_eff * h_UT_eff * f_c / 3e8;

d_3D = sqrt(d_2D^2 + (h_BS-h_UT)^2);

if d_2D < d_BP_eff
    pathloss_dB = 22*log10(d_3D) + 28 + 20*log10(f_c*1e-9);
else
    pathloss_dB = 40*log10(d_3D) + 28 + 20*log10(f_c*1e-9) - ...
        9*log10( (d_BP_eff^2 + (h_BS-h_UT)^2 ));
end

if d_2D > 5000 || d_2D < 10
    warning('2D distance no within bounds');
end

% Restore RNG
if exist('rand_state', 'var')
    rng(rand_state);
end

end