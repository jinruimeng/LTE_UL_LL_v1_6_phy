function pathloss_dB = pathloss_UMi_NLOS(f_c, h_BS, h_UT, d_2D, varargin)
%pathloss_UMi_NLOS - Calculates pathloss in dB for UMi NLOS
%TR 36.873 Table 7.2-1, Urban Micro Non Line Of Sight
%
% Syntax:  
%  pathloss_dB = pathloss_UMi_NLOS(f_c, h_BS, h_UT, d_2D)
%  pathloss_dB = pathloss_UMi_NLOS(f_c, h_BS, h_UT, d_2D, rand_seed)
%
% Inputs:
%  f_c - Carrier frequency [Hz]
%  h_BS - Height of base station [m]
%  h_UT - Height of user [m]
%  d_2D - 2D distance between base station and user [m]
%  rand_seed - Seed for the random number generator
%  h - Building height [m]
%  W - Street width [m]
%
% Outputs:
%  pathloss_dB - Pathloss [dB]
%
% Example: 
%  pathloss = pathloss_UMi_NLOS(2.1e9, 25, 1.5, 100)

% Author: Markus Gasser, markus.gasser@nt.tuwien.ac.at
% http://www.nt.tuwien.ac.at
% (c) 2016 ITC, TU Wien

% Set RNG
if numel(varargin) > 0
    rand_state = rng;
    rng(varargin{1});
end

% Calculate Pathloss
d_3D = sqrt(d_2D^2 + (h_BS-h_UT)^2);
pathloss_LOS = linklevel.pathloss.pathloss_UMi_LOS(f_c, h_BS, h_UT, d_2D);
pathloss_NLOS = 36.7*log10(d_3D) + 22.7 + 26*log10(f_c) - 0.3*(h_UT-1.5);
pathloss_dB = max(pathloss_LOS, pathloss_NLOS);

% Restore RNG
if exist('rand_state', 'var')
    rng(rand_state);
end

end