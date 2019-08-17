function pathloss_dB = pathloss_UMi_OtoI(f_c, h_BS, h_UT, d_2D, d_2D_in, is_LOS, varargin)
%pathloss_UMi_OtoI - Calculates pathloss in dB for UMi O to I
%TR 36.873 Table 7.2-1, Urban Micro Outdoor to Indoor
%
% Syntax:  
%  pathloss_dB = pathloss_UMi_OtoI(f_c, h_BS, h_UT, d_2D, d_2D_in)
%  pathloss_dB = pathloss_UMi_OtoI(f_c, h_BS, h_UT, d_2D, d_2D_in, rand_seed)
%
% Inputs:
%  f_c - Carrier frequency [Hz]
%  h_BS - Height of base station [m]
%  h_UT - Height of user [m]
%  d_2D - 2D distance between base station and user [m]
%  d_2D_in - Indoor 2D distance between user and wall [m] 
%  is_LOS - true means the connection is LOS, false means it is NLOS
%  rand_seed - Seed for the random number generator
%
% Outputs:
%  pathloss_dB - Pathloss [dB]
%
% Example: 
%  pathloss = pathloss_UMi_OtoI(2.1e9, 25, 1.5, 100, 10)

% Author: Markus Gasser, markus.gasser@nt.tuwien.ac.at
% http://www.nt.tuwien.ac.at
% (c) 2016 ITC, TU Wien

if is_LOS
    PL_b = linklevel.pathloss.pathloss_UMi_LOS(f_c, h_BS, h_UT, d_2D, varargin);
else
    PL_b = linklevel.pathloss.pathloss_UMi_NLOS(f_c, h_BS, h_UT, d_2D, varargin);
end
PL_tw = 20;
PL_in = 0.5 * d_2D_in;
pathloss_dB = PL_b + PL_tw + PL_in;