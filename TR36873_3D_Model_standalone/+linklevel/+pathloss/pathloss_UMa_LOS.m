function pathloss_dB = pathloss_UMa_LOS(f_c, h_BS, h_UT, d_2D, varargin)
%pathloss_UMa_LOS - Calculates pathloss in dB for UMa LOS
%TR 36.873 Table 7.2-1, Urban Macro Line Of Sight
%
% Syntax:  
%  pathloss_dB = pathloss_UMa_LOS(f_c, h_BS, h_UT, d_2D)
%  pathloss_dB = pathloss_UMa_LOS(f_c, h_BS, h_UT, d_2D, rand_seed)
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
%  pathloss = pathloss_UMa_LOS(2.1e9, 25, 1.5, 100)

% Author: Markus Gasser, markus.gasser@nt.tuwien.ac.at
% http://www.nt.tuwien.ac.at
% (c) 2016 ITC, TU Wien


% Check conditions
if h_BS ~= 25
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
C_prob = 1/(1+C_func(d_2D, h_UT));
if rand < C_prob
    h_E = 1;
else
    h_UT_round = max(12, round((h_UT-1.5)/3)*3);
    h_E_possible = 12:3:h_UT_round;
    h_E = h_E_possible(randi(length(h_E_possible)));
end

h_BS_eff = h_BS - h_E;
h_UT_eff = h_UT - h_E;
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

% Function C (and g) from Table 7.2-2
function C_val = C_func(d_2D, h_UT)

% Check Conditions
if h_UT > 23
    warning('UE height not within bounds');
end

% Calculate C
if h_UT < 13
    C_val = 0;
else
    if d_2D > 18
        g = (1.25*exp(1)-6) * d_2D^3 * exp(-d_2D/150);
    else
        g = 0;
    end
    C_val = ((h_UT-13)/10)^1.5 * g;
end
    
end