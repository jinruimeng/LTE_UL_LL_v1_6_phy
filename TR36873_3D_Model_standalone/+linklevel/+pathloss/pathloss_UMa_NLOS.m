function pathloss_dB = pathloss_UMa_NLOS(f_c, h_BS, h_UT, d_2D, varargin)
%pathloss_UMa_NLOS - Calculates pathloss in dB for UMa NLOS
%TR 36.873 Table 7.2-1, Urban Macro Non Line Of Sight
%
% Syntax:  
%  pathloss_dB = pathloss_UMa_NLOS(f_c, h_BS, h_UT, d_2D)
%  pathloss_dB = pathloss_UMa_NLOS(f_c, h_BS, h_UT, d_2D, rand_seed)
%  pathloss_dB = pathloss_UMa_NLOS(f_c, h_BS, h_UT, d_2D, rand_seed, h, W)
%  pathloss_dB = pathloss_UMa_NLOS(f_c, h_BS, h_UT, d_2D, h, W)
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
%  pathloss = pathloss_UMa_NLOS(2.1e9, 25, 1.5, 100)

% Author: Markus Gasser, markus.gasser@nt.tuwien.ac.at
% http://www.nt.tuwien.ac.at
% (c) 2016 ITC, TU Wien

% Set RNG
if numel(varargin) == 1 || numel(varargin) == 3
    rand_state = rng;
    rng(varargin{1});
end

% Set street width and building height
if numel(varargin) == 2 || numel(varargin) == 3
    if numel(varargin) == 2
        h = varargin{1};
        W = varargin{2};
    else
        h = varargin{2};
        W = varargin{3};
    end
    
    % Check conditions
    if h < 5 || h > 50
        warning('Building height not within bounds');
    end
    if W < 5 || W > 50
        warning('Street width not within bounds');
    end
else
    W = 20;
    h = 20;
end

% Calculate Pathloss
d_3D = sqrt(d_2D^2 + (h_BS-h_UT)^2);
pathloss_LOS = linklevel.pathloss.pathloss_UMa_LOS(f_c, h_BS, h_UT, d_2D);
pathloss_NLOS = 161.04 - 7.1*log10(W) + 7.5*log10(h) - ...
    (24.37-3.7*(h/h_BS)^2)*log10(h_BS) + ...
    (43.42-3.1*log10(h_BS))*(log10(d_3D)-3) + 20*log10(f_c*1e-9) - ...
    - (3.2*(log10(17.625))^2-4.97) - 0.6*(h_UT-1.5);
pathloss_dB = max(pathloss_LOS, pathloss_NLOS);

% Restore RNG
if exist('rand_state', 'var')
    rng(rand_state);
end

end