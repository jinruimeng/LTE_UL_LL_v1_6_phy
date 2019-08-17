function P_LOS = LOS_probability_UMa(d_2D_out, h_UT)
%LOS_probability_UMa - Calculates Line Of Sight probability for UMa
%TR 36.873 Table 7.2-2, Urban Macro
%
% Syntax:
%  P_LOS = LOS_probability_UMa(d_2D_out, h_UT)
%
% Inputs:
%  d_2D_out - 2D (outdoor) distance between base station and user [m]
%  h_UT - Height of user [m]
%
% Outputs:
%  P_LOS - Probability of the connection being LOS
%
% Example: 
%  P_LOS = LOS_probability_UMa(100, 1.5)

% Author: Markus Gasser, markus.gasser@nt.tuwien.ac.at
% http://www.nt.tuwien.ac.at
% (c) 2016 ITC, TU Wien
P_LOS = (min(18./d_2D_out, 1).*(1-exp(-d_2D_out./63))+ ...
    exp(-d_2D_out./63)) .* (1+C_func(d_2D_out, h_UT));
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
    g = (1.25*exp(1)-6) .* d_2D.^3 .* exp(-d_2D/150);
    g_act = (d_2D > 18);
    C_val = ((h_UT-13)/10).^1.5 .* g .* g_act;
end

end