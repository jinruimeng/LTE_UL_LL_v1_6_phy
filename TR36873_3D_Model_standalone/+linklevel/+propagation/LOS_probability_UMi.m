function P_LOS = LOS_probability_UMi(d_2D_out)
%LOS_probability_UMi - Calculates Line Of Sight probability for UMi
%TR 36.873 Table 7.2-2, Urban Micro
%
% Syntax:
%  P_LOS = LOS_probability_UMi(d_2D_out)
%
% Inputs:
%  d_2D_out - 2D (outdoor) distance between base station and user [m]
%
% Outputs:
%  P_LOS - Probability of the connection being LOS
%
% Example: 
%  P_LOS = LOS_probability_UMi(100)

% Author: Markus Gasser, markus.gasser@nt.tuwien.ac.at
% http://www.nt.tuwien.ac.at
% (c) 2016 ITC, TU Wien
P_LOS = (min(18./d_2D_out, 1).*(1-exp(-d_2D_out./36))+ exp(-d_2D_out./36));

end