% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function [ connection_table ] = generate_connection_table( nBS, nUE )
%generate_connection_table generates a standardised connection table
%
% inputs: nBS ... number of basestations
%         nUE ... number of users per basestation
%
% output: connection_table of the following structure:
%
%           1      ...      nUE      1    ...    nUE         ...          1  ...  nUE       
%  1        1       1        1       0     0      0           0               0     
%  2        0       0        0       1     1      1           0               0
%  ...
%  nBS      0       0        0       0     0      0           0           1   1    1
    
connection_table = kron(eye(nBS), ones(1, nUE));

end

