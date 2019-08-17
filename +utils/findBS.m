% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function [ bb ] = findBS( connection_table, uG )
%FINDBS returns the base station index for a global user
%   
% inputs: connection_table (LTE_params.connection_table)
%         uG ... global user (set of nBS*nUE users)
%
% output: bb ... index of base station the user is attached to 
%


bb = find(connection_table(:,uG));

end

