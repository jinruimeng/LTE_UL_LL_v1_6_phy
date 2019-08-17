% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function [ uG ] = globalToLocalUser( connection_table, bb, uL )
%GLOBALTOLOCALUSER returns the global user index 

% inputs: connection_table (LTE_params.connection_table)
%         bb ... basestation index
%         uL ... "local user index"
%
% output: uL ... the global user index


% extract the right row
connected_users = connection_table(bb,:);



% compute the relative user index
uG = find(cumsum(connected_users) == uL, 1, 'first');

end

