% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function [ uL ] = globalToLocalUser( connection_table, bb, uG )
%GLOBALTOLOCALUSER returns the user index inside a basestation
%   needed if a loop iterates over nUE*nBS to map the global index
%   uG out of nUE*nBS to uL out of 1:nUE according to the current
%   basestation.
% inputs: connection_table (LTE_params.connection_table)
%         bb ... basestation index
%         uG ... "global user index"
%
% output: uL ... the local user index
%

if connection_table(bb,uG) == 0 
    % this should never happen...
    error('user not connected to the selected basestation');
end

% extract the right row
connected_users = connection_table(bb,:);

% compute the relative user index
uL = sum(connected_users(1:uG));

end

