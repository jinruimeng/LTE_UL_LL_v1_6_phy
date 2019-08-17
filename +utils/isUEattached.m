% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function [ isAttached ] = isUEattached( connection_table, bb, uG )
%isUEattached checks if user uG is attached to basestation bb
%
% inputs: connection_table (LTE_params.connection_table)
%         bb ... basestation index
%         uG ... "global user index"
%
% output: uL ... the local user index
%

isAttached = false;

if connection_table(bb,uG) == 1
    isAttached = true;
end

end

