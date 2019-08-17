% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function [ pathloss_matrix ] = generate_pathloss_matrix(connection_table, loss)
%generate_pathloss_matrix generates a standardised pathloss_matrix
% the additional loss is always "loss" (dB)
%
% inputs: connection_table ... which UE is connected to what BS
%         loss             ... additional loss in dB to non connected BSs
%
% output: pathloss_matrix  ... specifies additional loss in dB

nBS = size(connection_table, 1);
nUE = size(connection_table, 2) / nBS;
pathloss_matrix = zeros(nBS, nUE);

for bb = 1:nBS
    for uu = 1:nBS*nUE
        if utils.isUEattached(connection_table, bb, uu) == false
            pathloss_matrix(bb,uu) = loss;
        end
    end
end

end

