function W = LTE_common_get_precoding_matrix(tx_mode,nAtPort,codebook_index,nLayers,LTE_params)
% This function returns the precoding matrix for a specific transmission mode
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

W = 0;

switch tx_mode
    case {2,3}   % transmit diversity and open loop spatial multiplexing    
        error('not standardized in uplink LTE');
    case 4   % closed loop spatial multiplexing 
        W=LTE_UL_Codebook(nLayers, codebook_index, nAtPort);
end

end