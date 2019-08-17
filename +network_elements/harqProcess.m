classdef harqProcess < handle
% HARQ process that waits for ACKs at the eNodeB
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

   properties
       id
       cqi            % CQI for the Tx and what should be used for the re-tx
       assigned_RBs   % assigned_RBs for the Tx and what should be used for the re-tx
       HARQ_tx_buffer % Buffer where the data bits are stored
       rv_idx         % last rv_idx used
       max_rv_idx     % maximum value rv_idx can take
   end

   methods
       % Class constructor
       function obj = harqProcess(id,maxStreams,max_rv_idx)
           obj.id        = id;
           obj.max_rv_idx = max_rv_idx;
       end
       
       % Reset the HARQ process
       function reset(obj)
           obj.cqi            = [];
           obj.assigned_RBs   = 0;
           obj.rv_idx         = 0;
           obj.HARQ_tx_buffer = [];
       end
   end
end 
