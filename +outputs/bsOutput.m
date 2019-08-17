classdef bsOutput < handle
% Output of the eNodeB. Contains data, signaling information and genie
% data also.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

   properties
   
       UE_signaling_UL  % Signaling for each UE, left in BS, because there is no scheduler in UE
       genie            % Genie information
   end

   methods
       function obj = bsOutput(nUE, nBS, nRB,maxStreams)
           obj.UE_signaling_UL = outputs.ueSignaling;
           obj.genie = outputs.genieInformation;
           for b_ = 1:nBS
               for u_=1:nUE
                   obj.UE_signaling_UL(b_, u_) = outputs.ueSignaling;
                   obj.genie(b_, u_) = outputs.genieInformation;
               end
           end
       end
       
       function clear(obj) 
           obj.y_tx = [];
           obj.y_ul_rx = [];
           for u_=1:length(obj.UE_signaling);
               obj.UE_signaling_UL(u_).clear;
               obj.genie(u_).clear;
           end
       end
   end
end 
