classdef ueSignaling < handle
% DL-SCH signaling for each UE.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
% UL-SCH signaling for UPLINK,
% Jan Prokopec, prokopec@feec.vutbr.cz

   properties
       % Channel coding
%       turbo_rate_matcher  % Signaling related to turbo rate matching
       turbo_rate_matcher_UL
%       TB_size             % TB size
       TB_size_UL
%       TB_segmentation     % Information about TB segmentation
       TB_segmentation_UL
%       turbo_encoder       % Turbo encoder related info
       turbo_encoder_UL
       
       % Scheduling
%	   MCS_and_scheduling % Information related to modulation and coding (CQI-related), scheduling and precoding
       MCS_and_scheduling_UL % UPLINK Information related to modulation and coding (CQI-related), scheduling and precoding

   end

   methods
       function clear(obj)
           obj.turbo_rate_matcher = [];
           obj.TB_size = [];
           obj.TB_segmentation = [];
           obj.turbo_encoder = [];
           obj.turbo_rate_matcher_UL = [];
           obj.TB_size_UL = [];
           obj.TB_segmentation_UL = []; 
           obj.turbo_encoder_UL = [];
%           obj.MCS_and_scheduling = [];
           obj.MCS_and_scheduling_UL = [];
       end
   end
end 
