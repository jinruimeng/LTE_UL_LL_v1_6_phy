classdef ueOutput < handle
% UE output, including ACK, subframe rv_idx, CQI...
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

   properties
       ACK
       UE_scheduled
       rv_idx
       CQI_feedback
       rx_data_bits
       rx_coded_bits       
       RI
       PMI
       CQI
       CQI_bar
       HARQ_process
       UL_HARQ_process % Prokopec
       channel_estimation_error
       channel_prediction_error
       % Kejik - B1/2
       PCFICH_rx 
       PCFICH_CFI 
       PHICH_HI_rx
       y_tx 
       UE_genie
       papr
       LayerMapping 
       
       %Traffic
       type
       data_buffer_left
       data_generated
       TTI_origin
       delay_TTI
       ID_count_current
       ID_count_next

   end

   methods
       function obj = ueOutput(nUE)
           obj.UE_genie = outputs.genieInformation;
%            for u_=1:nUE
%                obj.UE_genie(u_) = outputs.genieInformation;
%            end
       end  
           
       function clear(obj)
           obj.ACK = [];
           obj.UE_scheduled = [];
           obj.rv_idx = [];
           obj.CQI_feedback = [];
           obj.rx_data_bits = [];
           obj.rx_coded_bits = [];
           obj.RI = [];
           obj.PMI = [];
           obj.HARQ_process = [];
           obj.UL_HARQ_process = []; % prokopec
           obj.channel_estimation_error = [];
           obj.CQI = [];
           % Kejik - B2/2
           obj.PCFICH_rx = [];
           obj.PCFICH_CFI = [];
           obj.PHICH_HI_rx = [];
           % Kejik - E2/2
           
           %traffic
%            obj.data_buffer_left
%            obj.data_generated
%            obj.TTI_origin
%            obj.delay_TTI
%            obj.ID_count_current
%            obj.ID_count_next

       end
   end
end 
