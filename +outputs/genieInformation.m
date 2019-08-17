classdef genieInformation < handle
% Genie information. Including transmitted bits and so on.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

   properties
       data_bits
       bits_to_turboencode
       sent_bits
       y_tx_assembled       % Contains the modulated sent symbols before the padding and IFFT
       y_tx                 % contains the TX output after padding and IFFT
       v                    % added noise at RX-antennas (preFFT, y = H*x + v) with respect to Nfft and Ntot
       n                    % noise after FFT at the receiver (postFFT), where v = sqrt(LTE_params.Nfft/LTE_params.Ntot) * n
       rs                   % reference signal for uplink
       rs_base_seq          % reference signal base sequence
       CSFiDCIFormat
       CH_mapping           % allocation of UE data
       DM_mapping           % allocation of reference signals 
       channel_estimate     % save channel estimate to use it for interpolation in the next subframe
       channel_estimate_complete
       perfect_channel      % save perfect channel for link adaptation with delay
   end

   methods
       function clear(obj)
           obj.data_bits = [];
           obj.bits_to_turboencode = [];
           obj.sent_bits = [];
       end
   end
end 
