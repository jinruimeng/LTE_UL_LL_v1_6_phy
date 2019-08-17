classdef clock < handle
% This object represents a network-wide clock that all network elements
% share.
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       current_block  % Current TTI (an integer)
       BF_time      % How long does a TTI last (in seconds)
       time         % Actual time (seconds)
   end

   methods
       % Class constructor
       %   - Input
       function obj = clock(BF_time)
           obj.current_block = 0;
           obj.time        = 0;
           obj.BF_time    = BF_time;
       end
       
       % Advances network TTI by one
       function advance_1_TTI(obj)
           obj.current_block = obj.current_block + 1;
           obj.time = obj.time + obj.BF_time;
       end
   end
end 
