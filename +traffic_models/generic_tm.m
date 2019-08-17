classdef generic_tm < handle
% This is a generic class for all traffic models
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC

% Last changes 10/2013: Victor Sen, vicsen34@gmail.com

properties
    type
    data_buffer % number of bytes to be transmitted
    data_generated % saves the data generated
    TTI_origin
    delay_TTI

    ID_count_current    % Packet that is already being transmitted
    ID_count_next    % Next packet generate would have that ID
    Max_num_packet  % maximum number of packet to be simulated
    deadline    % deadline in multiples of TTI
    deadline_miss % stores how often the deadline is missed
end

methods
    function obj = generic_tm()
        obj.type = [];
        obj.Max_num_packet = 500;
        obj.data_buffer = zeros(1,obj.Max_num_packet);
        obj.data_generated = zeros(1,obj.Max_num_packet);
        obj.TTI_origin = zeros(1,obj.Max_num_packet);
        obj.delay_TTI = zeros(1,obj.Max_num_packet);
        obj.deadline = Inf;
    end
    
    function decrease_data_buffer(obj,N_data_bits,ACK,current_TTI,UE_output) % decrease the data buffer by the amount of acknowledged data
        if ACK
            obj.data_buffer(obj.ID_count_current) = obj.data_buffer(obj.ID_count_current)-floor(N_data_bits/8);
            if obj.data_buffer(obj.ID_count_current) <= 0 % packet already transmitted
                obj.data_buffer(obj.ID_count_current) = 0; % if due to some bug the buffer gets negative (which should not happen)
                 obj.delay_TTI(obj.ID_count_current) = current_TTI - obj.TTI_origin(obj.ID_count_current);
                obj.ID_count_current = obj.ID_count_current + 1;
            end
        end
        UE_output.type=obj.type;
        UE_output.data_buffer_left=obj.data_buffer;
        UE_output.data_generated=obj.data_generated;
        UE_output.TTI_origin=obj.TTI_origin;
        UE_output.delay_TTI=obj.delay_TTI;
        UE_output.ID_count_current=obj.ID_count_current;
        UE_output.ID_count_next=obj.ID_count_next;
    end
    
    function check_deadline(obj)    % check wheter the deadline is met
        if obj.deadline             % decreases the deadline every TTI
           obj.deadline = obj.deadline-1; 
        end
        if ~obj.deadline && obj.data_buffer(obj.ID_count_current) > 0
            obj.deadline_miss = obj.deadline_miss+1;
        end
    end
    
    function generate_packet(obj,dat_size,current_TTI)
        obj.data_buffer(obj.ID_count_next) = dat_size; % the size of a voice packet
        obj.data_generated(obj.ID_count_next) = dat_size;
        obj.TTI_origin(obj.ID_count_next) = current_TTI;
        obj.ID_count_next = obj.ID_count_next + 1;
    end
    
end

end