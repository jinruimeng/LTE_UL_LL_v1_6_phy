classdef gaming < traffic_models.generic_tm
% This class is used for gaming simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC
properties
    UDP_header = 2; % header overhead
    next_packet;  
    initial_packet_a = 0;
    initial_packet_b = 40;
    initial_packet_size_a = 45;
    initial_packet_size_b = 5.7;
    packet_time_a = 55;
    packet_time_b = 6;
    packet_size_a = 120;
    packet_size_b = 36;
    delay_constraint = 60;
%     arrival_rate = (120+36*-psi(1))*8/(55+6*-psi(1));
end

methods
    function obj = gaming
        obj = obj@traffic_models.generic_tm;
        obj.type = 'gaming';
        obj.ID_count_current = 1;
        obj.ID_count_next = 1;
        obj.next_packet = randi([obj.initial_packet_a,obj.initial_packet_b],1); % intial packet of the game 
        if obj.next_packet == 0
            coin_toss = rand;
            generate_packet(obj,round(obj.initial_packet_size_a-obj.initial_packet_size_b*log(-log(coin_toss))+obj.UDP_header),0);
            obj.generate_iat;
        end
    end
    
    function generate_iat(obj) % inter arrival time
        coin_toss = rand;
        obj.next_packet = round(obj.packet_time_a-obj.packet_time_b*log(-log(coin_toss)));
    end
    
    function check_deadline(obj,current_TTI)    % check wheter the deadline is met
        obj.next_packet = obj.next_packet-1;
        if ~obj.next_packet
            coin_toss = rand;
            generate_packet(obj,round(obj.packet_size_a-obj.packet_size_b*log(-log(coin_toss))+obj.UDP_header),current_TTI);
            obj.generate_iat;
        end
    end
end
end