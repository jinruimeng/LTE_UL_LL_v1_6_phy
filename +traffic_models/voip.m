classdef voip < traffic_models.generic_tm
% This class is used for voip simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC

properties
    state
    c = 0.01;  % data according to RAN R1-070674
    d = 0.99;
end

methods
    function obj = voip
        obj = obj@traffic_models.generic_tm;
        obj.type = 'voip';
        obj.data_buffer(1) = 40; % the size of a voice packet (assume we start in the active state)
        obj.TTI_origin(1) = 0;
        obj.ID_count_current = 1;
        obj.ID_count_next = 2;
        obj.state = true;
        obj.deadline = 20;
    end
    
    function check_deadline(obj,current_TTI)    % check wheter the deadline is met
        check_deadline@traffic_models.generic_tm(obj);
%         obj.deadline
        if (obj.deadline==0)    % here comes the two state markov model assumed for the speech process (state changes every 20 ms)
            coin_toss = rand;
            if obj.state    % active state
               if coin_toss < obj.d
                   obj.state = true;
               else
                   obj.state = false;
               end
            else        % inactive state
                if coin_toss < obj.c
                    obj.state = true;
                else
                    obj.state = false;
                end
            end
            
            if obj.state % generate a new speech packet in the active state 
                generate_packet(obj,40,current_TTI);
                obj.deadline = 20;
            else        % in the inactive state
                if ~obj.deadline % generate a new silence descriptor (SID) packet if the deadline is zero
                    generate_packet(obj,15,current_TTI);
                    obj.deadline = 160;
                end
            end
        end
    end
end

end