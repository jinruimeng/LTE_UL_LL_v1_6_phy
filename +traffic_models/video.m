classdef video < traffic_models.generic_tm
% This class is used for gaming simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC
properties
    slice_mean = 100;    % mean slice size
    slice_max = 250;    % max slice size in bytes
    slice_min = 20;    % min slice size
    inter_mean = 6;     % mean slice interarrival time (encoder delay)
    inter_max = 12.5;   % max slice interarrival time
    inter_min = 2.5;   % min slice interarrival time
    slice_x
    inter_x
    slices = zeros(8,1); % eight slices per frame
    inters = zeros(7,1);
    slice_cmf
    inter_cmf
    counter = 1;
end

methods
    function obj = video
        obj = obj@traffic_models.generic_tm;
        obj.deadline = 100;   % 100ms deadline
        obj.type = 'video';
        obj.ID_count_current = 1;
        obj.ID_count_next = 1;
        
        % Truncated Pareto distribution for the slice size
        a = 1.2;
        k = obj.slice_min;
        m = obj.slice_max;
        obj.slice_x = linspace(obj.slice_min,obj.slice_max,1000);
        pmf = a*k^a./obj.slice_x.^(a+1)*(obj.slice_x(end)-obj.slice_x(end-1));
        pmf(end) = pmf(end)+(k/m)^a;
        obj.slice_cmf = [0,cumsum(pmf)];
        obj.generate_slices;
        
        % Truncated Pareto distribution for the interarrival time
        a = 1.2;
        k = obj.inter_min;
        m = obj.inter_max;
        obj.inter_x = linspace(obj.inter_min,obj.inter_max,1000);
        pmf = a*k^a./obj.inter_x.^(a+1)*(obj.inter_x(end)-obj.inter_x(end-1));
        pmf(end) = pmf(end)+(k/m)^a;
        obj.inter_cmf = [0,cumsum(pmf)];
        obj.generate_inters;
        
        % fill the buffer
        generate_packet(obj,obj.slices(1),0);
        obj.slices(1) = 0; % data is now in the buffer
    end
    
    function generate_slices(obj)
        for i = 1:length(obj.slices) 
            obj.slices(i) = obj.eval_cmf(obj.slice_cmf,obj.slice_x);
        end
    end
    
    function generate_inters(obj)
        for i = 1:length(obj.inters) 
            obj.inters(i) = obj.eval_cmf(obj.inter_cmf,obj.inter_x);
        end
    end
    
    function check_deadline(obj,current_TTI)    % check wheter the deadline is met
        check_deadline@traffic_models.generic_tm(obj);
        obj.counter = obj.counter+1;
        if ~obj.deadline % deadline is over -> generate new data
            if (sum(obj.slices)~=0)     % if there are still slices left put them in the buffer now, because a new frame is generated.
                generate_packet(obj,sum(obj.slices),current_TTI);
            end
            obj.generate_slices;
            obj.generate_inters;
            generate_packet(obj,obj.slices(1),current_TTI);
            obj.slices(1) = 0; % data is now in the buffer
            obj.deadline = 100;
            obj.counter = 0;
        else % check whether a new slice must be sent
            iit = cumsum(obj.inters);
            ind = find(obj.counter == iit, 1);
            ind = ind+1; % first slice is not included in obj.inters
            if ~isempty(ind)
                generate_packet(obj,obj.slices(ind),current_TTI);
                obj.slices(ind) = 0;
            end
        end        
    end
    
    function X=eval_cmf(obj,cmf,x)
        coin_toss = rand;
        while coin_toss==1      % we have to add this loop bcs when coin_toss==1 larger is a vector of zeros and an error appears
            coin_toss = rand;
        end
        larger = cmf/coin_toss>ones(size(cmf));
        [C,I] = max(larger);
        toss = rand;
        if toss > 0.5   % instead of just rounding
            X = ceil(x(I-1)); 
        else
            X = floor(x(I-1));
        end
    end
end
end