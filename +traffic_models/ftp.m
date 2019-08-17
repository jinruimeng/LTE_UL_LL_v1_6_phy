classdef ftp < traffic_models.generic_tm
% This class is used for ftp simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC 

properties
    mean_file = 2*10^6; % mean file size before truncation bytes
    max_file = 5*10^6;  % maximum file size (truncation value)  bytes
    sigma_file = 0.722*10^6; % file size standard deviation before truncation bytes
    mean_reading_time = 180; % reading time in sec.
    data_cmf
    reading_cmf
    data_x
    reading_x
    waiting_time
    state
end

methods
    function obj = ftp
        obj = obj@traffic_models.generic_tm;
        obj.type = 'ftp';
        obj.ID_count_current = 1;
        obj.ID_count_next = 1;
        
        % Log normally distributed file size
        obj.data_x = linspace(obj.max_file/10000,obj.max_file,10000); 
        sigma1 = log(1+obj.sigma_file/obj.mean_file^2);
        mu1 = log(obj.mean_file)-0.5*log(1+obj.sigma_file/obj.mean_file^2);
        pmf = 1./sqrt(2*pi*sigma1*obj.data_x.^2).*exp(-(log(obj.data_x)-mu1).^2/(2*sigma1))*(obj.data_x(end)-obj.data_x(end-1)); % according to RAN R1-070674 (mean = 2 Mbyte)
        pmf(end) = pmf(end)+1-sum(pmf);  % renormalization due to truncation to max_file
        obj.data_cmf = [0,cumsum(pmf)];
        
        % Exponentially distributed reading time
        lambda = 1/obj.mean_reading_time;
        percentile9995 = 1/lambda + 1/lambda*10;
        obj.reading_x = linspace(percentile9995/10^4,percentile9995,10^4);
        pmf = lambda*exp(-lambda*obj.reading_x)*(obj.reading_x(end)-obj.reading_x(end-1));
        pmf(end) = pmf(end)+1-sum(pmf);  % renormalization due to truncation to max_file
        obj.reading_cmf = [0,cumsum(pmf)];
        
        obj.generate_file(0);
    end
    
    function generate_file(obj,current_TTI)     
        generate_packet(obj,obj.eval_cmf(obj.data_cmf,obj.data_x),current_TTI);
        obj.state = true; % we are in data transfer state
    end
    
    function check_deadline(obj,TTI)    % check wheter the deadline is met
        check_deadline@traffic_models.generic_tm(obj);
        if ~obj.state
            obj.waiting_time = obj.waiting_time-1;
            if obj.waiting_time == 0  % if the end of the reading state is reached
                obj.generate_file(TTI);
            end
        end
    end
    
     function decrease_data_buffer(obj,N_data_bits,ACK,current_TTI,UE_output) % decrease the data buffer by the amount of acknowledged data
        decrease_data_buffer@traffic_models.generic_tm(obj,N_data_bits,ACK,current_TTI);
        if obj.data_buffer(obj.ID_count_current) <= 0  &&  obj.state  % when all data is transmitted switch to reading state
            obj.waiting_time = obj.eval_cmf(obj.reading_cmf,obj.reading_x);
            obj.waiting_time = round(obj.waiting_time/10^-3); % transform from sec. to TTIs
            obj.state = false;  % enter the reading state
        end
        UE_output.type=obj.type;
        UE_output.data_buffer_left=obj.data_buffer;
        UE_output.data_generated=obj.data_generated;
        UE_output.TTI_origin=obj.TTI_origin;
        UE_output.delay_TTI=obj.delay_TTI;
        UE_output.ID_count_current=obj.ID_count_current;
        UE_output.ID_count_next=obj.ID_count_next;
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