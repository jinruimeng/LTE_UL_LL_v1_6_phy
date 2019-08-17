classdef eNodeB < handle
% Class that represents an eNodeB.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

    properties
        NIDcell                             % Cell identity
        UE_specific                         % UE specific data, such as the HARQ TX buffer
        user_count                          % How many users there are
        channel_estimation_method           % Channel estimation method
        channel_interpolation_method        % channel interpolation method
        channel_interpolation_past_points   % numberof previous points for interpolation
        autocorrelation_matrix_type         % type of autocorrelation amtrix('ideal','estiamted')
        realization_num                     % number of realizations of channel, used for averaging fo channel autocorrelation matrix
        HARQ_rx_soft_buffer                 % HARQ soft buffer. Stores codewords, not data bits!!!!
        demapping_method                    % Demapping method to use
        turbo_iterations                    % Number of iterations of the turbo decoder
        realization_num_total               % first xy number of channel realizations are used just for estimation of autocorrelation matrix
        nRX                                 % Number of receive antennas
        scheduler                           % Where we will store this eNodeB's scheduler
        receiver                            % receiver that the UE applys
        AtPort                              % supported antenna ports 
        nAtPort                             % number of antenna ports
        clock                               % So the eNodeB is aware in which TTI he is
        HARQ_process_count                  % Number of HARQ processes
    end

   methods
       % Class constructor
       function obj = eNodeB(NID,nRX,nUE,AtPort,nAtPort,maxStreams,HARQ_processes,max_rv_idx,BS_params)
           obj.NIDcell = NID;
           obj.UE_specific = network_elements.UeSpecificEnodebData(maxStreams,HARQ_processes,max_rv_idx);
           for uu=1:nUE
               obj.UE_specific(uu) = network_elements.UeSpecificEnodebData(maxStreams,HARQ_processes,max_rv_idx);
           end
           obj.nRX                                  = nRX;
           obj.user_count                           = nUE;
           obj.AtPort                               = AtPort;
           obj.nAtPort                              = nAtPort;
           obj.HARQ_process_count                   = HARQ_processes;
           obj.channel_estimation_method            = BS_params.channel_estimation_method;
           obj.channel_interpolation_method         = BS_params.channel_interpolation_method;
           obj.channel_interpolation_past_points    = BS_params.channel_interpolation_past_points;
           obj.autocorrelation_matrix_type          = BS_params.autocorrelation_matrix_type;
           obj.realization_num                      = BS_params.realization_num;
           obj.turbo_iterations                     = BS_params.turbo_iterations;
           obj.realization_num_total                = BS_params.realization_num_total;
           obj.receiver                             = BS_params.receiver;
       end
       
       % Advance the HARQ process that is being used
       function update_current_HARQ_process(obj,UE_output)
           
           for u_=1:obj.user_count
               
               max_codewords = size(obj.UE_specific(u_).HARQ_processes,1);
               
               % empty HARQ process means one of the first ones, where no ACKs are present.
               if isempty(UE_output(u_).HARQ_process)
                   % Set the current HARQ process as the first free one
                   for cw_=1:max_codewords
                       if obj.UE_specific(u_).current_HARQ_process(cw_).assigned_RBs > 0
                           HARQ_idx(cw_) = find(obj.UE_specific(u_).free_HARQ_processes(cw_,:),1,'first');
                           obj.UE_specific(u_).free_HARQ_processes(cw_,HARQ_idx(cw_)) = false;
                           obj.UE_specific(u_).current_HARQ_process(cw_) = obj.UE_specific(u_).HARQ_processes{cw_,HARQ_idx(cw_)};
                           obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx = 0;
                       else
                           HARQ_idx(cw_) = obj.UE_specific(u_).current_HARQ_process(cw_).id;
                           obj.UE_specific(u_).free_HARQ_processes(cw_,HARQ_idx(cw_)) = false;
                       end
                   end
               else
                   % Mark HARQ process from which the current ACK has been received as free
                   UE_scheduled      = UE_output(u_).UE_scheduled;
                   feedback_HARQ_idx = [ UE_output(u_).HARQ_process(:);zeros(max_codewords-length(UE_output(u_).HARQ_process(:)),1)]; % Processes from which the current ACK comes
                   
                   for cw_=1:max_codewords
                       
                       % Mark newly-free HARQ processes
                       if UE_scheduled
                           if feedback_HARQ_idx(cw_)~=0
                               obj.UE_specific(u_).free_HARQ_processes(cw_,feedback_HARQ_idx(cw_)) = true;
                           end
                       end
                       
                       % TX HARQ process
                       % If the current HARQ process was used (ie. scheduled), move the current HARQ process forward and update rv_idx
                       % Otherwise do nothing
                       if obj.UE_specific(u_).current_HARQ_process(cw_).assigned_RBs > 0
                           next_free_HARQ_idx(cw_) = find(obj.UE_specific(u_).free_HARQ_processes(cw_,:),1,'first');
                           obj.UE_specific(u_).current_HARQ_process(cw_) = obj.UE_specific(u_).HARQ_processes{cw_,next_free_HARQ_idx(cw_)};
                           current_rv_idx = obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx;
                           max_rv_idx     = obj.UE_specific(u_).current_HARQ_process(cw_).max_rv_idx;
                           if isempty(current_rv_idx)
                               obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx = 0;
                           else
                               if UE_output(u_).ACK(cw_)
                                   obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx = 0;
                               else
                                   obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx = mod(current_rv_idx+1,max_rv_idx+1);
                               end
                           end
                           obj.UE_specific(u_).free_HARQ_processes(cw_,next_free_HARQ_idx(cw_)) = false;
                       else
                           % Current HARQ process is still the same
                       end
                   end
               end
           end
       end
       
       % Reset the HARQ process index
       function reset_HARQ_process_index(obj)
           for u_=1:obj.user_count
               % Set the HARQ processes for all the codewords
               obj.UE_specific(u_).current_HARQ_process = [obj.UE_specific(u_).HARQ_processes{:,1}];
               % Reset the rest of the HARQ processes
               for cw_=1:size(obj.UE_specific(u_).HARQ_processes,1)
                   for harq_=1:size(obj.UE_specific(u_).HARQ_processes,2)
                       obj.UE_specific(u_).HARQ_processes{cw_,harq_}.reset;
                       obj.UE_specific(u_).free_HARQ_processes(cw_,:) = ones(size(obj.UE_specific(u_).free_HARQ_processes(cw_,:)));
                   end
               end
           end
       end
   end
end 
