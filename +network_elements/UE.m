classdef UE < handle
% Class that represents an LTE UE.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

properties
    user_speed                      % paramter for fast fading channel
    mode                            % DEFINED IN STANDARD 3GPP TS 36.213-820 Section 7.1, page 12
                                    % 1: Single Antenna, 4: Closed Loop SM
    nTX                             % number of transmit antennas at UE
    clock = 0;                      % So the UE is aware in which TTI he is (initialized to 0)
    channel_autocorrelation_matrix  % channel autocorrelation matrix
    channel_coef_rosa_zheng         % phi, theta, psi for Rosa Zhneg Channel Model
    previous_channels               % necessary for channel extrapolation
    PMI_fb_gran                     % PMI feedback granularity in multiples of resource blocks
    PMI_fb                          % wheter PMI feedback is activated
    RI_fb                           % wheter RI  feedback is activated
    predict                         % wheter channel prediction is activated
    CQI_fb_gran                     % CQI feedback granularity 
    CQI_fb                          % wheter CQI feedback is activated
    SINR_averager                   % defines the SINR averaging method used for feedback calculation
    MSE_buffer                      % channel estimation MSE buffer
    ch_est_buff                     % buffer of channel estimates
    rx_refsym_buff                  % buffer of received reference symbols
    refsym_buff                     % buffer of reference symbols
    traffic_model                   % traffic model used for this UE
    N_soft                          % Defines the total number of soft channel bits available for HARQ processing (TS 36.306 4.2.1.3)
    
    CSFiDCIFormat
    Isrs                            % SRS configuration index
    SRS_triggerType
    SRS_duration                    % single or periodic
    Tsrs                            % SRS periodicity
    Toffset                         % SRS subframe offset
    nCS_SRS                         % SRS cyclic shift
    Bsrs                            % SRS bandwidth
    b_hop                           % SRS hopping bandwidth
    nRRC                            % SRS frequency domain position
    kTC                             % SRS transmission comb
end

   methods
       % Class constructor, defining default values.
       % Get as input a struct defining the input variables
       function obj = UE(UE_params,Ntot,Nsub,nRX,fading,frameStructure,buffer_length)
           
           % Default values
           if ~isfield(UE_params,'turbo_iterations')
               UE_params.turbo_iterations = 8;
           end
           
           % Assign values
           obj.N_soft                       = UE_params.N_soft;
           obj.user_speed                   = UE_params.user_speed;
           obj.mode                         = UE_params.mode;
           obj.nTX                          = UE_params.nTX;
           if strcmp(fading,'BlockFading')
                obj.previous_channels       = zeros(Ntot,Nsub,10,nRX,UE_params.nTX);
           else
                obj.previous_channels        = zeros(Ntot,Nsub,1,nRX,UE_params.nTX);
           end
           obj.PMI_fb                       = UE_params.PMI_fb;
           obj.CQI_fb                       = UE_params.CQI_fb;
           obj.RI_fb                        = UE_params.RI_fb;
           
           obj.SRS_config(UE_params.Isrs,UE_params.SRS_triggerType,frameStructure);
           obj.SRS_duration                 = UE_params.SRS_duration;
           obj.nCS_SRS                      = UE_params.nCS_SRS;
           obj.Bsrs                         = UE_params.Bsrs;
           obj.b_hop                        = UE_params.b_hop;
           obj.nRRC                         = UE_params.nRRC;
           obj.kTC                          = UE_params.kTC;
           obj.MSE_buffer                   = zeros(1,buffer_length);
           

           switch UE_params.SINR_averaging.averager
               case 'EESM'
                   obj.SINR_averager        = network_elements.eesmAverager(UE_params.SINR_averaging.EESMbetas,UE_params.SINR_averaging.MCSs);
               case 'MIESM'
                   obj.SINR_averager        = network_elements.miesmAverager(UE_params.SINR_averaging.MIESMbetas,UE_params.SINR_averaging.MCSs);
               case 'HARM_MEAN'    
                   obj.SINR_averager        = network_elements.harmMeanAverager;
               otherwise
                   error('SINR averager not supported');
           end
           
           % Optional initialization 
           if isfield(UE_params,'clock')
               obj.clock                    = UE_params.clock;
           else
               obj.clock                    = 0;
           end
       end

        % see TS 136.213 section 8.2
        function SRS_config(obj,Isrs,SRS_triggerType,frameStructure)
            obj.Isrs             = Isrs;
            obj.SRS_triggerType  = SRS_triggerType;
           
            if Isrs < 0
                error('SRS index cannot be negative');
            else
                if frameStructure == 1       % FDD
                    if Isrs < 2 
                        obj.Tsrs = 2;
                        obj.Toffset = Isrs;
                    elseif Isrs < 7
                        obj.Tsrs = 5;
                        obj.Toffset = Isrs - 2;
                    elseif Isrs < 17
                        obj.Tsrs = 10;
                        obj.Toffset = Isrs - 7;
                    elseif SRS_triggerType == 0
                        if Isrs < 37
                            obj.Tsrs = 20;
                            obj.Toffset = Isrs - 17;
                        elseif Isrs < 77
                            obj.Tsrs = 40;
                            obj.Toffset = Isrs - 37;
                        elseif Isrs < 157
                            obj.Tsrs = 80;
                            obj.Toffset = Isrs - 77;
                        elseif Isrs < 317
                            obj.Tsrs = 160;
                            obj.Toffset = Isrs - 157;
                        elseif Isrs < 637
                            obj.Tsrs = 320;
                            obj.Toffset = Isrs - 317;
                        else
                            error('SRS index not supported');
                        end
                    else
                        error('SRS index not supported');
                    end
                else                        %TDD
                    if Isrs < 0
                        error('SRS index cannot be negative');
                    elseif Isrs < 10
                        obj.kSRS = 0:9;
                        obj.Tsrs = 2;
                        if Isrs == 0
                            obj.Toffset = [0 1];
                        elseif Isrs == 1
                            obj.Toffset = [0 2];
                        elseif Isrs == 2
                            obj.Toffset = [1 2];
                        elseif Isrs == 3
                            obj.Toffset = [0 3];
                        elseif Isrs == 4
                            obj.Toffset = [1 3];
                        elseif Isrs == 5
                            obj.Toffset = [0 4];
                        elseif Isrs == 6
                            obj.Toffset = [1 4];
                        elseif Isrs == 7
                            obj.Toffset = [2 3];
                        elseif Isrs == 8
                            obj.Toffset = [2 4];
                        elseif Isrs == 9
                            obj.Toffset = [3 4];
                        end
                    elseif Isrs < 15
                        obj.Tsrs = 5;
                        obj.Toffset = Isrs - 10;
                    elseif SRS_triggerType == 0
                        if Isrs < 25
                            obj.Tsrs = 10;
                            obj.Toffset = Isrs - 15;
                        elseif Isrs < 45
                            obj.Tsrs = 20;
                            obj.Toffset = Isrs - 25;
                        elseif Isrs < 85
                            obj.Tsrs = 40;
                            obj.Toffset = Isrs - 45;
                        elseif Isrs < 165
                            obj.Tsrs = 80;
                            obj.Toffset = Isrs - 85;
                        elseif Isrs < 325
                            obj.Tsrs = 160;
                            obj.Toffset = Isrs - 165;
                        elseif Isrs < 645
                            obj.Tsrs = 320;
                            obj.Toffset = Isrs - 325;
                        else
                            error('SRS index not supported');
                        end
                    else
                        error('SRS index not supported');
                    end
                end
           end
        end      
   end
end 
