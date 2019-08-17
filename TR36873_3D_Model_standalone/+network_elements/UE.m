classdef UE < handle
    % This object represents the UE with properties
    % id                                        ... user id
    % pos                                       ... user position [x y]
    % pos_pixel                                 ... pixel position [x y]
    % attached_site                             ... site to where this user is attached
    % attached_sector_idx                       ... sector index to which the UE is attached
    % attached_eNodeB                           ... eNodeB to which the UE is attached to
    % channels                                  ... a  list of teh desired channel (small scale parameters)
    % nRX                                       ... number of antenna ports
    % all_large_scale_params                    ... a vector of seven large scale parameters listed as matrix: 
    %                                               each row determines the large scale parameters for serving eNodeB (1st row) 
    %                                               and interfering eNodeBs (2nd row -to the end)
    %                                               a vector contains the large scale parameters in the following order [SF, K-factor, DS, ASD, ASA, ZSD, ZSA]
    % all_ZOD_params                            ... offset parameters for Zenith Spread of Departure (ZSD) Table 7.3-7 and 7.3-8 listed as matrix
    %                                               each row determines the offset parameters for serving eNodeB (1st row) 
    %                                               and interfering eNodeBs (2nd row -to the end)                   
    %                                               a vector contains
    %                                               [mean, std, offset_value] of the ZSD
    % rx_height                                 ... user height in [m]
    % is_LOS                                    ... boolean type: 1-LOS; 0-NLOS
    % is_indoor                                 ... boolean type: 1-indoor; 0-outdoor
    % dist_indoor                               ... distance indoor in [m]
    % direction                                 ... user pointing direction in [°]
    % H_0_channel_trace                         ... channel impulse response of the desired link -the output from TR 36.873 eq(7.21/7.26)
    % sampled_channel_H_0                       ... sampled channel impulse response of the desired link
    % H_0_final                                 ... the channel transfer function of the desired link (after FFT of the sampled version)
    % H_i_channel_trace                         ... channel impulse response of the interfering links -the output from TR 36.873 eq(7.21/7.26)
    % sampled_channel_H_i                       ... sampled channel impulse response of the interfering links
    % H_i_after_fft                             ... the channel transfer function of the interfering links (after FFT of the sampled version)
    % H_i_full_final                            ... the channel tranfer function of the interfering links after applying the pathloss and shadow fading
    % TTI_of_smallscale_fading_recalculation    ...
    % TTI_of_smallscale_fading_recalculation_i  ...
    % recalculate_3D_smallscale_fading          ...
    % recalculate_3D_smallscale_fading_i        ...
    %
    % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016
    
    properties
        id
        pos
        pos_pixel
        attached_site         
        attached_sector_idx   
        attached_eNodeB       
        channels
        nRX
        all_large_scale_params
        all_ZOD_params 
        rx_height
        is_LOS
        is_indoor 
        dist_indoor 
        direction
        H_0_channel_trace 
        sampled_channel_H_0
         H_0_final
        H_i_channel_trace
        sampled_channel_H_i 
        H_i_after_fft 
        H_i_full_final
        TTI_of_smallscale_fading_recalculation
        TTI_of_smallscale_fading_recalculation_i
        recalculate_3D_smallscale_fading
        recalculate_3D_smallscale_fading_i
        azimuth_angles_of_departure_all 
        azimuth_angles_of_arrival_all
      zenith_angles_of_departure_all  
    zenith_angles_of_arrival_all
    cluster_power_per_ray
    
        

    end
      
end

