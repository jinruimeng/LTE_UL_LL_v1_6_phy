function output = LTE_common_turbo_rate_matcher_circular_buffer(input,interleave_flag)
% Circular buffer part of the LTE Turbo Code Rate Matcher, as of TS 36.212, Section 5.1.4.1.
% [output] = LTE_common_turbo_rate_matcher_circular_buffer(input,interleave_flag)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input :   input_bits        ... if interleave_flag==1: Output of the
%                                  sub-block interleaver (3xN matrix)
%                                 if interleave_flag!=1: row vector
%           interleave_flag   ... 1 to interleave, 0 to deinterleave
% output:   output_bits       ... if interleave_flag==1: row vector (output
%                                  of the buffer)
%                                 if interleave_flag!=1: 3xN matrix that
%                                  can be inputted to the sub-block
%                                  deinterleaver
%
% date of creation: 2008/08/11
% last changes:

% When interleaving we interleave uint8s (bits). When deinterleaving,
% what we deinterleave are soft bits (prior to turbo decoding).
if interleave_flag==1
    interleaving_function = @LTE_common_bit_interleaver;
    K_pi = size(input,2);    % input is 3xK_pi matrix
    data_to_process = [ input(1,:) input(2,:) input(3,:) ];
else
    interleaving_function = @LTE_common_soft_bit_interleaver;
    K_pi = size(input,2)/3;  % input is row vector
    data_to_process = input;
end

circular_buffer_mapping = LTE_common_turbo_rate_matching_circular_buffer_mapping(K_pi);
processed_data = interleaving_function(data_to_process,circular_buffer_mapping,interleave_flag);

if interleave_flag==1
    output = processed_data;
else
    output(1,:) = processed_data(1:K_pi);
    output(2,:) = processed_data((K_pi+1):(2*K_pi));
    output(3,:) = processed_data((2*K_pi+1):(3*K_pi));
end

