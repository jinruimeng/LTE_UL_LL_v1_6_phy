function [ y_perm Nd R K_pi ] = LTE_common_subblock_interleaver(LTE_params,input,interleave_flag,varargin)
% LTE Turbo Sub-Block interleaver/deinterleaver, as of TS 36.212, Section 5.1.4.1.
% [y_perm Nd R K_pi] = LTE_common_subblock_interleaver(input,interleave_flag,varargin)
% Author: Josep Colom Ikuno, josep.colom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input :   input             ... dk, as of TS 36.212, Section 5.1.4.1 if you
%                                 are interleaving (logical array). The interleaved 
%                                 sequence (LLRs, soft bits) in the other
%                                 case. This is a matrix of 3xN bits/soft
%                                 bits.
%           interleave_flag   ... 1 to interleave, 0 to deinterleave
%           padding_bits      ... when deinterleaving, the number of
%                                 padding dummy bits is required to delete them.
% output:   output            ... interleaved codeword (logicals) when interleave==1
%                             ... deinterleaved codeword (LLRs) when
%                                 interleave==0
%           Nd                ... number of padding bits used for
%                                 interleaving. When deinterleaving, its 
%                                 the same value as the padding_bits input.
%                                 Not used when deinterleaving (set to 0)
%           R                 ... number of rows used in the sub-block
%                                 interleaver. Not used when deinterleaving
%                                 (set to 0)
%           K_pi              ... length of the resulting interleaved
%                                 sequence (actually R*32). Not used when
%                                 deinterleaving (set to 0)
%
% date of creation: 2008/08/11

% Interleaving
if interleave_flag==1
    [ y_perm(1,:) Nd R K_pi ] = LTE_common_subblock_interleaver_block(LTE_params,input(1,:),interleave_flag,0);
      y_perm(2,:)             = LTE_common_subblock_interleaver_block(LTE_params,input(2,:),interleave_flag,1);
      y_perm(3,:)             = LTE_common_subblock_interleaver_block(LTE_params,input(3,:),interleave_flag,2);

% Deinterleaving
else
    Nd   = 0;
    R    = 0;
    K_pi = 0;
    y_perm(1,:) = LTE_common_subblock_interleaver_block(LTE_params,input(1,:),interleave_flag,0,varargin{1});
    y_perm(2,:) = LTE_common_subblock_interleaver_block(LTE_params,input(2,:),interleave_flag,1,varargin{1});
    y_perm(3,:) = LTE_common_subblock_interleaver_block(LTE_params,input(3,:),interleave_flag,2,varargin{1});
end

function [ y_perm Nd R K_pi ] = LTE_common_subblock_interleaver_block(LTE_params,input,interleave_flag,dk_index,varargin)
% LTE Turbo Sub-Block interleaver/deinterleaver, as of TS 36.212, Section 5.1.4.1.
% (c) Josep Colom Ikuno, ITC
% josep.colom@nt.tuwien.ac.at
% www.nt.tuwien.ac.at
%
% input :   input             ... dk(i), as of TS 36.212, Section 5.1.4.1 if you
%                                 are interleaving (logical array). The interleaved 
%                                 sequence (LLRs, soft bits) in the other case.
%           interleave_flag   ... 1 to interleave, 0 to deinterleave
%           dk_index          ... the interleaving is different for
%                                 {dk(0),dk(1)} than for {dk(2)}. This flag
%                                 signals which interleaving should be
%                                 applied. 0 for dk(0), 1 for dk(1) and 2
%                                 for dk(2).
%           padding_bits      ... when deinterleaving, the number of
%                                 padding bits is required to delete them.
%           R_tc              ... number of rows use din the sub-block
%                                 interleaver
%           K_pi              ... sub-block interleaver size
% output:   output            ... interleaved codeword (logicals) when interleave==1
%                             ... deinterleaved codeword (LLRs) when interleave==0

% When interleaving we interleave logicals (bits). When deinterleaving,
% what we deinterleave are soft bits (prior to turbo decoding).
if interleave_flag==1
    interleaving_function = @LTE_common_bit_interleaver;
else
    interleaving_function = @LTE_common_soft_bit_interleaver;
end

D = length(input); % Number of bits
C = 32;            % Number of columns of the permutation matrix
R = ceil(D/32);    % Number of rows of the permutation matrix
Nd = R*C-D;        % Number of padding bits

% Needed to treat the Nd <NULL> padding bits. In our convention:
%   0 --> 0 bit
%   1 --> 1 bit
%   3 --> <NULL> bit

% If we are interleaving, the padding bits are inputted at the beginning of
% the sequence.
if interleave_flag==1
    data_to_interleave = [3*ones(1,Nd,'uint8') uint8(input)];
else
    Nd_to_delete = varargin{1};
    data_to_interleave = input;
end

% For the first and second streams {dk(0),dk{1}) use the column
% permutations. for dk(2), perform a simple plain interleaving
% See TS 36.212, Section 5.1.4.1.1 for more details.
if dk_index==0 || dk_index==1
    if interleave_flag==1
        % Put the data row by row into the permutation matrix
        data_to_interleave = reshape(data_to_interleave,32,[])';
        % Permute the columns
        y_perm             = data_to_interleave(:,LTE_params.sub_block_interleaver_permutation_pattern+1);
        % Read the permuted matrix column-wise. This is the output
        y_perm             = reshape(y_perm,1,[]);
    else
        % Put the data column by column into the permutation matrix. Data
        % is multiple of 32, so no need to check that
        y = reshape(data_to_interleave',[],32);
        y_perm = zeros(size(y));
        % De-interleave columns
        y_perm(:,LTE_params.sub_block_interleaver_permutation_pattern+1) = y;
        y_perm = reshape(y_perm',1,[]);
    end
    K_pi = D + Nd; % Size of the interleaver: data size + padding
% Interleaving for dk(2)
else
    K_pi = D + Nd; % Size of the interleaver: data size + padding
    k = 0:K_pi-1;
    P = LTE_params.sub_block_interleaver_permutation_pattern(floor(k/R)+1); % From TS 36.212, Table 5.1.4-1
    k_mod_R = mod(k,R);
    interleaving_map = mod(P+C*k_mod_R+1,K_pi);
    y_perm = interleaving_function(data_to_interleave,interleaving_map,interleave_flag);
end

% Take out the padding bits when deinterleaving
if interleave_flag~=1
    y_perm = y_perm(1+Nd_to_delete:end);
end


