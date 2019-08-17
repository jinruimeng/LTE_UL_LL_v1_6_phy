function [output] = LTE_UL_tx_code_block_segmentation(LTE_params,transport_block,UE_signaling,stream_index)
% Performs code block segmentation and CRC attachment according to TS 36.212, Section 5.1.2
% [output signaling_info] = LTE_tx_code_block_segmentation(transport_block)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input:    transport_block     ... logical vector representing a transport
%                                   block
%
% output:   code_blocks         ... cell array containing the several code
%                                   blocks. code_block{1} is the first one
%           signaling_info      ... info to be signalled (numbe rof filler
%                                   bits, etc)
%
% date of creation: 2008/08/11
% last changes:

%% Segmentation logic
% Number of bits of the Transport Block
B = length(transport_block);
% Given max code block size
Z =6144;

%transport_block = uint8(transport_block);

if B<40
    % add filler bits and return. They are <NULL> bits (technically a
    % uint(3), but since the encoder later replaces this with 0s, we will just use 0s
    F = 40-B;
    output{1}(1:F) = false;
    output{1} = [output{1} transport_block];
    C = 1;
    signaling_info.K_plus  = 40;
    signaling_info.K_minus = 0;
    signaling_info.F = F;
    signaling_info.C = C;
    
    UE_signaling.TB_segmentation(stream_index).F = signaling_info.F;
    UE_signaling.TB_segmentation(stream_index).C = signaling_info.C;
    UE_signaling.TB_segmentation(stream_index).K_plus  = signaling_info.K_plus;
    UE_signaling.TB_segmentation(stream_index).K_minus = signaling_info.K_minus;
    UE_signaling.TB_segmentation(stream_index).Z = Z;
    UE_signaling.TB_segmentation(stream_index).B = B;
    
    return
end

if  B<=Z
	L = 0;
	C = 1;
    B_p = B;
else
	L = 24;
	C = ceil(B/(Z-L));
    B_p = B + C*L;
end

%% Get K+ and K-
% Calculate the index in the table from the Block size
avg_code_block_size = (B_p/C);
index_number = avg_code_block_size/8 - 4;

if index_number > 60 && index_number <= 124
    real_index = avg_code_block_size/16 + 28;
elseif index_number > 124 && index_number <= 252
    real_index = avg_code_block_size/32 + 60;
elseif index_number > 252
    real_index = avg_code_block_size/64 + 92;
else
    real_index = index_number;
end

% First segmentation size:   = minimum K in table 5.1.3-3 such that   C K \geq B'
% Second segmentation size:  = maximum K in table 5.1.3-3 such that K_- < K_+
K_plus = LTE_params.turbo_interleaver_table(ceil(real_index),1);
if C==1
    K_minus = 0;
else
    K_minus = LTE_params.turbo_interleaver_table(ceil(real_index)-1,1);
end

%% Determine the segmentation
% Calculate number of blocks of each size
if C==1
    C_plus = 1;
    K_minus = 0;
    C_minus = 0;
else
    diff_K = K_plus - K_minus;
    
    C_minus = floor((C*K_plus-B_p)/diff_K);
    C_plus = C - C_minus;
end

% Calculate filler bits
F = C_plus*K_plus + C_minus*K_minus - B_p;

%% Generate ouput (cell array of code blocks)
% Add filler bits. They are <NULL> bits (technically a
% uint(3), but since the encoder later replaces this with 0s, we will just use 0s
output{1}(1:F) = false;
% Add the rest of bits
k = F+1; %output position index
s = 1;   %input position index
for r=1:C
    if r<=C_minus
        K_r = K_minus;
    else
        K_r = K_plus;
    end
    bits_to_write = K_r-L;
    if r==1 %The first code block may have filler bits
        bits_to_write = bits_to_write-F;
    end
    output{r}(k:(k+bits_to_write-1)) = transport_block(s:(s+bits_to_write-1));
    %Update position of indexes
    s = s + bits_to_write;
    k = 1;
end

% Add CRC to the code blocks when necessary
if C>1
    for i=1:C
        output{i} = LTE_tx_append_crc(output{i},'24b');
    end
end
% UE_signaling.TB_segmentation(stream_index).F
% Generate the data that will be signalled
UE_signaling.TB_segmentation_UL(stream_index).F = F;
UE_signaling.TB_segmentation_UL(stream_index).C = C;
UE_signaling.TB_segmentation_UL(stream_index).K_plus  = K_plus;
UE_signaling.TB_segmentation_UL(stream_index).K_minus = K_minus;
UE_signaling.TB_segmentation_UL(stream_index).Z = Z;
UE_signaling.TB_segmentation_UL(stream_index).B = B;
for i=1:length(output)
    UE_signaling.TB_segmentation_UL(stream_index).CB_sizes(i) = length(output{i});
end