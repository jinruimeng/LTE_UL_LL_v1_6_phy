function f=LTE_UL_tx_ULSCH_encode(LTE_params,a,UE_signaling,UE_genie,UE,stream_index)

UE_genie.data_bits{stream_index} = a;

%% Transport block CRC, Code block Segmentation
b = LTE_tx_append_crc(a,'24a');
UE_signaling.TB_size_UL(stream_index) = length(b);

c = LTE_UL_tx_code_block_segmentation(LTE_params,b,UE_signaling,stream_index);
UE_genie.bits_to_turboencode{stream_index} = c;

%% Channel Coding, Rate Matching, Code Block concatenation
d = LTE_UL_tx_turbo_encode(LTE_params,c,UE_signaling,stream_index);

e = LTE_UL_tx_turbo_rate_matcher2(LTE_params,d,UE_signaling,UE,stream_index);

f = LTE_tx_code_block_concatenation(e);
UE_genie.sent_bits{stream_index} = f;