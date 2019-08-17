stream_index = 1;
a = tx_data_bits;
UE_output(1,uplink_user).a_tx = a;
BS_output.genie(1,uplink_user).data_bits{1} = tx_data_bits;

%% Transport block CRC, Code block Segmentation
b = LTE_tx_append_crc(a,'24a');
UE_signaling.TB_size(stream_index) = length(b);

%%
BS_output.UE_signaling(uu).turbo_rate_matcher(cw).G = N_coded_bits; % How many bits the we are allowed to transmit. Decided by the scheduler(uu)
BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 1;

BS_output.UE_signaling(uu).turbo_rate_matcher(cw).rv_idx = BS.UE_specific(uu).current_HARQ_process(cw).rv_idx;

UE_signaling = BS_output.UE_signaling;

c = LTE_UL_tx_code_block_segmentation(LTE_params,b,BS_output.UE_signaling,stream_index,uu);
UE_output(1,uplink_user).c_tx = c;
UE_output(1,uplink_user).UE_genie.bits_to_turboencode = c;



%% Channel Coding, Rate Matching, Code Block concatenation
d = LTE_UL_tx_turbo_encode(LTE_params,c,BS_output.UE_signaling,stream_index,uu);
UE_output(1,uplink_user).d = d{1};

% rate matcher fixed, Prokopec
e = LTE_tx_turbo_rate_matcher(LTE_params,d,UE_signaling,UE,uu,stream_index);


UE_output(1,uplink_user).e = e{1};
f = LTE_tx_code_block_concatenation(e);
UE_output(1,uplink_user).f = f;
UE_genie.sent_bits{stream_index} = f;