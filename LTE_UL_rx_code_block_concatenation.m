function ek = LTE_UL_rx_code_block_concatenation(fk,BS_signaling,stream_index)
% LTE code block concatenation (TS36.212, subclause 5.1.5)
% [ek] = LTE_rx_code_block_concatenation(fk,BS_signaling)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input:    fk   ... concatenated code blocks
% output:   ek   ... code blocks to concatenate
%
% date of creation: 2008/08/11
% last changes: 2008/09/15  Bosanska    added input [1 x 1]struct BS_signaling (corresponds to changed BS_output struct)


for i=1:BS_signaling.TB_segmentation_UL(stream_index).C
    if i==1
        index_begin = 1;
        index_end = BS_signaling.turbo_rate_matcher_UL(stream_index).ek_sizes(1);
    else
        index_begin = index_end+1;
        index_end = index_begin+BS_signaling.turbo_rate_matcher_UL(stream_index).ek_sizes(i)-1;
    end
    
    ek{i} = fk(index_begin:index_end);
end