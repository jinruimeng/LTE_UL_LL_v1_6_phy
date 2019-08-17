function fk = LTE_tx_code_block_concatenation(ek)
% LTE code block concatenation (TS36.212, subclause 5.1.5)
% [fk] = LTE_tx_code_block_concatenation(ek)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% input :   ek   ... code blocks to concatenate
% output:   fk   ... concatenated code blocks
%
% date of creation: 2008/08/11
% last changes:

% After checking Matlab's cell2mat, I think this implementation should be faster.
fk = ek{1};
for i=2:length(ek)
    fk=[fk ek{i}];
end