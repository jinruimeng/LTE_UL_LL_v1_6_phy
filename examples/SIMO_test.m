%%

doc_simo = multisim('examples/doc_simo');

multisim.sim_throughput_user( doc_simo , 'no_uncoded');
title('');

set(gcf, 'PaperPosition', [-0.2 -0.3 6 5]); 
set(gcf, 'PaperSize', [5.5 4.6]); 

saveas(gcf, 'MIMO_test', 'pdf')
