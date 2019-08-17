

doc_mimo = multisim('examples/doc_mimo');

multisim.sim_throughput_user( doc_mimo , 'no_uncoded');
title('');

set(gcf, 'PaperPosition', [-0.2 -0.3 6 5]); 
set(gcf, 'PaperSize', [5.5 4.6]); 

saveas(gcf, 'MIMO_test', 'pdf')