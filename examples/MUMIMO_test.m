

doc_mumimo1 = multisim('examples/doc_mumimo');

multisim.sim_throughput_cell( doc_mumimo1 , 'no_uncoded');
title('');
ylabel('cell throughput [Mbit/s]');
xlabel('SNR [dB]');

set(gcf, 'PaperPosition', [-0.2 -0.3 6 5]); 
set(gcf, 'PaperSize', [5.5 4.6]); 

saveas(gcf, 'MUMIMO_test', 'pdf')