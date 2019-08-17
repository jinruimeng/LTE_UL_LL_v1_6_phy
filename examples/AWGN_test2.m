

awgn_test2 = multisim('examples/doc_awgn2');




multisim.sim_throughput_user( awgn_test2 , 'no_uncoded');

lh=findall(gcf,'tag','legend');
set(lh,'location','northeastoutside');
title('AWGN throughput 1.4 MHz, 2000 subframes');
%

set(gcf, 'PaperPosition', [0 0 9 7]); 
set(gcf, 'PaperSize', [9 7]); 


pause(1)
saveas(gcf, 'AWGN_test2_throughput', 'pdf')

%%

multisim.sim_bler_user( awgn_test2 );
title('AWGN BLER 1.4 MHz, 2000 subframes');
ylim([1e-3,1]);


set(gcf, 'PaperPosition', [0 0 9 7]); 
set(gcf, 'PaperSize', [9 7]); 

lh=findall(gcf,'tag','legend');
set(lh,'location','northeastoutside');
pause(1)
print('AWGN_test2_bler', '-dpdf', '-r300');
