%%
%cd ..
%%

awgn_test1 = multisim('examples/doc_awgn1');


multisim.sim_throughput_user( awgn_test1 , 'no_uncoded');
title('AWGN throughput 1.4 MHz, 2000 subframes');


set(gcf, 'PaperPosition', [0 -0.5 7 5]);
set(gcf, 'PaperSize', [6 4.5]); 


lh=findall(gcf,'tag','legend');
set(lh,'location','northeastoutside');


print('AWGN_test1_throughput', '-dpdf', '-r300');

%%
multisim.sim_bler_user( awgn_test1 );
title('AWGN BLER 1.4 MHz, 2000 subframes');
ylim([1e-3,1]);

set(gcf, 'PaperPosition', [0 0 9 7]); 
set(gcf, 'PaperSize', [9 7]); 

pause(1);
print('AWGN_test1_bler', '-dpdf', '-r300');


