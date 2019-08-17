%% config script for subsimulations

simulation.name = 'awgn doc1';

simulation.parameter_scripts = {'examples/load_params_doc'};

% global settings
simulation.DEBUG_LEVEL = 1;
simulation.cqi_i = 1;   % irrelevant if scheduler.assignment = 'dynamic'
simulation.SNR_vec = linspace(-10,20,120);  % SNR points
simulation.N_subframes = 100;  
simulation.show_plot = false;



% different subsimulations
s01 = multisim.Config('CQI 1');
s01.setCQI(1);
simulation.add(s01);

s02 = multisim.Config('CQI 2');
s02.setCQI(2);
simulation.add(s02);

s03 = multisim.Config('CQI 3');
s03.setCQI(3);
simulation.add(s03);

s04 = multisim.Config('CQI 4');
s04.setCQI(4);
simulation.add(s04);

s05 = multisim.Config('CQI 5');
s05.setCQI(5);
simulation.add(s05);

s06 = multisim.Config('CQI 6');
s06.setCQI(6);
simulation.add(s06);

s07 = multisim.Config('CQI 7');
s07.setCQI(7);
simulation.add(s07);

s08 = multisim.Config('CQI 8');
s08.setCQI(8);
simulation.add(s08);

s09 = multisim.Config('CQI 9');
s09.setCQI(9);
simulation.add(s09);

s10 = multisim.Config('CQI 10');
s10.setCQI(10);
simulation.add(s10);

s11 = multisim.Config('CQI 11');
s11.setCQI(11);
simulation.add(s11);

s12 = multisim.Config('CQI 12');
s12.setCQI(12);
simulation.add(s12);

s13 = multisim.Config('CQI 13');
s13.setCQI(13);
simulation.add(s13);

s14 = multisim.Config('CQI 14');
s14.setCQI(14);
simulation.add(s14);

s15 = multisim.Config('CQI 15');
s15.setCQI(15);
simulation.add(s15);


