%%

clear; clc;

load('./examples/plot_much_tti.mat')

% plot mimo for documentation...%%

%load

% plot to file 1 - 4

% SISO
multisim.plot_simulation(simulation, 'subsimulations', {'SISO'}, 'mode', 'tikz', 'tikzoutput', './documentation/graphics/siso.tex');

% SISO, CLSM 2x2
multisim.plot_simulation(simulation, 'subsimulations', {'SISO', 'CLSM 2x2'}, 'mode', 'tikz', 'tikzoutput', './documentation/graphics/mimo_2x2.tex');

% SISO, CLSM 4x4
multisim.plot_simulation(simulation, 'subsimulations', {'SISO', 'CLSM 4x4'}, 'mode', 'tikz', 'tikzoutput', './documentation/graphics/mimo_4x4.tex');

% SISO, CLSM 2x2, CLSM 2x4, CLSM 4x4
multisim.plot_simulation(simulation, 'subsimulations', {'SISO', 'CLSM 2x2', 'CLSM 2x4', 'CLSM 4x4'}, 'mode', 'tikz', 'tikzoutput', './documentation/graphics/mimo_all.tex');