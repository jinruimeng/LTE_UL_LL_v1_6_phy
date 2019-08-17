%% plot simply tries to plot simulations creaded with multisim

if exist('simulation', 'var')
    multisim.plot_simulation2(simulation);
else
    files = dir('multisim_results/*.mat');
    [~, ind] = sort([files.datenum],'descend');
    sim = load(['multisim_results/',files(ind(1)).name]);

    multisim.plot_simulation2(sim.simulation);
end       


