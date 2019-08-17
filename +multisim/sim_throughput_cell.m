% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function sim_throughput_cell( simulation, varargin )
%sim_throuphput_user plots each users throughput for each simulation

N_sub = length(simulation.subsimulations.keys);
N_snr = size(simulation.SNR_vec,2);

if ~isempty(simulation.labels_in_order)
    sim_labels = simulation.labels_in_order;
else
    sim_labels = simulation.subsimulations.keys;
end


throughput_coded = cell(N_sub,1);
throughput_uncoded = cell(N_sub,1);

throughput_coded_confidence = cell(N_sub,1);
throughput_uncoded_confidence = cell(N_sub,1);

plot_uncoded = true;

for k = 1:length(varargin)
    if strcmp(varargin{k}, 'no_uncoded');
        plot_uncoded = false;
    end
end

for n = 1:N_sub
    
    cell_spec = simulation.subsimulations(sim_labels{n}).simulation_results.cell_specific;
    N_BS = length(cell_spec);
    
    throughput_coded{n}.data = zeros(N_BS, N_snr);
    throughput_uncoded{n}.data = zeros(N_BS, N_snr);

    throughput_coded_confidence{n} = struct;
    throughput_coded_confidence{n}.lower = zeros(N_BS, N_snr);
    throughput_coded_confidence{n}.upper = zeros(N_BS, N_snr);
    if plot_uncoded
        throughput_uncoded_confidence{n} = struct;
        throughput_uncoded_confidence{n}.lower = zeros(N_BS, N_snr);
        throughput_uncoded_confidence{n}.upper = zeros(N_BS, N_snr);
    end
    

    for bb = 1:N_BS
        if simulation.subsimulations(sim_labels{n}).simulation_results.conf_interval_probability == 0 || isempty(fieldnames(cell_spec(bb).confidence))
            coded_tmp = cell_spec(bb).throughput_coded;
            throughput_coded{n}.data(bb, :) = mean(sum(coded_tmp,3),1) / 1e3; % tp in Mbit/s
            if plot_uncoded
                uncoded_tmp = cell_spec(bb).throughput_uncoded;
                throughput_uncoded{n}.data(bb, :) = mean(sum(uncoded_tmp,3),1) / 1e3; % tp in Mbit/s
            end
        else
            throughput_coded{n}.data(bb, :) = cell_spec(bb).confidence.throughput_coded(1,:) / 1e3;
            throughput_coded_confidence{n}.lower(bb, :) = cell_spec(bb).confidence.throughput_coded(2,:) / 1e3;
            throughput_coded_confidence{n}.upper(bb, :) = cell_spec(bb).confidence.throughput_coded(3,:) / 1e3;
            if plot_uncoded
                throughput_uncoded{n}.data(bb, :) = cell_spec(bb).confidence.throughput_uncoded(1,:) / 1e3;
                throughput_uncoded_confidence{n}.lower(bb, :) = cell_spec(bb).confidence.throughput_uncoded(2,:) / 1e3;
                throughput_uncoded_confidence{n}.upper(bb, :) = cell_spec(bb).confidence.throughput_uncoded(3,:) / 1e3;
            end
        end
    end
end

fig = figure('Name', ['N_sub = ', num2str(simulation.N_subframes)]);
cmp = colormap(jet(N_sub));
markers = plots.getCustomMarkers();
legend_entries = cell(1, N_sub);

for n = 1:N_sub
    % plot curves 
    simulation_results = simulation.subsimulations(sim_labels{n}).simulation_results;
    LTE_params = simulation.subsimulations(sim_labels{n}).LTE_params;
    cell_spec = simulation.subsimulations(sim_labels{n}).simulation_results.cell_specific;
    N_BS = length(cell_spec);
    
    for bb = 1:N_BS
        marker = plots.getMarker(bb, markers);
        show_in_legend = 'off';
        % add legend entry
        if bb == 1
            legend_entries{n} = sim_labels{n};
            show_in_legend = 'on';
        end

        plot( simulation_results.SNR_vector, throughput_coded{n}.data(bb, :), ...
            'LineStyle', '-', 'Marker', marker, 'Color', cmp(n,:), 'HandleVisibility',show_in_legend );
        hold on;
        if plot_uncoded
            plot( simulation_results.SNR_vector, throughput_uncoded{n}.data(bb, :), ...
                'LineStyle', '--', 'Marker', marker, 'Color', cmp(n,:), 'HandleVisibility','off');
        end
        

        if simulation_results.conf_interval_probability > 0 % confidence intervals computed
            for s = 1:N_snr
                plot( simulation_results.SNR_vector(s) * ones(1,2), ...
                    [throughput_coded_confidence{n}.lower(bb, s), throughput_coded_confidence{n}.upper(bb, s)], ...
                     'LineStyle', '-', 'Color', cmp(n,:), 'HandleVisibility','off');
                if plot_uncoded
                    plot( simulation_results.SNR_vector(s) * ones(1,2), ...
                    [throughput_uncoded_confidence{n}.lower(bb, s), throughput_uncoded_confidence{n}.upper(bb, s)], ...
                     'LineStyle', '--', 'Color', cmp(n,:), 'HandleVisibility','off');
                end
            end
        end

    end
end

% hack to add an extra legend entry for uncoded
if plot_uncoded
    plot(NaN, NaN, 'LineStyle', '--', 'Color', ones(1,3)*0.5);
end
hold off;

if plot_uncoded
    legend_entries{N_sub+1} = sprintf('uncoded');
end
legend(legend_entries, 'Location', 'best');



xlabel('SNR [dB]');
ylabel('throughput [Mbit/s]');
title('throughput per cell');
grid on;

if LTE_params.save_plots 
    print('report/throughput_cell','-dpng')
end

end
