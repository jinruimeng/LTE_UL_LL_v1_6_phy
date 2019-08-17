% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function sim_throughput_user( simulation, varargin )
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
    
    UE_spec = simulation.subsimulations(sim_labels{n}).simulation_results.UE_specific;
    N_UE = length(UE_spec);
    
    throughput_coded{n}.data = zeros(N_UE, N_snr);
    throughput_uncoded{n}.data = zeros(N_UE, N_snr);

    throughput_coded_confidence{n} = struct;
    throughput_coded_confidence{n}.lower = zeros(N_UE, N_snr);
    throughput_coded_confidence{n}.upper = zeros(N_UE, N_snr);
    if plot_uncoded
        throughput_uncoded_confidence{n} = struct;
        throughput_uncoded_confidence{n}.lower = zeros(N_UE, N_snr);
        throughput_uncoded_confidence{n}.upper = zeros(N_UE, N_snr);
    end
    

    for u = 1:N_UE
        if simulation.subsimulations(sim_labels{n}).simulation_results.conf_interval_probability == 0 || isempty(fieldnames(UE_spec(u).confidence))
            coded_tmp = UE_spec(u).throughput_coded;
            throughput_coded{n}.data(u, :) = mean(sum(coded_tmp,3),1) / 1e3; % tp in Mbit/s
            if plot_uncoded
                uncoded_tmp = UE_spec(u).throughput_uncoded;
                throughput_uncoded{n}.data(u, :) = mean(sum(uncoded_tmp,3),1) / 1e3; % tp in Mbit/s
            end
        else
            throughput_coded{n}.data(u, :) = UE_spec(u).confidence.throughput_coded(1,:) / 1e3;
            throughput_coded_confidence{n}.lower(u, :) = UE_spec(u).confidence.throughput_coded(2,:) / 1e3;
            throughput_coded_confidence{n}.upper(u, :) = UE_spec(u).confidence.throughput_coded(3,:) / 1e3;
            if plot_uncoded
                throughput_uncoded{n}.data(u, :) = UE_spec(u).confidence.throughput_uncoded(1,:) / 1e3;
                throughput_uncoded_confidence{n}.lower(u, :) = UE_spec(u).confidence.throughput_uncoded(2,:) / 1e3;
                throughput_uncoded_confidence{n}.upper(u, :) = UE_spec(u).confidence.throughput_uncoded(3,:) / 1e3;
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
    UE_spec = simulation.subsimulations(sim_labels{n}).simulation_results.UE_specific;
    N_UE = length(UE_spec);
    
    for u = 1:N_UE
        marker = plots.getMarker(u, markers);
        show_in_legend = 'off';
        % add legend entry
        if u == 1
            legend_entries{n} = sim_labels{n};
            show_in_legend = 'on';
        end

        plot( simulation_results.SNR_vector, throughput_coded{n}.data(u, :), ...
            'LineStyle', '-', 'Marker', marker, 'Color', cmp(n,:), 'HandleVisibility',show_in_legend );
        hold on;
        if plot_uncoded
            plot( simulation_results.SNR_vector, throughput_uncoded{n}.data(u, :), ...
                'LineStyle', '--', 'Marker', marker, 'Color', cmp(n,:), 'HandleVisibility','off');
        end
        

        if simulation_results.conf_interval_probability > 0 % confidence intervals computed
            for s = 1:N_snr
                plot( simulation_results.SNR_vector(s) * ones(1,2), ...
                    [throughput_coded_confidence{n}.lower(u, s), throughput_coded_confidence{n}.upper(u, s)], ...
                     'LineStyle', '-', 'Color', cmp(n,:), 'HandleVisibility','off');
                if plot_uncoded
                    plot( simulation_results.SNR_vector(s) * ones(1,2), ...
                    [throughput_uncoded_confidence{n}.lower(u, s), throughput_uncoded_confidence{n}.upper(u, s)], ...
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
title('throughput per user');
grid on;

if LTE_params.save_plots 
    print('report/throughput_user','-dpng')
end

end
