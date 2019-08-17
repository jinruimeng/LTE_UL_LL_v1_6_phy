% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function sim_bler_user( simulation )
%sim_bler_user plots each users bler for each simulation

N_sub = length(simulation.subsimulations.keys);
N_snr = size(simulation.SNR_vec,2);

if ~isempty(simulation.labels_in_order)
    sim_labels = simulation.labels_in_order;
else
    sim_labels = simulation.subsimulations.keys;
end

bler = cell(N_sub,1);


bler_confidence = cell(N_sub,1);

for n = 1:N_sub
    
    UE_spec = simulation.subsimulations(sim_labels{n}).simulation_results.UE_specific;
    N_UE = length(UE_spec);
    
    bler{n}.data = zeros(N_UE, N_snr);

    bler_confidence{n} = struct;
    bler_confidence{n}.lower = zeros(N_UE, N_snr);
    bler_confidence{n}.upper = zeros(N_UE, N_snr);
    
    

    for u = 1:N_UE
        if simulation.subsimulations(sim_labels{n}).simulation_results.conf_interval_probability == 0 || isempty(fieldnames(UE_spec(u).confidence))
            bler{n}.data(u, :) = UE_spec(u).BLER_overall;

        else
            bler{n}.data(u, :) = UE_spec(u).confidence.BLER_coded(1,:);
            bler_confidence{n}.lower(u, :) = UE_spec(u).confidence.BLER_coded(2,:);
            bler_confidence{n}.upper(u, :) = UE_spec(u).confidence.BLER_coded(3,:);
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

        semilogy( simulation_results.SNR_vector, bler{n}.data(u, :), ...
            'LineStyle', '-', 'Marker', marker, 'Color', cmp(n,:), 'HandleVisibility',show_in_legend );
        hold on;

        

        if simulation_results.conf_interval_probability > 0 % confidence intervals computed
            for s = 1:N_snr
                semilogy( simulation_results.SNR_vector(s) * ones(1,2), ...
                    [bler_confidence{n}.lower(u, s), bler_confidence{n}.upper(u, s)], ...
                     'LineStyle', '-', 'Color', cmp(n,:), 'HandleVisibility','off');
            end
        end

    end
end

% hack to add an extra legend entry for uncoded
hold off;

legend(legend_entries, 'Location', 'best');



xlabel('SNR [dB]');
ylabel('BLER');
title('BLER per user');
grid on;

% if LTE_params.save_plots 
%     print('report/ber_user','-dpng')
% end

end
