% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function sim_ber_user( simulation )
%sim_ber_user plots each users ber for each simulation

N_sub = length(simulation.subsimulations.keys);
N_snr = size(simulation.SNR_vec,2);

if ~isempty(simulation.labels_in_order)
    sim_labels = simulation.labels_in_order;
else
    sim_labels = simulation.subsimulations.keys;
end

ber_coded = cell(N_sub,1);
ber_uncoded = cell(N_sub,1);

ber_coded_confidence = cell(N_sub,1);
ber_uncoded_confidence = cell(N_sub,1);

for n = 1:N_sub
    
    UE_spec = simulation.subsimulations(sim_labels{n}).simulation_results.UE_specific;
    N_UE = length(UE_spec);
    
    ber_coded{n}.data = zeros(N_UE, N_snr);
    ber_uncoded{n}.data = zeros(N_UE, N_snr);

    ber_coded_confidence{n} = struct;
    ber_coded_confidence{n}.lower = zeros(N_UE, N_snr);
    ber_coded_confidence{n}.upper = zeros(N_UE, N_snr);

    ber_uncoded_confidence{n} = struct;
    ber_uncoded_confidence{n}.lower = zeros(N_UE, N_snr);
    ber_uncoded_confidence{n}.upper = zeros(N_UE, N_snr);
    
    

    for u = 1:N_UE
        if simulation.subsimulations(sim_labels{n}).simulation_results.conf_interval_probability == 0 || isempty(fieldnames(UE_spec(u).confidence))
            coded_tmp = UE_spec(u).BER_coded;
            ber_coded{n}.data(u, :) = mean(sum(coded_tmp,3),1);

            uncoded_tmp = UE_spec(u).BER_uncoded;
            ber_uncoded{n}.data(u, :) = mean(sum(uncoded_tmp,3),1);
        else
            ber_coded{n}.data(u, :) = UE_spec(u).confidence.BER_coded(1,:);
            ber_coded_confidence{n}.lower(u, :) = UE_spec(u).confidence.BER_coded(2,:);
            ber_coded_confidence{n}.upper(u, :) = UE_spec(u).confidence.BER_coded(3,:);

            ber_uncoded{n}.data(u, :) = UE_spec(u).confidence.BER_uncoded(1,:);
            ber_uncoded_confidence{n}.lower(u, :) = UE_spec(u).confidence.BER_uncoded(2,:);
            ber_uncoded_confidence{n}.upper(u, :) = UE_spec(u).confidence.BER_uncoded(3,:);
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

        semilogy( simulation_results.SNR_vector, ber_coded{n}.data(u, :), ...
            'LineStyle', '-', 'Marker', marker, 'Color', cmp(n,:), 'HandleVisibility',show_in_legend );
        hold on;
        semilogy( simulation_results.SNR_vector, ber_uncoded{n}.data(u, :), ...
            'LineStyle', '--', 'Marker', marker, 'Color', cmp(n,:), 'HandleVisibility','off');

        

        if simulation_results.conf_interval_probability > 0 % confidence intervals computed
            for s = 1:N_snr
                semilogy( simulation_results.SNR_vector(s) * ones(1,2), ...
                    [ber_coded_confidence{n}.lower(u, s), ber_coded_confidence{n}.upper(u, s)], ...
                     'LineStyle', '-', 'Color', cmp(n,:), 'HandleVisibility','off');
                semilogy( simulation_results.SNR_vector(s) * ones(1,2), ...
                    [ber_uncoded_confidence{n}.lower(u, s), ber_uncoded_confidence{n}.upper(u, s)], ...
                     'LineStyle', '--', 'Color', cmp(n,:), 'HandleVisibility','off');
            end
        end

    end
end

% hack to add an extra legend entry for uncoded
plot(NaN, NaN, 'LineStyle', '--', 'Color', ones(1,3)*0.5);
hold off;

legend_entries{N_sub+1} = sprintf('uncoded');
legend(legend_entries, 'Location', 'best');



xlabel('SNR [dB]');
ylabel('BER');
title('ber per user');
grid on;

% if LTE_params.save_plots 
%     print('report/ber_user','-dpng')
% end

end
