% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function throughput_user( simulation_results, LTE_params )
%throuphput_user plots each users throughput

UE_spec = simulation_results.UE_specific;
N_UE = length(UE_spec);
N_snr = size(simulation_results.SNR_vector,2);

% calculate throughput per user
throughput_coded = zeros(N_UE, N_snr);
throughput_uncoded = zeros(N_UE, N_snr);

throughput_coded_confidence = struct;
throughput_coded_confidence.lower = zeros(N_UE, N_snr);
throughput_coded_confidence.upper = zeros(N_UE, N_snr);

throughput_uncoded_confidence = struct;
throughput_uncoded_confidence.lower = zeros(N_UE, N_snr);
throughput_uncoded_confidence.upper = zeros(N_UE, N_snr);

for u = 1:N_UE
    if simulation_results.conf_interval_probability == 0 || isempty(fieldnames(UE_spec(u).confidence))
        coded_tmp = UE_spec(u).throughput_coded;
        throughput_coded(u, :) = mean(sum(coded_tmp,3),1) / 1e3; % tp in Mbit/s

        uncoded_tmp = UE_spec(u).throughput_uncoded;
        throughput_uncoded(u, :) = mean(sum(uncoded_tmp,3),1) / 1e3; % tp in Mbit/s
    else
        throughput_coded(u, :) = UE_spec(u).confidence.throughput_coded(1,:) / 1e3;
        throughput_coded_confidence.lower(u, :) = UE_spec(u).confidence.throughput_coded(2,:) / 1e3;
        throughput_coded_confidence.upper(u, :) = UE_spec(u).confidence.throughput_coded(3,:) / 1e3;
        
        throughput_uncoded(u, :) = UE_spec(u).confidence.throughput_uncoded(1,:) / 1e3;
        throughput_uncoded_confidence.lower(u, :) = UE_spec(u).confidence.throughput_uncoded(2,:) / 1e3;
        throughput_uncoded_confidence.upper(u, :) = UE_spec(u).confidence.throughput_uncoded(3,:) / 1e3;
    end
end


fig = figure('Name', ['N_sub = ', num2str(LTE_params.N_subframes)]);
cmp = colormap(jet(N_UE));
markers = plots.getCustomMarkers();
legend_entries = cell(1, N_UE);

% plot curves 
for u = 1:N_UE
    marker = plots.getMarker(u, markers);
    
    plot( simulation_results.SNR_vector, throughput_coded(u, :), ...
        'LineStyle', '-', 'Marker', marker, 'Color', cmp(u,:) );
    hold on;
    plot( simulation_results.SNR_vector, throughput_uncoded(u, :), ...
        'LineStyle', '--', 'Marker', marker, 'Color', cmp(u,:), 'HandleVisibility','off');
    
    % add legend entry
    if LTE_params.nBS > 1
        b = utils.findBS( LTE_params.connection_table, u );
        legend_entries{u} = sprintf('UE %d (BS %d)', u, b);
    else
        legend_entries{u} = sprintf('UE %d', u);
    end
    
    if simulation_results.conf_interval_probability > 0 % confidence intervals computed
        for s = 1:N_snr
            plot( simulation_results.SNR_vector(s) * ones(1,2), ...
                [throughput_coded_confidence.lower(u, s), throughput_coded_confidence.upper(u, s)], ...
                 'LineStyle', '-', 'Color', cmp(u,:), 'HandleVisibility','off');
            plot( simulation_results.SNR_vector(s) * ones(1,2), ...
                [throughput_uncoded_confidence.lower(u, s), throughput_uncoded_confidence.upper(u, s)], ...
                 'LineStyle', '--', 'Color', cmp(u,:), 'HandleVisibility','off');
        end
    end
    
end

% hack to add an extra legend entry for uncoded
plot(NaN, NaN, 'LineStyle', '--', 'Color', ones(1,3)*0.5);
hold off;

legend_entries{N_UE+1} = sprintf('uncoded');
legend(legend_entries, 'Location', 'best');

%if simulation_results.conf_interval_probability
%    annotation('textbox',[.78 .8 .1 .1],...
%        'String',sprintf('Confidence %d%%',simulation_results.conf_interval_probability*100), 'FitBoxToText','on', 'HorizontalAlignment', 'right');
%end

xlabel('SNR [dB]');
ylabel('throughput [Mbit/s]');
title('throughput per user');
grid on;

if LTE_params.save_plots 
    print('report/throughput_user','-dpng')
end

end
