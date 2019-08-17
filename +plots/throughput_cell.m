% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function throughput_cell( simulation_results, LTE_params )
%throuphput_cell plots each cells throughput

cell_spec = simulation_results.cell_specific;
N_BS = length(cell_spec);
N_snr = size(simulation_results.SNR_vector,2);

% calculate throughput per cell
throughput_coded = zeros(N_BS, N_snr);
throughput_uncoded = zeros(N_BS, N_snr);

throughput_coded_confidence = struct;
throughput_coded_confidence.lower = zeros(N_BS, N_snr);
throughput_coded_confidence.upper = zeros(N_BS, N_snr);

throughput_uncoded_confidence = struct;
throughput_uncoded_confidence.lower = zeros(N_BS, N_snr);
throughput_uncoded_confidence.upper = zeros(N_BS, N_snr);

for b = 1:N_BS
    if simulation_results.conf_interval_probability == 0 || isempty(fieldnames(cell_spec(b).confidence))
        coded_tmp = cell_spec(b).throughput_coded;
        throughput_coded(b, :) = mean(sum(coded_tmp,3),1) / 1e3; % tp in Mbit/s
        
        uncoded_tmp = cell_spec(b).throughput_uncoded;
        throughput_uncoded(b, :) = mean(sum(uncoded_tmp,3),1) / 1e3; % tp in Mbit/s
        
    else
        throughput_coded(b, :) = cell_spec(b).confidence.throughput_coded(1,:) / 1e3;
        throughput_coded_confidence.lower(b, :) = cell_spec(b).confidence.throughput_coded(2,:) / 1e3;
        throughput_coded_confidence.upper(b, :) = cell_spec(b).confidence.throughput_coded(3,:) / 1e3;
        
        throughput_uncoded(b, :) = cell_spec(b).confidence.throughput_uncoded(1,:) / 1e3;
        throughput_uncoded_confidence.lower(b, :) = cell_spec(b).confidence.throughput_uncoded(2,:) / 1e3;
        throughput_uncoded_confidence.upper(b, :) = cell_spec(b).confidence.throughput_uncoded(3,:) / 1e3;
    end
        
    
end


fig = figure('Name', ['N_sub = ', num2str(LTE_params.N_subframes)]);
cmp = colormap(winter(N_BS)); 
markers = plots.getCustomMarkers();
legend_entries = cell(1, N_BS);

SNR_vec = simulation_results.SNR_vector(1, :);

% plot curves 
for b = 1:N_BS
    marker = plots.getMarker(b, markers);
    
    plot( SNR_vec, throughput_coded(b, :), ...
        'LineStyle', '-', 'Marker', marker, 'Color', cmp(b,:) );
    hold on;
    plot( SNR_vec, throughput_uncoded(b, :), ...
        'LineStyle', '--', 'Marker', marker, 'Color', cmp(b,:), 'HandleVisibility','off');
    
    if simulation_results.conf_interval_probability > 0 % confidince intervals computed
        for s = 1:N_snr
            plot( SNR_vec(s) * ones(1,2), ...
                [throughput_coded_confidence.lower(b, s), throughput_coded_confidence.upper(b, s)], ...
                 'LineStyle', '-', 'Color', cmp(b,:), 'HandleVisibility','off');
            plot( SNR_vec(s) * ones(1,2), ...
                [throughput_uncoded_confidence.lower(b, s), throughput_uncoded_confidence.upper(b, s)], ...
                 'LineStyle', '--', 'Color', cmp(b,:), 'HandleVisibility','off');
        end
    end
    
    % add legend entry
    legend_entries{b} = sprintf('BS %d', b);

    
end

% hack to add an extra legend entry for uncoded
plot(NaN, NaN, 'LineStyle', '--', 'Color', ones(1,3)*0.5);
hold off;

legend_entries{N_BS+1} = sprintf('uncoded');
legend(legend_entries, 'Location', 'best');

if simulation_results.conf_interval_probability
    annotation('textbox',[.78 .8 .1 .1],...
        'String',sprintf('Confidence %d%%',simulation_results.conf_interval_probability*100), 'FitBoxToText','on', 'HorizontalAlignment', 'right');
end

xlabel('SNR [dB]');
ylabel('throughput [Mbit/s]');
title('throughput per cell');
grid on;


if LTE_params.save_plots 
    print('report/throughput_cell','-dpng')
end
end
