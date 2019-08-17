% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function mse_user( simulation_results, LTE_params )
%mse_user plots each users estimation mse

UE_spec = simulation_results.UE_specific;
N_UE = length(UE_spec);
N_snr = size(simulation_results.SNR_vector,2);

% calculate mse per user
ce = zeros(N_UE, N_snr);

ce_confidence = struct;
ce_confidence.lower = zeros(N_UE, N_snr);
ce_confidence.upper = zeros(N_UE, N_snr);

for u = 1:N_UE
    if simulation_results.conf_interval_probability == 0 || isempty(fieldnames(UE_spec(u).confidence))
        ce_tmp = UE_spec(u).channel_error;
        ce(u, :) = mean(ce_tmp,1); 

    else
        ce(u, :) = UE_spec(u).confidence.channel_error(1,:);
        ce_confidence.lower(u, :) = UE_spec(u).confidence.channel_error(2,:);
        ce_confidence.upper(u, :) = UE_spec(u).confidence.channel_error(3,:);
        

    end
end


fig = figure('Name', ['N_sub = ', num2str(LTE_params.N_subframes)]);
cmp = colormap(jet(N_UE));
markers = plots.getCustomMarkers();
legend_entries = cell(1, N_UE);

% plot curves 
for u = 1:N_UE
    marker = plots.getMarker(u, markers);
    
    semilogy( simulation_results.SNR_vector(:), ce(u, :), ...
        'LineStyle', '-', 'Marker', marker, 'Color', cmp(u,:) );
    hold on;

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
                [ce_confidence.lower(u, s), ce_confidence.upper(u, s)], ...
                 'LineStyle', '-', 'Color', cmp(u,:), 'HandleVisibility','off');
            plot( simulation_results.SNR_vector(s) * ones(1,2), ...
                [ce_confidence.lower(u, s), ce_confidence.upper(u, s)], ...
                 'LineStyle', '--', 'Color', cmp(u,:), 'HandleVisibility','off');
        end
    end
    
end

hold off;

legend(legend_entries, 'Location', 'best');


xlabel('SNR [dB]');
ylabel('MSE');
title('estimation MSE per user');
grid on;


if LTE_params.save_plots 
    print('report/ber_user','-dpng')
end

end
