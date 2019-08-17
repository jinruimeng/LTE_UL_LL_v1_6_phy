% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function cqi_user( simulation_results, LTE_params , uu)
%cqi_user plots each users CQI


cqi = simulation_results.UE_specific(uu).used_cqi;
n_streams = size(cqi, 3);

%cqi = cqi + rand(size(cqi))*0.3-0.15;


% UE_spec = simulation_results.UE_specific;
% N_UE = length(UE_spec);
N_snr = size(simulation_results.SNR_vector,2);
snr = repmat(simulation_results.SNR_vector,[LTE_params.N_subframes,1]);
SNR = simulation_results.SNR_vector;
delta_SNR = (SNR(end)-SNR(1))/length(SNR);
if delta_SNR == 0 %if only one snr point
    delta_SNR = SNR(end)/5;
end

fig = figure('Name', ['N_sub = ', num2str(LTE_params.N_subframes)]);

cqi1 = cqi(:,:,1);
if n_streams > 1  % how many are possible 2 or 4?
    cqi2 = cqi(:,:,2);
end

%include_legend = true;
handles = cell(1);

for aa = 1:n_streams
    CQI = cqi(:,:,aa);
    for nn = 1:N_snr
        bins = 1:1:10;
        classes = bins.^2/max(bins)^2;

        [bin_num,bin_label]=hist(CQI(:,nn),unique(CQI(:,nn)));
        relative_bins = bin_num/sum(bin_num);
        for ll = 1:length(bin_label)
            mysize = find(classes < relative_bins(ll), 1, 'last');
            if isempty(mysize)
                mysize = 1;

            end

            if bin_num(ll) ~= 0
                if aa == 1
                    
                    if n_streams > 1 
                        deltay = -0.15;
                    else
                        deltay = 0;
                    end
                    
                    handles{1} = plot(simulation_results.SNR_vector(nn), bin_label(ll)+deltay, 'ro', 'MarkerSize', mysize, 'MarkerFaceColor', 'r');
                else
                    handles{2} = plot(simulation_results.SNR_vector(nn), bin_label(ll)+0.15, 'bo', 'MarkerSize', mysize, 'MarkerFaceColor', 'b');
                end
                    %fprintf('%d/%d: %d \n', simulation_results.SNR_vector(nn), bin_label(ll), mysize);
            end
            hold on;

        end
    end
    %include_legend = true;
end

%if n_streams > 1
%    hold on;
%    scatter(snr(:), cqi2(:), 'ro');
%end

xlabel('SNR [dB]');
ylabel('CQI');
title('CQI per user');
grid on;
set(gca,'YTick',0:1:15);

xlim([SNR(1)-delta_SNR*0.5, SNR(end)+delta_SNR])
ylim([-0.5, 15.5]);
hold off;

if n_streams == 1
    legend([handles{1}], {'codeword 1'}, 'Location', 'SouthEast')
else
    legend([handles{1}, handles{2}], {'codeword 1', 'codeword 2'}, 'Location', 'SouthEast')
end



% 
% % calculate throughput per user
% BLER = zeros(N_UE, N_snr);
% 
% BLER_confidence = struct;
% BLER_confidence.lower = zeros(N_UE, N_snr) + eps;
% BLER_confidence.upper = zeros(N_UE, N_snr) + eps;
% 
% 
% 
% 
% % save the BLER per user
% BLER = zeros(N_UE, N_snr);
% 
% for u = 1:N_UE
%     if simulation_results.conf_interval_probability == 0 || isempty(fieldnames(UE_spec(u).confidence))
%         BLER(u, :) = UE_spec(u).BLER_overall;
%     else
%         BLER(u, :) = UE_spec(u).confidence.BLER_coded(1,:);
%         BLER_confidence.lower(u, :) = UE_spec(u).confidence.BLER_coded(2,:);
%         BLER_confidence.upper(u, :) = UE_spec(u).confidence.BLER_coded(3,:);
%     end
% end
% 
% 
% fig = figure('Name', ['N_sub = ', num2str(LTE_params.N_subframes)]);
% cmp = colormap(jet(N_UE));
% markers = plots.getCustomMarkers();
% legend_entries = cell(1, N_UE);
% 
% % plot curves 
% for u = 1:N_UE
%     marker = plots.getMarker(u, markers);
%     
%     semilogy( simulation_results.SNR_vector, BLER(u, :), ...
%         'LineStyle', '-', 'Marker', marker, 'Color', cmp(u,:) );
%     hold on;
%     
%     % add legend entry
%     if LTE_params.nBS > 1
%         b = utils.findBS( LTE_params.connection_table, u );
%         legend_entries{u} = sprintf('UE %d (BS %d)', u, b);
%     else
%         legend_entries{u} = sprintf('UE %d', u);
%     end
%     
%     if simulation_results.conf_interval_probability > 0 % confidence intervals computed
%         for s = 1:N_snr
%             plot( simulation_results.SNR_vector(s) * ones(1,2), ...
%                 [BLER_confidence.lower(u, s), BLER_confidence.upper(u, s)], ...
%                  'LineStyle', '-', 'Color', cmp(u,:), 'HandleVisibility','off');
%         end
%     end
%     
% end
% 
% hold off;
% 
% legend(legend_entries, 'Location', 'best');
% 
% 
% xlabel('SNR [dB]');
% ylabel('BLER');
% title('BLER per user');
% grid on;
% 
% if LTE_params.save_plots 
%     print('report/bler_user','-dpng')
% end

end
