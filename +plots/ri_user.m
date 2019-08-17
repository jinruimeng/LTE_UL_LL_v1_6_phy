% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function ri_user( simulation_results, LTE_params , uu)
%cqi_user plots each users RI


% TODO

ri = simulation_results.UE_specific(uu).used_RI;

N_snr = size(simulation_results.SNR_vector,2);
SNR = simulation_results.SNR_vector;
delta_SNR = (SNR(end)-SNR(1))/length(SNR);

if delta_SNR == 0 %if only one snr point
    delta_SNR = SNR(end)/5;
end

fig = figure('Name', ['N_sub = ', num2str(LTE_params.N_subframes)]);


handles = cell(1);

for nn = 1:N_snr
    bins = 1:1:10;
    classes = bins.^2/max(bins)^2;

    [bin_num,bin_label]=hist(ri(:,nn),unique(ri(:,nn)));
    relative_bins = bin_num/sum(bin_num);
    for ll = 1:length(bin_label)
        mysize = find(classes < relative_bins(ll), 1, 'last');
        if isempty(mysize)
            mysize = 1;

        end

        if bin_num(ll) ~= 0
            handles{1} = plot(simulation_results.SNR_vector(nn), bin_label(ll), 'ro', 'MarkerSize', mysize, 'MarkerFaceColor', 'r');
                %fprintf('%d/%d: %d \n', simulation_results.SNR_vector(nn), bin_label(ll), mysize);
        end
        hold on;

    end
end



xlabel('SNR [dB]');
ylabel('RI');
title(['RI of user ',  num2str(uu)]);
grid on;
set(gca,'YTick',1:1:4);

xlim([SNR(1)-delta_SNR*0.5, SNR(end)+delta_SNR])
ylim([0.5, 4.5]);
hold off;


%legend([handles{1}], {'codeword 1'}, 'Location', 'SouthEast')


