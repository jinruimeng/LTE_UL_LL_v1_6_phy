% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function papr_user( simulation_results , LTE_params)
%papr plots the papr (of the first user)

N_UE = LTE_params.nUE * LTE_params.nBS;

fig = figure('Name', ['N_sub = ', num2str(LTE_params.N_subframes)]);
cmp = colormap(jet(N_UE));
legend_entries = cell(1, N_UE);

for u = 1:N_UE
    if ~isempty(simulation_results.UE_specific(u).papr)
        papr1 = simulation_results.UE_specific(u).papr(:,1,:); % take always the first SNR vector (shouldn't matter).
        papr2 = papr1(:);
        [f, x] = ecdf(papr2);
        semilogy(x, 1-f,  '-', 'LineWidth', 1,  'Color', cmp(u,:));

        hold on;
    end
    
    if LTE_params.nBS > 1
        b = utils.findBS( LTE_params.connection_table, u );
        legend_entries{u} = sprintf('UE %d (BS %d)', u, b);
    else
        legend_entries{u} = sprintf('UE %d', u);
    end
end    
    


xlabel('PAPR [dB]');
ylabel('P (PAPR^ > PAPR)');
legend(legend_entries, 'Location', 'best');
grid on;

title('PAPR');



if LTE_params.save_plots 
    print('report/papr_user','-dpng')
end

end


