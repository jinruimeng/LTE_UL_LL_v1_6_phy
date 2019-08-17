% Author: Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

function LTE_UL_plot_results( varargin ) %simulation_results, LTE_params )
%LTE_UL_plot_results plots a summary of the simulation_results

if nargin < 2
    error('LTE_UL_plot_results needs at least two arguments.');
else
    simulation_results = varargin{1};
    LTE_params = varargin{2};
    
    if nargin == 2
        % plot everything
        plots.throughput_user(simulation_results, LTE_params);
        plots.throughput_cell(simulation_results, LTE_params);
        plots.bler_user(simulation_results, LTE_params);
        plots.ber_user(simulation_results, LTE_params);
        plots.papr_user(simulation_results, LTE_params);
        if ~strcmp(LTE_params.BS_config.channel_estimation_method,'PERFECT') 
            plots.mse_user(simulation_results, LTE_params);
        end
        if LTE_params.BS_config.channel_prediction
            plots.pred_mse_user(simulation_results, LTE_params);
        end
    else
        % assume that the user wants to select plots
        to_plot = varargin{3};
        for i = 1:length(to_plot)
            switch to_plot{i}
                case 'throughput_user'
                    plots.throughput_user(simulation_results, LTE_params);
                case 'throughput_cell'
                    plots.throughput_cell(simulation_results, LTE_params);
                case 'bler_user'
                    plots.bler_user(simulation_results, LTE_params);
                case 'ber_user'
                    plots.ber_user(simulation_results, LTE_params);
                case 'mse_user'
                    if ~strcmp(LTE_params.BS_config.channel_estimation_method,'PERFECT') 
                        plots.mse_user(simulation_results, LTE_params);
                    else
                        fprintf('no MSE plot generated because "channel_estimation_method" is "PERFECT"\n');
                    end
                case 'pred_mse_user'
                    if LTE_params.BS_config.channel_prediction 
                        plots.pred_mse_user(simulation_results, LTE_params);
                    end
                case 'papr'
                    plots.papr_user(simulation_results, LTE_params);
                otherwise
                    warning('unknown argument "%s"', to_plot{i});
            end
        end
    end
end


end

