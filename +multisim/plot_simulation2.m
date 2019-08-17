function plot_simulation2( varargin )
% plot_simulation plots the results of a simulation


simulation = varargin{1};


% default parameters
%type = 'throughput_cell';



% if N_arg > 1
%     if mod(N_arg, 2) == 0
%         error(['plot_simulation always needs an odd number of arguments:' ...
%                'simulation, key1, value1, key2, value2, ...']);
%     end
%     
%     N_last = N_arg;
%     current_ = 2;
%     
%     while (current_ <= N_last)
%         key = varargin{current_};
%         value = varargin{current_+1};
%         
%         switch (key)
%             
%             % select only specific subsims to be shown (default is show
%             % all subsimulations included in the simulation)
%             case 'subsimulations'
%                 select_subsims = true;
%                 labels = value;
%                 
%             % change the plot's title
%             case 'title'
%                 simulation.name = value;
%                 
%             % set the type of the plot
%             case 'type'
%                 type = value;
%        
%             otherwise
%                 error('undefined option %s', key); 
%             
%         end
%    
%         current_ = current_ + 2;
%     end
% end
%     


subsims = simulation.subsimulations;

multisim.sim_throughput_user(simulation);
multisim.sim_ber_user(simulation);
multisim.sim_bler_user(simulation);
multisim.sim_throughput_cell(simulation);

tmp = simulation.subsimulations.values;

if ~strcmp(tmp{1}.LTE_params.BS_config.channel_estimation_method,'PERFECT') 
    multisim.sim_mse_user(simulation);
end







end












