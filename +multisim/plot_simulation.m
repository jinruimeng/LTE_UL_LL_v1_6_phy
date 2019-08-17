function plot_simulation( varargin )
% plot_simulation plots the results of a simulation
% modes: matlab, matlabandtikz, tikz

%close all;

simulation = varargin{1};

N_arg = length(varargin);

% default parameters
type = 'throughput_cell';
mode = 'matlab';
matlab_plot = true;
tikz_plot = false;
tikz_output_file = 'multisim_plots/plot1.tex';
select_subsims = false;
tikz_markerlist  = {'o', 'square*', 'triangle*', 'diamond*'};
tikz_linelist = {'solid', 'dashed', 'loosely dashed'};
tikz_ylimit = [];

matlab_markerlist = {'o', 'x', 'd', '^', 'v', '*' };
matlab_linelist = {'-', '--' };
% options: south west|south east|north west|north east
tikz_legend_style = 'legend style={at={(0.05,0.95)}, anchor=north west}, ';
coded_or_uncoded = 0;

if N_arg > 1
    if mod(N_arg, 2) == 0
        error(['plot_simulation always needs an odd number of arguments:' ...
               'simulation, key1, value1, key2, value2, ...']);
    end
    
    N_last = N_arg;
    current_ = 2;
    
    while (current_ <= N_last)
        key = varargin{current_};
        value = varargin{current_+1};
        
        switch (key)
            % choose the desired output format. either 'matlab' or 'tikz'
            % or 'matlabandtikz'
            case 'mode'
                available_modes = {'matlab', 'matlabandtikz', 'tikz'};
                mode = validatestring(value, available_modes); 
                if strcmp(mode, 'matlab')
                    matlab_plot = true;
                    tikz_plot = false;
                elseif strcmp(mode, 'matlabandtikz')
                    matlab_plot = true;
                    tikz_plot = true;
                elseif strcmp(mode, 'tikz')
                    matlab_plot = false;
                    tikz_plot = true;
                end
            
            % select only specific subsims to be shown (default is show
            % all subsimulations included in the simulation)
            case 'subsimulations'
                select_subsims = true;
                labels = value;
                
            % change the plot's title
            case 'title'
                simulation.name = value;
                
            % set the type of the plot
            case 'type'
                type = value;
                
            case 'tikz_ylimit'
                tikz_ylimit = value;
            
            % selects to which file the tikz code is written
            case 'tikzoutput'
                tikz_output_file = value;
                
            % change the legend position
            case 'tikzlegend'
                tikz_legend_style = ['legend pos = ' value, ', '];
            
            % decide if ber should be shown coded or uncoded
            % expects containers.Map of the form { name -> 'c', name -> 'u', name -> 'uc}
            % plot for subsim: c-> coded, u-> uncoded, uc -> both
            case 'coded_uncoded'
                coded_or_uncoded = value;
                 
            otherwise
                error('undefined option %s', key); 
            
        end
   
        current_ = current_ + 2;
    end
end
    

if ~select_subsims
    labels = simulation.subsimulations.keys;
end

subsims = simulation.subsimulations;
N_subsims = length(labels);
% enter coded uncoded decision here <-
if coded_or_uncoded == 0 % this is the default, no option is selected.
    coded_or_uncoded = containers.Map;
    for i = 1:N_subsims
        coded_or_uncoded(labels{i}) = 'u';
    end
end


legend_names = cell(1,N_subsims); % for the legend % *2 hack for ber

SNR = simulation.SNR_vec;

plotcolors = {'red', 'blue', 'green', 'magenta' };

if isempty(simulation.sweep_name) % not in sweep mode
    if matlab_plot
        figure();
    end
    
    if tikz_plot
        tikzlines = {};

        tikzlines{end+1} = '\begin{tikzpicture}[scale=\tikzscale]';
        tikzlines{end+1} = '\begin{axis}[';
        switch (type)
            case 'throughput_cell'
                tikzlines{end+1} =  'xlabel={SNR (dB)},ylabel={cell throughput (Mbit/s)}, ';
            case 'ber_cell'
                tikzlines{end+1} =  'xlabel={SNR (dB)},ylabel={BER},ymode=log,';
            case 'papr'
                tikzlines{end+1} =  'xlabel={PAPR (dB)},ylabel={eccdf},ymode=log,';
            otherwise
                error('undefined plot type %s', type);
        end
        
        if ~isempty(tikz_ylimit)
            tikzlines{end+1} =  sprintf('ymin=%f, ymax=%f,', tikz_ylimit(1), tikz_ylimit(2) ); 
        end
        
        tikzlines{end+1} =  tikz_legend_style;
        tikzlines{end+1} =  'grid=both]';
    end
    
    i_ = 1; % second index need in case there are more than one plot per subsim (eg.: BER coded/uncoded)
    j_ = 1; % same as above
    for i = 1:N_subsims
  
        if matlab_plot
            switch (type)
                case 'throughput_cell'
                    ptp = plot( SNR, subsims(labels{i}).simulation_results.cell_specific.confidence.throughput_coded(1,:)/1e3, ...
                    [matlab_markerlist{mod(i-1,length(matlab_markerlist))+1}, '-']);
                    legend_names{i} = labels{i};
                    
                case 'ber_cell'
                    if strcmp(coded_or_uncoded(labels{i}), 'u') % uncoded
                        ptp = semilogy( SNR, subsims(labels{i}).simulation_results.cell_specific.BER_uncoded, ...
                        [matlab_markerlist{mod(i_-1,length(matlab_markerlist))+1}, '-']);
                        
                        legend_names{i_} = ['uncoded ' labels{i}];
                        i_ = i_ + 1;
                    elseif strcmp(coded_or_uncoded(labels{i}), 'c') % coded
                        ptp = semilogy( SNR, subsims(labels{i}).simulation_results.cell_specific.BER_coded, ...
                        [matlab_markerlist{mod(i_-1,length(matlab_markerlist))+1}, '-']);
                        
                        legend_names{i_} = ['coded ', labels{i}];
                        i_ = i_ + 1;
                    else % coded + uncoded
                        ptp = semilogy( SNR, subsims(labels{i}).simulation_results.cell_specific.BER_uncoded, ...
                        [matlab_markerlist{mod(i_-1,length(matlab_markerlist))+1}, '-']);
                        legend_names{i_} = ['uncoded ', labels{i}];
                        i_ = i_ + 1;
                        hold on;
                        
                        ptp2 = semilogy( SNR, subsims(labels{i}).simulation_results.cell_specific.BER_coded, ...
                        [matlab_markerlist{mod(i_-1,length(matlab_markerlist))+1}, '-']);
                        
                        legend_names{i_} = ['coded ', labels{i}];
                        i_ = i_ + 1;
                        hasbehavior(ptp2,'legend',true);
                    end
                case 'papr'
                    papr1 = subsims(labels{i}).simulation_results.UE_specific.papr(:,1,:); % take always the first SNR vector (shouldn't matter). what is the third dimension??
                    papr2 = papr1(:);
                    [f, x] = ecdf(papr2);
                    ptp = semilogy(x, 1-f,  matlab_linelist{mod(i-1,length(matlab_linelist))+1}, 'LineWidth', 3);
                    legend_names{i} = labels{i};
                    
                    
                otherwise
                    error('undefined plot type %s', type);     
            end
            
            hasbehavior(ptp,'legend',true);
            hold on;
        end
        
        if tikz_plot 
            switch (type)
                % TODO optimise this copy paste code...
                case 'throughput_cell'
                    tikzlines{end+1} =  sprintf('\\addplot[smooth,mark=%s,%s, error bars/.cd, y dir=both, y explicit] plot coordinates {',...
                        tikz_markerlist{mod(i-1,length(tikz_markerlist))+1},plotcolors{mod(i-1,length(plotcolors))+1});

                    tp = subsims(labels{i}).simulation_results.cell_specific.confidence.throughput_coded(1,:)/1e3; 
                    conf = subsims(labels{i}).simulation_results.cell_specific.confidence.throughput_coded(2:3,:)/1e3;

                    for j=1:length(SNR)
                        tikzlines{end+1} = sprintf('(%f, %f) += (0,%f) -= (0,%f)', SNR(j), tp(j), abs(tp(j) - conf(1,j)),abs(tp(j) - conf(2,j))) ; 
                    end
                
                    tikzlines{end+1} =  '};';
                    tikzlines{end+1} = sprintf('\\addlegendentry{%s}', legend_names{i});
                    tikzlines{end+1} =  '';
            
                case 'ber_cell'
                    if strcmp(coded_or_uncoded(labels{i}), 'u')
                        tikzlines{end+1} =  sprintf('\\addplot[smooth,mark=%s,%s, error bars/.cd, y dir=both, y explicit] plot coordinates {',...
                        tikz_markerlist{mod(j_-1,length(tikz_markerlist))+1},plotcolors{mod(j_-1,length(plotcolors))+1});
                    
                        ber = subsims(labels{i}).simulation_results.cell_specific.BER_uncoded;
                        
                        for j=1:length(SNR)
                            tikzlines{end+1} = sprintf('(%f, %f)', SNR(j), ber(j)) ; 
                        end
                        tikzlines{end+1} =  '};';
                        tikzlines{end+1} = sprintf('\\addlegendentry{uncoded %s}', labels{i});
                        tikzlines{end+1} =  '';
                        j_ = j_ + 1;
                        
                    elseif strcmp(coded_or_uncoded(labels{i}), 'c')
                        tikzlines{end+1} =  sprintf('\\addplot[smooth,mark=%s,%s, error bars/.cd, y dir=both, y explicit] plot coordinates {',...
                        tikz_markerlist{mod(j_-1,length(tikz_markerlist))+1},plotcolors{mod(j_-1,length(plotcolors))+1});
                    
                        ber = subsims(labels{i}).simulation_results.cell_specific.BER_coded;
                        
                        for j=1:length(SNR)
                            tikzlines{end+1} = sprintf('(%f, %f)', SNR(j), ber(j)) ; 
                        end
                        
                        tikzlines{end+1} =  '};';
                        tikzlines{end+1} = sprintf('\\addlegendentry{coded %s}', labels{i});
                        tikzlines{end+1} =  '';
                        j_ = j_ + 1;
                    else % uncoded and coded
                        ber_uncoded = subsims(labels{i}).simulation_results.cell_specific.BER_uncoded;
                        ber_coded = subsims(labels{i}).simulation_results.cell_specific.BER_coded;
                        
                        for a_=1:2
                            if a_ == 1 %uncoded
                                tikzlines{end+1} =  sprintf('\\addplot[smooth,mark=%s,%s, error bars/.cd, y dir=both, y explicit] plot coordinates {',...
                                tikz_markerlist{mod(j_-1,length(tikz_markerlist))+1},plotcolors{mod(j_-1,length(plotcolors))+1});
                            else % coded
                                tikzlines{end+1} =  sprintf('\\addplot[smooth,mark=%s,%s, error bars/.cd, y dir=both, y explicit] plot coordinates {',...
                                tikz_markerlist{mod(j_-1,length(tikz_markerlist))+1},plotcolors{mod(j_-1,length(plotcolors))+1});
                            end
                            
                            for j=1:length(SNR)
                                if a_ == 1 %uncoded
                                    tikzlines{end+1} = sprintf('(%f, %f)', SNR(j), ber_uncoded(j)) ; 
                                else % coded
                                    tikzlines{end+1} = sprintf('(%f, %f)', SNR(j), ber_coded(j)) ; 
                                end
                            end
                            tikzlines{end+1} =  '};';
                            if a_ == 1 %uncoded
                                tikzlines{end+1} = sprintf('\\addlegendentry{uncoded %s}', labels{i});
                            else  % coded
                                tikzlines{end+1} = sprintf('\\addlegendentry{coded %s}', labels{i});
                            end
                            tikzlines{end+1} =  '';
                            j_ = j_ + 1;
                        end
                    end
                case 'papr'
                    papr1 = subsims(labels{i}).simulation_results.UE_specific.papr(:,1,:); % take always the first SNR vector (shouldn't matter)
                    papr2 = papr1(:);
                    [f, x] = ecdf(papr2);
                    f1 = 1-f;
                    downsampling = 30;
                    tikz_thicknesses = {'very thick', 'thick', 'semithick', 'thin' };
                    x = x(1:downsampling:end);
                    f1 = f1(1:downsampling:end);
                    
                    
                    
                    tikzlines{end+1} =  sprintf('\\addplot[smooth,%s,%s,%s] plot coordinates {',...
                        tikz_linelist{mod(i-1,length(tikz_linelist))+1}, plotcolors{mod(i-1,length(plotcolors))+1}, tikz_thicknesses{mod(i-1,length(tikz_thicknesses))+1});

                    for j=1:length(x)
                        tikzlines{end+1} = sprintf('(%f, %f)', x(j), f1(j)) ; 
                    end
                
                    tikzlines{end+1} =  '};';
                    tikzlines{end+1} = sprintf('\\addlegendentry{%s}', legend_names{i});
                    tikzlines{end+1} =  '';
                    
                otherwise
                    error('undefined plot type %s', type);  
            end
        end

        % plot confidence intervals
        switch (type)
            case 'throughput_cell'
                for j = 1:length(SNR)
                     if matlab_plot
                         pl = plot([SNR(j) SNR(j)],subsims(labels{i}).simulation_results.cell_specific.confidence.throughput_coded(2:3,j)/1e3 ,'k');
                         hasbehavior(pl,'legend',false);
                     end
                end
            case 'ber_cell'
                %
            case 'papr'
                %
            otherwise
                error('undefined plot type %s', type);
        end

    end
    
    if tikz_plot
        tikzlines{end+1} = '\end{axis}';
        tikzlines{end+1} = '\end{tikzpicture}';
    end

    if matlab_plot
        legend(legend_names, 'Location','best');
        title(simulation.name);
        xlabel('SNR (dB)');
        switch (type)
            case 'throughput_cell'
                ylabel('cell throughput (Mbit/s)');
            case 'ber_cell'
                ylabel('BER');
            case 'papr'
                ylabel('P (PAPR^ > PAPR)');
            otherwise
                    error('undefined plot type %s', type);
        end
        grid on;
        hold off;
    end
    
else % simulation is parameter sweep
    if matlab_plot
        figure();
    end
    
    if tikz_plot
        tikzlines = {};

        tikzlines{end+1} = '\begin{tikzpicture}';
        tikzlines{end+1} = '\begin{axis}[';
        tikzlines{end+1} =  sprintf('xlabel={%s},ylabel={cell throughput (Mbit/s)},', simulation.sweep_label);
        tikzlines{end+1} = tikz_legend_style;
        
        
        if ~isempty(tikz_ylimit)
            tikzlines{end+1} =  sprintf('ymin=%f, ymax=%f,', tikz_ylimit(1), tikz_ylimit(2) ); 
        end
        
        tikzlines{end+1} = 'grid=both]';
    end
    
    %todo: check if sweep really works with specific subsims
    N_sweep = length(simulation.sweep_values);
    for i = 1:N_subsims
        legend_names{i} = labels{i};
        
         if tikz_plot
            tikzlines{end+1} =  sprintf('\\addplot[smooth,mark=o,%s, error bars/.cd, y dir=both, y explicit] plot coordinates {',plotcolors{mod(i-1,3)+1});
         end
        
        raw_tp = zeros(size(simulation.sweep_values));
        
        tmp = subsims(labels{i});
        
        for s_ = 1:N_sweep
            current_sim = tmp(s_);
            raw_tp(s_) = current_sim.simulation_results.cell_specific.confidence.throughput_coded(1,:)/1e3; % assume only one SNR value
            conf_int = current_sim.simulation_results.cell_specific.confidence.throughput_coded(2:3,1)/1e3;
            if matlab_plot
                pl = plot(simulation.sweep_values(s_)*ones(1,2), conf_int, 'k');
                hasbehavior(pl,'legend',false);
                hold on;
            end
            
            if tikz_plot
                tikzlines{end+1} = sprintf('(%f, %f) += (0,%f) -= (0,%f)', simulation.sweep_values(s_), raw_tp(s_), ...
                abs(conf_int(2) - raw_tp(s_)), abs(conf_int(1) - raw_tp(s_))); 
            end
            
        end
        if matlab_plot
            ptp = plot(simulation.sweep_values, raw_tp, [matlab_markerlist{mod(i-1,length(matlab_markerlist))+1}, '-']);
            hasbehavior(ptp,'legend',true);
            hold on;
        end
        
        if tikz_plot
            tikzlines{end+1} =  '};';
            tikzlines{end+1} = sprintf('\\addlegendentry{%s}', legend_names{i});
            tikzlines{end+1} =  '';
        end

    end
    
    if tikz_plot
        tikzlines{end+1} = '\end{axis}';
        tikzlines{end+1} = '\end{tikzpicture}';
    end

    if matlab_plot
        legend(legend_names, 'Location','best');
        title(simulation.name);
        xlabel(simulation.sweep_label);
        switch (type)
            case 'throughput_cell'
                ylabel('cell throughput (Mbit/s)');
            case 'ber_cell'
                ylabel('BER');
            otherwise
                    error('undefined plot type %s', type);
        end
        grid on;
        hold off;
    end

     
end
    

% if isempty(simulation.sweep_name) % only if one simulation
%     % only show MSE plot if channel est. is not perfect
%     if strcmp(subsims{1}.LTE_params.BS_config.channel_estimation_method, 'PERFECT') == 0 
%         figure(2);
% 
%         for i = 1:N_subsims
%             legend_names{i} = labels{i};
% 
% 
%             MSE = subsims{i}.simulation_results.cell_specific.MSE_overall;
% 
%             ptp = semilogy(SNR,MSE, 'o-');
%             hasbehavior(ptp,'legend',true);
%             hold on;
% 
%             % todo confidence intervals for MSE
%         end
%     end
% 
% 
%     legend(legend_names, 'Location','northwest');
%     title(simulation.name);
%     xlabel('SNR (dB)');
%     ylabel('MSE');
%     grid on;
%     hold off;
%     % 
% 
% end


if tikz_plot
    fileID = fopen(tikz_output_file, 'w');

    [nrows,ncols] = size(tikzlines);
    formatSpec = '%s\n';

    for row = 1:nrows
        fprintf(fileID,formatSpec,tikzlines{row,:});
    end

    fclose(fileID);
    fprintf('wrote tikz code to "%s"\n', tikz_output_file);
end


end












