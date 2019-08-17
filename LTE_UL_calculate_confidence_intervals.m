function [ simulation_results ] = LTE_UL_calculate_confidence_intervals(simulation_results, conf_probability)
% Calculates the confidence intervals
% Lukas Nagel, lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

if simulation_results.N_subframes < 3 % too short to compute intervals with bootci
    warning('no confidence intervals computed (too short)');
    return
end

alpha = 1-conf_probability;


for bb = 1:simulation_results.nBS
    simulation_results.cell_specific(bb).confidence = struct;
    
    tc = sum(simulation_results.cell_specific(bb).throughput_coded,3);
    simulation_results.cell_specific(bb).confidence.throughput_coded = [mean(tc,1); bootci(2000, {@mean, tc}, 'alpha',alpha)];
    
    tu = sum(simulation_results.cell_specific(bb).throughput_uncoded,3);
    simulation_results.cell_specific(bb).confidence.throughput_uncoded = [mean(tu,1); bootci(2000, {@mean, tu}, 'alpha',alpha)];
end

for uu = 1:(simulation_results.nUE * simulation_results.nBS)
    simulation_results.UE_specific(uu).confidence = struct;
    
    tc = sum(simulation_results.UE_specific(uu).throughput_coded,3);
    simulation_results.UE_specific(uu).confidence.throughput_coded = [mean(tc,1); bootci(2000, {@mean, tc}, 'alpha',alpha)];
    
    tu = sum(simulation_results.UE_specific(uu).throughput_uncoded,3);
    simulation_results.UE_specific(uu).confidence.throughput_uncoded = [mean(tu,1); bootci(2000, {@mean, tu}, 'alpha',alpha)];
    
    ber_func = @(x, y) sum(x,1)./sum(y,1);
    
    tf = sum(simulation_results.UE_specific(uu).FER_coded,3);
    tcode = sum(simulation_results.UE_specific(uu).used_codewords,3);
    
    if ~any(max(sum(simulation_results.UE_specific(uu).used_codewords,3)) ==0)
        simulation_results.UE_specific(uu).confidence.BLER_coded = [ber_func(tf, tcode); bootci(2000, {ber_func, tf, tcode},'alpha',alpha)];
    else
        simulation_results.UE_specific(uu).confidence.BLER_coded = zeros(3,size(tf,2));
    end
        
    tber_coded = sum(double(simulation_results.UE_specific(uu).biterrors_coded),3);
    tber_uncoded = sum(double(simulation_results.UE_specific(uu).biterrors_uncoded),3);
    
    tblock_coded = sum(double(simulation_results.UE_specific(uu).blocksize_coded),3); 
    tblock_uncoded = sum(double(simulation_results.UE_specific(uu).blocksize_uncoded),3); 
    
    if ~any(max(sum(simulation_results.UE_specific(uu).blocksize_coded,3)) == 0)
        simulation_results.UE_specific(uu).confidence.BER_coded = [ber_func(tber_coded, tblock_coded); bootci(2000, {ber_func, tber_coded, tblock_coded},'alpha',alpha)];
        simulation_results.UE_specific(uu).confidence.BER_uncoded = [ber_func(tber_uncoded, tblock_uncoded); bootci(2000, {ber_func, tber_uncoded, tblock_uncoded},'alpha',alpha)];
    else
        
    end
    channel_error = simulation_results.UE_specific(uu).channel_error;
    simulation_results.UE_specific(uu).confidence.channel_error = [mean(channel_error,1); bootci(2000, {@mean, channel_error}, 'alpha',alpha)];
    
    channel_pred_error = simulation_results.UE_specific(uu).channel_pred_error;
    simulation_results.UE_specific(uu).confidence.channel_pred_error = [mean(channel_pred_error,1); bootci(2000, {@mean, channel_pred_error}, 'alpha',alpha)];
    
end


end




