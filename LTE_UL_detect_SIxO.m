function [LLR_SD,M] = LTE_UL_detect_SIxO(MCS_and_scheduling,filtering,rx_user_symbols,H_est,LTE_params,receiver,sigma_n2,UE_output,uu, bb)
% SIxO detection
%
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

nLayers         = MCS_and_scheduling.nLayers;
M               = MCS_and_scheduling.CQI_params.modulation_order;
nRX             = LTE_params.BS_config.nRX;
bittable        = LTE_params.bittable{M};           % Bitmapping table
symbol_alphabet = LTE_params.SymbolAlphabet{M}.';   % Symbol alphabet

% generate UE specific mapping tables
CH_mapping_UE_spec = UE_output(uu).UE_genie.CH_mapping(logical(kron(MCS_and_scheduling.UE_mapping,ones(LTE_params.Nsc,LTE_params.Ns))));
CH_mapping_UE_spec = repmat(reshape(CH_mapping_UE_spec, [], LTE_params.Nsub),[1 1 nRX]);

% if nRX > 1
%     error('not working properly with more than one receive antenna; the implemented receivers assume single antenna');
% end

if (strcmp(filtering,'BlockFading'))
    H_complete = squeeze(H_est(:,1,:,:));
    Nsymbols = size(rx_user_symbols,1)*size(rx_user_symbols,2);
    user_symbols = reshape(rx_user_symbols,Nsymbols,[]);
    Hg_full = zeros(1,Nsymbols);
    switch receiver
        case 'ZF'
            inv_temp = conj(H_complete)./repmat(sum(H_complete.*conj(H_complete),2),1,nRX);
            Hg = 1;
        case 'MMSE'
            inv_temp = repmat((sum(abs(H_complete).^2,2) + sigma_n2).^(-1),[1 size(H_complete,2)]) .* conj(H_complete);
            Hg = mean(sum(inv_temp.*H_complete,2));
            Hg_full = sum(inv_temp.*H_complete,2).';
    end    
    
    inv_temp = repmat(inv_temp,size(rx_user_symbols,2),1);
    rx_layer_x = sum(inv_temp.*user_symbols,2);
    rx_layer_x = reshape(rx_layer_x,size(rx_user_symbols,1),[]);
    
     if LTE_params.DFT_spreading_off    % OFDMA
        rx_layer_x = rx_layer_x;
        noise_enhancement = sigma_n2*repmat(transpose(sum(abs(inv_temp).^2,2)),M,1);
     else                               % SC-FDMA
        rx_layer_x = ifft(rx_layer_x,size(rx_user_symbols,1))*sqrt(size(rx_user_symbols,1));
        noise_enhancement = sigma_n2*mean(sum(abs(inv_temp).^2,2),1)*ones(M,Nsymbols);
     end
    
     if LTE_params.DFT_spreading_off && ~strcmp(receiver,'ZF')  
        LLR_SD = zeros(size(noise_enhancement));
        s1 = size(rx_layer_x,1);
        s2 = size(rx_layer_x,2);
        s3 = size(rx_layer_x,3);
        for ii = 1:s1
            LLR_temp = LTE_demapper(reshape(rx_layer_x(ii,:,:),s2,s3).',symbol_alphabet,bittable,nLayers,M,Hg_full(:,ii),noise_enhancement(:,ii:s1:end),receiver);  
            LLR_SD(:,ii:s1:end) = LLR_temp;
        end
     else
        rx_layer_x = reshape(rx_layer_x,1,[]);
        LLR_SD = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement,receiver);
     end
    
    rx_layer_x = reshape(rx_layer_x,size(rx_user_symbols,1),size(rx_user_symbols,2));
    
else % FastFading
    LLR_SD = zeros(M,length(rx_user_symbols));   % Log likelihood Ratios of the Sphere decoder
    index = 0;
    rx_layer_x = nan(size(rx_user_symbols,1),size(rx_user_symbols,2));
    for ctr = 1:size(rx_user_symbols,2)
        H_est_user=reshape(H_est(CH_mapping_UE_spec),size(H_est,1),[],nRX);
        H_complete = squeeze(H_est_user(:,ctr,:,:));
        temp = conj(H_complete);
        user_symbols = squeeze(rx_user_symbols(:,ctr,:));
        Nsymbols = size(rx_user_symbols,1);
        
        switch receiver
            case 'ZF'
                inv_temp = temp./repmat(sum(H_complete.*temp,2),1,nRX);
                Hg = 1;
            case 'MMSE'
                inv_temp = temp./repmat(sum(H_complete.*temp+sigma_n2,2),1,nRX);
                Hg = mean(sum(inv_temp.*H_complete,2));
        end

        rx_layer_x_temp = transpose(sum(inv_temp.*user_symbols,2));

        if LTE_params.DFT_spreading_off
            rx_layer_x_temp = rx_layer_x_temp;
            noise_enhancement = sigma_n2*repmat(transpose(sum(abs(inv_temp).^2,2)),M,1);
        else
            rx_layer_x_temp = ifft(rx_layer_x_temp,Nsymbols)*sqrt(Nsymbols); 
            noise_enhancement = sigma_n2*mean(sum(abs(inv_temp).^2,2),1)*ones(M,Nsymbols);
        end

        LLR_SD(:,index+1:index+Nsymbols) = LTE_demapper(rx_layer_x_temp,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);%,receiver);
        index = index+Nsymbols;

        rx_layer_x(:,ctr)=rx_layer_x_temp;
        rx_layer_x = reshape(rx_layer_x,size(rx_user_symbols,1),size(rx_user_symbols,2));
    
    end
    
end
% 
% scatter(real(rx_layer_x(:)),imag(rx_layer_x(:)));
% grid on;
