function [LLR_SD,M] = LTE_UL_detect_MUMIMO(MCS_and_scheduling,filtering,rx_user_symbols,H_est,LTE_params,receiver,sigma_n2,MSE,UE_output,uu,H_int)
% MU MIMI detection
% Author 
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
nUE = LTE_params.nUE;
nLayers = MCS_and_scheduling.nLayers;
M = MCS_and_scheduling.CQI_params.modulation_order;
nRX = LTE_params.BS_config.nRX;
bittable = LTE_params.bittable{M}; % Bitmapping table
symbol_alphabet = LTE_params.SymbolAlphabet{M}.'; % Symbol alphabet
freq_temp = unique(MCS_and_scheduling.freq_indices);

if (strcmp(filtering,'BlockFading'))

    H_complete = squeeze(H_est(:,1,:,:)); %h error
    Nsymbols = size(rx_user_symbols,1);
    user_symbols = size(rx_user_symbols,1);
    rx_layer_x = zeros(user_symbols,size(rx_user_symbols,2));

    switch receiver
        case 'ZF'  
            noise_enhancement_tmp = 0;
            for ctr2 = 1:size(rx_user_symbols,1)
                H_temp = reshape(H_int(freq_temp(ctr2),1,:,:,:),size(H_int,3),size(H_int,4)*size(H_int,5));
                col_norm = sum(abs((H_temp)'*H_temp),1);
                H_temp = H_temp(:,col_norm ~= 0);
                temp2 =  H_temp*(H_temp'*H_temp)^(-1);
                inv_chan_int = temp2*H_temp'; 
                temp_id = eye(size(inv_chan_int));
                Gk = temp_id-inv_chan_int;
                rx_temp = permute(reshape(rx_user_symbols(ctr2,:,:),size(rx_user_symbols,2),size(rx_user_symbols,3)),[2,1]);
                H_user = reshape(H_complete(ctr2,:,:),size(H_complete,2),size(H_complete,3));
                temp_sign = H_user'*(Gk)*rx_temp; %r_k
                c_temp = H_user'*Gk'*H_user;
                inv_temp_fin =(1/c_temp)*H_user'*(Gk'); %1*4)
                rx_layer_x(ctr2,:) = (1/c_temp)*temp_sign;
                if LTE_params.DFT_spreading_off
                    noise_enhancement_tmp(:,ctr2) = sum(abs(inv_temp_fin).^2,2);
                else
                    noise_enhancement_tmp = noise_enhancement_tmp + sum(abs(inv_temp_fin).^2,2);
                end       

            end                 
            Hg = 1;
    end

    if LTE_params.DFT_spreading_off
        rx_layer_x = rx_layer_x;
    else
        rx_layer_x = ifft(rx_layer_x,Nsymbols)*sqrt(Nsymbols);
    end
%     plot(real(rx_layer_x),imag(rx_layer_x),'bx')
%     hold on
%     grid on
       
    Nsymbols = size(rx_user_symbols,1)*size(rx_user_symbols,2);
    noise_enhancement = zeros(sum(M),Nsymbols);
    
    for ii=1:nLayers

        if LTE_params.DFT_spreading_off
            noise_enhancement(sum(M(1:ii-1)) + (1:M(ii)),:) = kron(ones(1,size(rx_user_symbols,2)),repmat(noise_enhancement_tmp(ii,:),length(1:M(ii)),1));
        else
            noise_enhancement(sum(M(1:ii-1)) + (1:M(ii)),:) = noise_enhancement_tmp(ii);
        end

    end
    
    
    if LTE_params.DFT_spreading_off
            noise_enhancement = noise_enhancement;
    else
            noise_enhancement =  1/Nsymbols*noise_enhancement;
    end
   
    rx_layer_x = transpose(reshape(rx_layer_x,Nsymbols,[])); 
    LLR_SD = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement,receiver);    
    
    
else 
    error('fast fading not working')
end

