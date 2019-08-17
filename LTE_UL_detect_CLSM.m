function [LLR_SD,M] = LTE_UL_detect_CLSM(MCS_and_scheduling,rx_user_symbols,H_est,LTE_params,filtering,receiver,sigma_n2,UE_output,uu, bb)
% Open-loop Spatial multiplexing detection.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

nLayers = MCS_and_scheduling.nLayers;
W = UE_output(1,uu).LayerMapping;

% generate UE specific mapping tables
CH_mapping_UE_spec = UE_output(uu).UE_genie.CH_mapping(logical(kron(MCS_and_scheduling.UE_mapping,ones(LTE_params.Nsc,LTE_params.Ns))));
CH_mapping_UE_spec = reshape(CH_mapping_UE_spec, [], LTE_params.Nsub);


M = [MCS_and_scheduling.CQI_params.modulation_order];
switch nLayers % when layer number is unequal to codeword number we need to do something
    case 2
        if(MCS_and_scheduling.nCodewords == 1) % 1 codewords, 2 layers
            M = [M,M];
        end
    case 3  % 2 codewords, 3 layers
        M = [M(1),M(2),M(2)];
    case 4  % 2 codewords, 4 layers
        M = [M(1),M(1),M(2),M(2)];
end
bittable = false(sum(M(1:nLayers)),2^max(M));
symbol_alphabet = zeros(nLayers,2^max(M));
for i = 1:nLayers
    bittable(sum(M(1:i-1))+(1:M(i)),1:2^M(i))=LTE_params.bittable{M(i)}; % Bitmapping table
    symbol_alphabet(i,1:2^M(i))=LTE_params.SymbolAlphabet{M(i)}.'; % Symbol alphabet
end
LLR_SD = zeros(sum(M),length(rx_user_symbols));   % Log likelihood Ratios of the Spere decoder

l = 1:length(rx_user_symbols);
p = mod(l-1,nLayers)+1;                 % needed for the choice of the precoding matrix (see 3GPP TS 36.213 section 7.1.3)
k = mod(floor((l-1)/nLayers),4)+1;

index = 0;
nRX = size(H_est,3);
nTX = size(H_est,4);

%===========================
if (strcmp(filtering,'BlockFading'))
    H_complete = squeeze(H_est(:,1,:,:));
    Nsymbols = size(rx_user_symbols,1);
    rx_layer_x = zeros(size(rx_user_symbols,1),size(rx_user_symbols,2),nLayers);
    
   
    if LTE_params.DFT_spreading_off
        noise_enhancement_tmp = zeros(nLayers,Nsymbols);
    else
        noise_enhancement_tmp = zeros(nLayers,1);
    end
    
    switch receiver
        case 'ZF'
            
            Hg = eye(nLayers);
            
            for ctr2 = 1:Nsymbols
                inv_temp = pinv(reshape(squeeze(H_complete(ctr2,:,:)),[],nLayers));

                rx_layer_x(ctr2,:,:) = squeeze(rx_user_symbols(ctr2,:,:))*transpose(inv_temp);

                if LTE_params.DFT_spreading_off
                    noise_enhancement_tmp(:,ctr2) = sum(abs(inv_temp).^2,2);
                else
                    noise_enhancement_tmp = noise_enhancement_tmp + sum(abs(inv_temp).^2,2);
                end
                
            end
            
            
        case 'MMSE'
            Hg = zeros(nLayers);
            if LTE_params.DFT_spreading_off
                Hg_full = zeros(nLayers,nLayers,Nsymbols);
            end
            for ctr2 = 1:Nsymbols
                Hsc = reshape(squeeze(H_complete(ctr2,:,:)),[],nLayers);        %ezoechma: reshape as for the ZF
%                 temp = Hsc';
%                 inv_temp = temp/(temp*Hsc+sigma_n2*nTX*eye(nRX));
                inv_temp = (Hsc'*Hsc+sigma_n2*eye(nLayers))^(-1)*Hsc';
                rx_layer_x(ctr2,:,:) = squeeze(rx_user_symbols(ctr2,:,:))*transpose(inv_temp);
                
                if LTE_params.DFT_spreading_off
                    noise_enhancement_tmp(:,ctr2) = sum(abs(inv_temp).^2,2);
                else
                    noise_enhancement_tmp = noise_enhancement_tmp + sum(abs(inv_temp).^2,2);
                end
                Hg_full(:,:,ctr2) = inv_temp*Hsc;
                Hg = Hg + inv_temp*Hsc;                
            end
            Hg = Hg/Nsymbols;
    end
    
    
    

    if LTE_params.DFT_spreading_off
        rx_layer_x = rx_layer_x;
    else
        rx_layer_x = ifft(rx_layer_x,Nsymbols)*sqrt(Nsymbols);
    end
    
    Nsymbols = size(rx_user_symbols,1)*size(rx_user_symbols,2);
    noise_enhancement = zeros(sum(M),Nsymbols);
    
    for ii=1:nLayers
        if LTE_params.DFT_spreading_off
            noise_enhancement(sum(M(1:ii-1)) + (1:M(ii)),:) = kron(ones(1,size(rx_user_symbols,2)),repmat(noise_enhancement_tmp(ii,:),length(1:M(ii)),1));
        else
            noise_enhancement(sum(M(1:ii-1)) + (1:M(ii)),:) = noise_enhancement_tmp(ii);
        end

    end
    
    %     noise_enhancement_tmp = sigma_n2/Nsymbols*noise_enhancement_tmp;
    
    if LTE_params.DFT_spreading_off
            noise_enhancement = sigma_n2*noise_enhancement;
    else
            noise_enhancement =  sigma_n2/Nsymbols*noise_enhancement;
    end
    
    
    if LTE_params.DFT_spreading_off && ~strcmp(receiver,'ZF')   
        LLR_SD = zeros(size(noise_enhancement));
        s1 = size(rx_layer_x,1);
        s2 = size(rx_layer_x,2);
        s3 = size(rx_layer_x,3);
        for ii = 1:s1
            LLR_temp = LTE_demapper(reshape(rx_layer_x(ii,:,:),s2,s3).',symbol_alphabet,bittable,nLayers,M,Hg_full(:,:,ii),noise_enhancement(:,ii:s1:end),receiver);  
            LLR_SD(:,ii:s1:end) = LLR_temp;
        end
    else
        rx_layer_x = transpose(reshape(rx_layer_x,Nsymbols,[]));
        LLR_SD = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement,receiver);        
    end
        
%     plot(real(rx_layer_x(1,:)),imag(rx_layer_x(1,:)),'bx','linewidth',2)
%     hold on
%     grid on
    
    
    
%======================================
else  %FastFading
    for ctr = 1:size(rx_user_symbols,2)                                         %perform decoding for each (time) symbol
        %         error('H_est_user -> H_est');
        for tt=1:size(H_est,4)
            for rr=1:size(H_est,3)
                H_temp=squeeze(H_est(:,:,rr,tt));
%                 H_est_user(:,:,rr,tt)=reshape(H_temp(UE_output.UE_genie.CH_mapping),size(H_est,1),[]);
                H_est_user(:,:,rr,tt)=reshape(H_temp(CH_mapping_UE_spec),size(H_est,1),[]);

            end
        end
        H_complete = squeeze(H_est_user(:,ctr,:,:));                   % H_complete including precoder!
        user_symbols = squeeze(rx_user_symbols(:,ctr,:));
        Nsymbols = size(H_complete,1);                                          %NSYMBOLS = how many carriers are in one time instant
        rx_layer_x = zeros(Nsymbols,nLayers);
        
        if LTE_params.DFT_spreading_off
            noise_enhancement_tmp = zeros(nLayers,Nsymbols);
        else
            noise_enhancement_tmp = zeros(nLayers,1);
        end
  
        switch receiver
            case 'ZF'
                Hg =eye(nLayers);
                for ctr2 = 1:Nsymbols                %for each subcarrier
                    inv_temp = pinv(reshape(squeeze(H_complete(ctr2,:,:)),[],nLayers)); %ezoechma: reshape as above
                    rx_layer_x(ctr2,:) = user_symbols(ctr2,:)*transpose(inv_temp);                    
                    
                    if LTE_params.DFT_spreading_off
                        noise_enhancement_tmp(:,ctr2) = sum(abs(inv_temp).^2,2);
                    else
                        noise_enhancement_tmp = noise_enhancement_tmp + sum(abs(inv_temp).^2,2);
                    end                
                    
                end
            case 'MMSE'
                Hg = zeros(nLayers);
                for ctr2 = 1:Nsymbols
                    Hsc = reshape(squeeze(H_complete(ctr2,:,:)),[],nLayers);    %ezoechma
                    temp = Hsc';
                    inv_temp = temp/(Hsc*temp+nTX*sigma_n2*eye(nRX));
                    rx_layer_x(ctr2,:) = user_symbols(ctr2,:)*transpose(inv_temp);
     
                    if LTE_params.DFT_spreading_off
                        noise_enhancement_tmp(:,ctr2) = sum(abs(inv_temp).^2,2);
                    else
                        noise_enhancement_tmp = noise_enhancement_tmp + sum(abs(inv_temp).^2,2);
                    end
                    
                    Hg = Hg + inv_temp*Hsc;
                end
                Hg = Hg/Nsymbols;
        end
        

        noise_enhancement = zeros(sum(M),Nsymbols);
        
        for ii=1:nLayers
            
            if LTE_params.DFT_spreading_off
                noise_enhancement(sum(M(1:ii-1)) + (1:M(ii)),:) = kron(ones(M(ii),1),noise_enhancement_tmp(ii,:));
            else
%                 noise_enhancement(sum(M(1:ii-1)) + (1:M(ii)),:) = kron(ones(M(ii),LTE_params.Ntot),noise_enhancement_tmp(ii));              
                noise_enhancement(sum(M(1:ii-1)) + (1:M(ii)),:) = kron(ones(M(ii),Nsymbols),noise_enhancement_tmp(ii));
            end
            
        end
        
        %===========
        %sigma_n2 in noise enhancement included!        
        %===========
        
        
        if LTE_params.DFT_spreading_off
            rx_layer_x = transpose(rx_layer_x);
            noise_enhancement = sigma_n2*noise_enhancement;
        else
            rx_layer_x = transpose(ifft(rx_layer_x,Nsymbols))*sqrt(Nsymbols);
            noise_enhancement = sigma_n2/Nsymbols*noise_enhancement;
        end        
        
        
        LLR_SD(:,index+1:index+Nsymbols)= LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
        index = index+Nsymbols;
        
        
    end
end