function [LLR_SD,M] = LTE_detecting_ICI(BS_output, UE_output, symbols_per_UE, H_est_per_UE, LTE_params, sigma_n2, UE_MCS_and_scheduling_info)
% Detection.
% Author: Lukas Nagel lnagel@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

nBS = LTE_params.nBS;
nUE = LTE_params.nUE;
nRX = LTE_params.BS_config.nRX;

scheduled_RBs_UE = cell(1, nBS*nUE);
detected_symbols_UE = cell(1, nBS*nUE);
M = cell(1, nBS*nUE);
symbol_alphabet = cell(1, nBS*nUE);
bittable = cell(1, nBS*nUE);
index = cell(1, nBS*nUE);
LLR_SD = cell(1, nBS*nUE);

for bb = 1:nBS
    for ul = 1:nUE
        uu = utils.localToGlobalUser( LTE_params.connection_table, bb, ul );
        scheduled_RBs_UE{uu} = UE_MCS_and_scheduling_info(bb,ul).UE_mapping;
        M{uu} = UE_MCS_and_scheduling_info(bb,ul).CQI_params.modulation_order;
        symbol_alphabet{uu} = LTE_params.SymbolAlphabet{M{uu}}.';
        bittable{uu} = LTE_params.bittable{M{uu}};
        
        index{uu} = 0;
        LLR_SD{uu} = zeros(M{uu},size(symbols_per_UE{uu},1)*size(symbols_per_UE{uu},2));
    end
end

transmit_power_factor = 10.^(-LTE_params.pathloss_matrix./10);


tmp_used = UE_output(1).UE_genie.CH_mapping(1,:);
tmp_sel = 1:length(tmp_used);
used_n = tmp_sel(tmp_used);
Ns = length(used_n);

%% TODO :)
switch LTE_params.UE_config.mode 
    case 1
        if isequal(LTE_params.ChanMod_config.filtering, 'FastFading')

            nl = 0; % second relative time index
            for nn = used_n  
                nl = nl + 1;
                if nn <= 7
                    s_slot = 1;
                else
                    s_slot = 2;
                end
                
                for uu = 1:nBS*nUE % for all users
                    scheduled_subc = sum(scheduled_RBs_UE{uu}(:,1))*12; % assume schedule is valid for two slots..
                    detected_symbols = zeros(scheduled_subc, 1);
                    noise_enhancement_h = zeros(scheduled_subc, 1);
                    
                    
                    if scheduled_subc > 0 % only detect if the user is scheduled
                        base_k = (find(scheduled_RBs_UE{uu}(:,s_slot)~=0,1,'first')-1)*12;
                        
                        for kk = 1:scheduled_subc 
                            

                            H1 = squeeze(H_est_per_UE{uu}{uu}(base_k+kk,nn,:)); % nRX x 1 

                            H2 = zeros(nRX,nRX); % nRX x nRX replace constants

                            for u_ = 1:nBS*nUE
                                if u_ ~= uu
                                    if scheduled_RBs_UE{u_}(floor((base_k+kk-1)/12)+1, s_slot) % check rb calculation!
                                        % build H matrix
                                        Ht = squeeze(H_est_per_UE{uu}{u_}(base_k+kk,nn,:));

                                        % find the right factor
                                        bb = find(LTE_params.connection_table(:,uu)==1);
                                        alpha = transmit_power_factor(bb, u_);

                                        H2 = H2 + alpha*Ht*Ht'; 
                                    end
                                end
                            end
                            
                            G = H1'* inv(H1*H1' + H2 + sigma_n2 * eye(nRX) );
                            
                            noise_enhancement_h(kk) = mean(abs(G).^2);
                            
                            yu1 = squeeze(symbols_per_UE{uu}(:,nl,:));
                           
                            yu = yu1(kk,:);

                            xu = G * yu.';

                            detected_symbols(kk) = xu;


                        end
                        
                        N_sym = length(detected_symbols);
                        if LTE_params.DFT_spreading_off % OFDM
                            noise_enhancement = ones(M{uu},1)*sigma_n2*noise_enhancement_h.';
                        else % SCFDM
                            detected_symbols = ifft(detected_symbols,N_sym)*sqrt(N_sym); 
                            noise_enhancement = ones(M{uu},1)*sigma_n2*noise_enhancement_h.';
                        end

                        
                        Nsymbols = size(symbols_per_UE{uu},1);
                        LLR_SD{uu}(:,index{uu}+1:index{uu}+Nsymbols) = LTE_demapper(detected_symbols.',symbol_alphabet{uu},bittable{uu},1,M{uu},1,noise_enhancement);%
                        index{uu} = index{uu} + Nsymbols;
                        
                    end % if scheduled
                end
            end
            
            
        else % BLOCKFADING
            
            nl = 0; % second relative time index
            nl = nl + 1;

            G = cell(nBS*nUE, 1);
            noise_enhancement_h = cell(nBS*nUE, 1);


            for uu = 1:nBS*nUE % for all users
                scheduled_subc = sum(scheduled_RBs_UE{uu}(:,1))*12; % assume schedule is valid for two slots..
                detected_symbols = zeros(scheduled_subc, 1);
                noise_enhancement_h{uu} = zeros(scheduled_subc, 1);
                G{uu} = zeros(scheduled_subc, nRX);


                if scheduled_subc > 0 % only detect if the user is scheduled
                    base_k = (find(scheduled_RBs_UE{uu}(:,1)~=0,1,'first')-1)*12;

                    for kk = 1:scheduled_subc 

                        H1 = squeeze(H_est_per_UE{uu}{uu}(base_k+kk,1,:)); % nRX x 1 

                        H2 = zeros(nRX,nRX); % nRX x nRX replace constants

                        for u_ = 1:nBS*nUE
                            if u_ ~= uu
                                if scheduled_RBs_UE{u_}(floor((base_k+kk-1)/12)+1, 1) 
                                    % build H matrix
                                    Ht = squeeze(H_est_per_UE{uu}{u_}(base_k+kk,1,:));

                                    % find the right factor
                                    bb = find(LTE_params.connection_table(:,uu)==1);
                                    alpha = transmit_power_factor(bb, u_);

                                    H2 = H2 + alpha*Ht*Ht'; 
                                end
                            end
                        end
                        
                        G{uu}(kk, :) = H1'* inv(H1*H1' + H2 + sigma_n2 * eye(nRX) );

                        noise_enhancement_h{uu}(kk) = mean(abs(G{uu}(kk, :)).^2);
                    end
                end
            end
            
            for uu = 1:nBS*nUE
                scheduled_subc = sum(scheduled_RBs_UE{uu}(:,1))*12; 
                nl = 0;
                for nn = used_n
                    detected_symbols = zeros(scheduled_subc, 1);
                    nl = nl + 1;
                    if scheduled_subc > 0 % only detect if the user is scheduled
                        for kk = 1:scheduled_subc 
                            yu1 = squeeze(symbols_per_UE{uu}(:,nl,:));

                            detected_symbols(kk) = G{uu}(kk,:) * yu1(kk,:).';
                        end


                        N_sym = length(detected_symbols);
                        if LTE_params.DFT_spreading_off % OFDM
                            noise_enhancement{uu} = ones(M{uu},1)*sigma_n2*noise_enhancement_h{uu}.';
                        else % SCFDM
                            detected_symbols = ifft(detected_symbols,N_sym)*sqrt(N_sym); 
                            noise_enhancement{uu} = ones(M{uu},1)*sigma_n2*noise_enhancement_h{uu}.';
                        end


                        Nsymbols = size(symbols_per_UE{uu},1);
                        LLR_SD{uu}(:,index{uu}+1:index{uu}+Nsymbols) = LTE_demapper(detected_symbols.',symbol_alphabet{uu},bittable{uu},1,M{uu},1,noise_enhancement{uu});%
                        index{uu} = index{uu} + Nsymbols;

                    end % if scheduled
                end
            end   
        end
        
        
    case 4
        error('currently only SISO is implemented for ICI');
        
    otherwise
        error('This transmission mode is not standardized in Uplink LTE');
end

end




