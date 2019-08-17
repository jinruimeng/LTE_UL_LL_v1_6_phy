function LTE_UL_TX(LTE_params, UE, BS, subframe_i, UE_output, BS_output, uu, bb, srs_subframe)
% physical layer functions for uplink
% based on LTE_TX ver. 1.6
% partial code inserted
% shared functions from downlink used
% PUSCH channel
% AtPort is not included, only MU-MIMO is defined
global DEBUG_LEVEL

% find the right basestation
attached_BS = BS(logical(LTE_params.connection_table(:,uu))');
local_uu = utils.globalToLocalUser(LTE_params.connection_table, bb, uu);

nCodewords   = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.nCodewords;
nLayers      = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.nLayers;
assigned_RBs = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.assigned_RBs;
UE_mapping   = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.UE_mapping;
tx_mode      = LTE_params.UE_config.mode;
DM_SRS_allocation = LTE_params.DM_SRS_allocation;


codeBookIndex = 0;
if LTE_params.UE_config.mode == 4
%     codeBookIndex = mode(mode(UE_output.PMI));  % pick the most freqent value of PMIs
    codeBookIndex = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.PMI;
end

Ns = LTE_params.Ns;     % number of symbols in one slot
Nsc = LTE_params.Nsc;   % number of subcarriers
Nrb = LTE_params.Nrb;   % number of resource blocks
Ntot = LTE_params.Ntot; % total number of subcarriers
Nsub = LTE_params.Nsub; % number of symbols in one subframe
MscPUSCH = assigned_RBs*Nsc/max(sum(UE_mapping,2));     % assigned_RBs*Nsc is the number of subcarrier for TWO RBs

% define where DMRS and SRS are in a subframe
switch LTE_params.CyclicPrefix
    case 'normal'   % symbols 0 to 13
        DMRS_index = 4;     % 3rd symblo per slot
        SRS_index = 14;     % 13th symbol in a subframe
    case 'extended' % smybols 0 to 11
        DMRS_index = 3;     % 2nd symbol per slot
        SRS_index = 12;     % 11th symbol per subframe
    otherwise
        error('CP not supported');
end

% reference symbol mapping
tx = zeros(nLayers,Ntot,Nsub);
tx_temp = zeros(Ntot,Nsub);
DM_tmp = false(Nsc,Ns);
DM_tmp(:,DMRS_index) = true;
DM_mapping = logical(kron(UE_mapping,DM_tmp));
data_mapping = logical(kron(UE_mapping,true(Nsc,Ns)));
UE_allocation = xor(data_mapping,DM_mapping);

if srs_subframe
    UE_allocation(:,SRS_index) = false;
end

%% Calculate the correct number of subframe from 1 - 10
subframe_corr = mod(subframe_i,10);
if(subframe_corr == 0)
    subframe_corr = 10;
end

%% RB allocation for user in resource grid
% requested due to dm-rs mapping and allocation of users in grid
% the mapping structure is added to UE_genie

% insertion of the user data indices for receiver
BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.slot_indices = [];
BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.freq_indices = [];
RB_indices = find(UE_mapping == 1);
for ww = 1:assigned_RBs
    zero_temp = zeros(size(UE_mapping));
    zero_temp(RB_indices(ww)) = 1;
    
    freq_tmp = mod(find(kron(zero_temp.*UE_mapping, ones(Nsc,Ns)).*(UE_allocation)), LTE_params.Ntot);
    BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.freq_indices = [BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.freq_indices;freq_tmp];
    BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.slot_indices = [BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.slot_indices;(floor((RB_indices(ww)-1)/LTE_params.Nrb)+1)*ones(size(freq_tmp))];
    
end
BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.freq_indices(~BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.freq_indices) = LTE_params.Ntot; % where the mod operation delivers a zero there should be Ntot

UE_output.HARQ_process = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.HARQ_process_id;

%% Generate data
layer_x = [];
for cw = 1:nCodewords

    
    % Get the number of data and coded bits from the scheduler
    N_coded_bits = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.N_coded_bits(cw);
    N_data_bits  = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.N_data_bits(cw);
    
    tx_rv_idx = attached_BS.UE_specific(local_uu).current_HARQ_process(cw).rv_idx;
    % Generate data bits
    if tx_rv_idx==0
        if LTE_params.simulate_with_all_zero_sequences
            tx_data_bits = false(1,N_data_bits);
        else
            tx_data_bits = logical(randi(LTE_params.data_RandStream,[0,1],[1,N_data_bits]));
        end
        attached_BS.UE_specific(local_uu).current_HARQ_process(cw).HARQ_tx_buffer = tx_data_bits;
    else
        % Retransmission, use previously stored bits in the current HARQ process
        tx_data_bits = attached_BS.UE_specific(local_uu).current_HARQ_process(cw).HARQ_tx_buffer;
    end
    %% Coding of the bits, as of TS 36.212
    BS_output.UE_signaling_UL(bb, local_uu).turbo_rate_matcher_UL(cw).G = N_coded_bits;
    % How many bits the we are allowed to transmit. Decided by the scheduler(uu)
    % According to TS 36.212, subclause 5.1.4.1.2. N_l is:
    % (equal to the number of layers a transport block is mapped onto)
    % equal to one for blocks mapped onto one layer,
    % equal to two for blocks mapped onto two or four layers
    switch nLayers
        case 1
            BS_output.UE_signaling_UL(bb,local_uu).turbo_rate_matcher_UL(cw).N_l = 1;
        case 2
            if nCodewords == 1
                BS_output.UE_signaling_UL(bb,local_uu).turbo_rate_matcher_UL(cw).N_l = 1;
            else
                BS_output.UE_signaling_UL(bb,local_uu).turbo_rate_matcher_UL(cw).N_l = 2;
            end
        case 3
            if nCodewords == 1
                BS_output.UE_signaling_UL(bb,local_uu).turbo_rate_matcher_UL(cw).N_l = 1;
            else
                if (cw == 1)
                    BS_output.UE_signaling_UL(bb,local_uu).turbo_rate_matcher_UL(cw).N_l = 1;
                else
                    BS_output.UE_signaling_UL(bb,local_uu).turbo_rate_matcher_UL(cw).N_l = 2;
                end
            end
        case 4
            BS_output.UE_signaling_UL(bb,local_uu).turbo_rate_matcher_UL(cw).N_l = 2;
    end
    
    BS_output.UE_signaling_UL(bb,local_uu).turbo_rate_matcher_UL(cw).rv_idx = tx_rv_idx; % redundancy version index
    
    %% Channel coding
    tx_coded_bits = logical(LTE_UL_tx_ULSCH_encode(LTE_params,tx_data_bits,BS_output.UE_signaling_UL(bb,local_uu),UE_output.UE_genie,UE,cw));
    BS_output.genie(bb,local_uu).sent_bits{cw} = tx_coded_bits;
    
    %% Bit scrambling
    % like in downlink, placeholders for ARQ-ACK and RI not yet implemented
    % see TS 36.211 V11.4.0 (2013-09) Section 5.3.1
    
    tx_scrambled_bits = LTE_common_scrambling(tx_coded_bits,attached_BS.NIDcell,uu,subframe_corr,cw,'scramble');
    
    %% Symbol mapping
    M = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.CQI_params(cw).modulation_order;
    powers = 2.^(M-1:-1:0);
    tx_scrambled_bits = reshape(tx_scrambled_bits,M,[]);
    symbols_int = powers*tx_scrambled_bits;
    
    tx_user_symbols = LTE_params.SymbolAlphabet{M}(symbols_int+1).';
    
    %% Layer mapping
    layer_x = LTE_common_layer_mapping(tx_mode,tx_user_symbols,nLayers,layer_x,nCodewords,cw,attached_BS.nAtPort);
end


for ii=1:nLayers
    %% Transform Precoding, defined TS 36.211 V11.4.0 section 5.3.3
    % DFT, SC-FDMA symbol is one column
    % number of subcarriers from LTE_params
    % number of symbols from scheduling
    
    if LTE_params.DFT_spreading_off
        z_ul = reshape(layer_x(ii,:),MscPUSCH,[]);   % to get OFDM no DFT spread!
    else
        z_ul = fft(reshape(layer_x(ii,:),MscPUSCH,[]),MscPUSCH)/sqrt(MscPUSCH); % check this for the amplitudes
    end
    
    
    
    %% Resource element Mapping, before iFFT
    % user is mapped independently
    % SRS, generation is made separately
    % located in 13th OFDM symbol
    % DMRS in 3rd symbol in a slot (0  1  2  3  4  5  6
    %                               7  8  9 10 11 12 13)
    
    dmrs = LTE_UL_dm_srs_mapping(LTE_params, BS(logical(LTE_params.connection_table(:,uu))'), BS_output, UE_output, UE, bb, uu, 'dmrs', ii); % generate the DMRS
    
    tx_temp(DM_mapping(:,1:Nsub-1)) = dmrs(:,DM_SRS_allocation(1:2));                       % multiplex DMRS into a 'Ntot x Nsub' grid  (all subcarriers and one subframe)
    tx_temp(UE_allocation) = z_ul;                                                          % multiplex precoder output into the same grid
    tx(ii,:,:) = tx_temp;                                                                   % do it for every layer
end

% UE_output.UE_genie.y_tx_assembled = tx;         % stored for receiver
UE_output.UE_genie.CH_mapping = UE_allocation;  % user allocation of RE
UE_output.UE_genie.DM_mapping = DM_mapping;     % DMRS allocation

%% precoding, defined TS 36.211 v11.1.0 section 5.3.3A

W = LTE_UL_Codebook(nLayers,codeBookIndex,UE.nTX);
UE_output.LayerMapping = W;

sizes = size(tx);
sizes(1) = size(W,1);
tx = reshape(tx,size(W,2),[]);
tx = W*tx;
tx = reshape(tx,sizes);

% multiplex SRS after precoding
if DM_SRS_allocation(3)
    srs = transpose(LTE_UL_dm_srs_mapping(LTE_params, BS(logical(LTE_params.connection_table(:,uu))'), BS_output, UE_output, UE, bb, uu, 'srs', ii));    % generate SRS
    tx(:,:,Nsub) = srs;                                                                                 % multiplex them
end

UE_output.UE_genie.y_tx_assembled = tx;

for ii=1:UE.nTX
    tx_antenna = squeeze(tx(ii,:,:));
    
    %% SC-FDMA signal generation, defined TS 36.211 v. 8.9.0
    % Padding, iFFT, Cyclic prefix insertion
    
    % zero padding
    tx_padded = ([tx_antenna; zeros(LTE_params.Nfft - LTE_params.Ntot,Nsub)]);
    
    % cyclic shift to baseband
    tx_shifted = (circshift(tx_padded,-(LTE_params.Ntot)/2));
    
    % IFFT
    tx_symb = ifft(tx_shifted)*sqrt(LTE_params.Nfft);
    
    % CP insertion
    UE_output.y_tx(:,ii) = [    tx_symb(LTE_params.Index_TxCyclicPrefix{1},1);...
        reshape(tx_symb(LTE_params.Index_TxCyclicPrefix{2},2:Ns),[],1);...
        tx_symb(LTE_params.Index_TxCyclicPrefix{1},Ns+1);...
        reshape(tx_symb(LTE_params.Index_TxCyclicPrefix{2},Ns+2:end),[],1) ];
    
    %% PAPR calculation
    % according (6) of "An Overview: Peak-to-Average Power RatioReduction
    % Techniques for OFDM Signals"
    
    L=4;    
    
    %  exclude DMRS+SRS
    datapos=UE_output.UE_genie.CH_mapping;
    tx_papr = reshape(tx_antenna(datapos),sum(datapos(:,1)),[]);
    
    % zero padding + cyclic shift
    tx_papr = ([tx_papr; zeros(LTE_params.Nfft - LTE_params.Ntot,size(tx_papr,2))]);
    tx_papr = (circshift(tx_papr,-(LTE_params.Ntot)/2));

    %one PAPR value per OFDM symbol 
    y_papr(:,:,ii) = ifft(tx_papr, LTE_params.Nfft*L)*sqrt(LTE_params.Nfft); 
end

UE_output.UE_genie.y_tx = UE_output.y_tx;
                                   
UE_output.papr = transpose(10*log10(max(max(abs(y_papr).^2,[],3))/ mean(mean(mean(abs(y_papr).^2)))));      %PAPR vector



% function LTE_UL_TX_plot(input_var,x_title,y_title,main_title)
% figure
% image((abs(squeeze(input_var(1,:,:))))*100);
% grid on
% set(gca,'xcolor','w','ycolor','w');
% set(gca,'Fontsize',12,'Linewidth',1.5,'XTick',[0.5:7:14.5],'YTick',[0.5:12:128.5],'YTicklabel',[1,13,25,37,49,61,73,85,97,109,121]-1,'XTicklabel',[1,8,15]-1);
% xlabel(x_title);
% ylabel(y_title);
% title(main_title)
