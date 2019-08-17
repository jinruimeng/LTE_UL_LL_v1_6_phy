function [LTE_params,BS_output]=LTE_tx_PHICH(LTE_params,subframe_corr,BS,UE,b_,BS_output)
%% PHICH processing 
% TS 36.211 V8.9.0, Section 6.9
% TS 36.212 V8.8.0, Section 5.35
% Author: Petr Kej�k, xkejik00@stud.feec.vutbr.cz
% 2010
% 
% input :   Nrb       ... [1 x 1]double number of resource blocks
%           Nsc       ... [1 x 1]double number - resource block size
%           Nsub      ... [1 x 1]double number of OFDM symbols in subframe
%           Ns        ... [1 x 1]double number of OFDM symbols in slot
%           NoData    ... reserved space for synchronization signal
 
%% Warnings, comments, possible changes
% The PHICH shall be transmitted on the same set of antenna ports as the PBCH.
% 4 subframes latter % PDCCH - format 0 with DCI
N_g = 1; % TS 36.101 V8.11.0 Section 8.5.1 (test parameters for PHICH) N_g = {1/6, 1/2, 1, 2}
 
%% Number of PHICH groups
% TS 36.211 V8.9.0 Section 6.9
AtPorts = LTE_params.BS_config.nAtPort;
switch LTE_params.CyclicPrefix         
     case 'normal'
         N_g_phich = ceil(N_g*(LTE_params.Nrb/8));
         N_SF = 4;
         w_phich = [ +1,  +1,  +1,  +1;
                     +1,  -1,  +1,  -1;
                     +1,  +1,  -1,  -1;
                     +1,  -1,  -1,  +1;
                    +1j, +1j, +1j, +1j;
                    +1j, -1j, +1j, -1j;
                    +1j, +1j, -1j, -1j;
                    +1j, -1j, -1j, +1j];
     case 'extended'
         N_g_phich = 2*ceil(N_g*(LTE_params.Nrb/8));
         N_SF = 2;
         w_phich = [ +1,  +1;
                     +1,  -1;
                    +1j, +1j;
                    +1j, -1j];
     otherwise
         error('Wrong cyclic prefix (LTE-tx-PHICH)')
end
 
%% Number of PHICH processes
N_phich_proc = LTE_params.HARQ_processes; 
% random generation of ACK(70%)/NACK(30%) 
HIgroup = randsrc(N_phich_proc,1,[0 1; 0.3 0.7]);

BS_output(b_).PHICH_HI_tx = HIgroup(1,1);
 
if N_phich_proc > 8
    error('PHICH groups error - not implemented yet (LTE-common-gen-PHICH)')
end
 
N_phich_proc = N_phich_proc/N_g_phich;
HI_y_tx = zeros(AtPorts,12,N_g_phich);
 
for i0 = 1:N_g_phich % number of PHICH groups
  HI_y_tx_el = zeros(AtPorts,12,N_phich_proc);

  for n_g_PH = 1:N_phich_proc % number of PHICH processes in one group

    %% Channel coding
    % TS 36.212 V8.8.0 Section 5.3.5
    HI_b = [HIgroup(n_g_PH,1) HIgroup(n_g_PH,1) HIgroup(n_g_PH,1)];
 
    %% Modulation
    % TS 36.211 V8.9.0 Section 6.9.1
    HI_z = zeros(1,3);
    nibble = HI_b(1,:);
    HI_z(1,:) = LTE_params.SymbolAlphabet{LTE_params.PHICH_common_param.modulation_order}(nibble+1).';
 
    %% Spreading and scrambling
    % TS 36.211 V8.9.0 Section 6.9.1
    c_init = (floor((2*subframe_corr-2)/2)+1)*(2*BS.NIDcell+1)*(2^9)+BS.NIDcell; % initialization value
    pn_seq = LTE_common_gen_gold_sequence_CCH(3*N_SF,c_init); % pn generation
 
    HI_d = zeros(1,3*N_SF);
    for i2 = 1:3*N_SF % index i            -1      +1       --- zero excluded ---       -1       +1
        HI_d(1,i2) = w_phich(n_g_PH,mod(i2-1,N_SF)+1)*(1-2*pn_seq(i2))*HI_z(1,floor((i2-1)/N_SF)+1); 
    end
    % check if the modulation and spreading process performs correctly 
    HI_d_help = [HI_z(1,1)*w_phich(n_g_PH,:) HI_z(1,2)*w_phich(n_g_PH,:) HI_z(1,3)*w_phich(n_g_PH,:)];
    HI_d_spread_help = (1-2*pn_seq).*HI_d_help;
    if sum(HI_d ~= HI_d_spread_help) > 0
        error('Wrong PHICH scrambling or spreading (LTE-tx-PHICH)' )
    end
     
    %% Resource group alignment
    % TS 36.211 V8.9.0 Section 6.9.2
    HI_d_al = zeros(1,3*N_SF);
    switch LTE_params.CyclicPrefix
        case 'normal'
            HI_d_al = HI_d ;
        case 'extended'
            if mod(n_g_PH-1,2) == 0
                 for i3 = 1:(3*N_SF/2)-1 +1
                     HI_d2(1,(i3-1)*4+1:(i3-1)*4+4) = [HI_d(1,2*i3-1) HI_d(1,2*i3) 0 0];
                 end
             elseif mod(n_g_PH-1,2) == 1
                 for i3 = 1:(3*N_SF/2)-1 +1
                     HI_d2(1,(i3-1)*4+1:(i3-1)*4+4) = [0 0 HI_d(1,2*i3-1) HI_d(1,2*i3)];
                 end
             end           
             HI_d_al = HI_d2;
        otherwise
            error('Wrong cyclic prefix (LTE-tx-PHICH)')
    end
 
    %% Layer mapping 
    % TS 36.211 V8.9.0 Section 6.9.2
    nLayers = LTE_params.BS_config.nAtPort; % Section 6.3.3.3
    switch nLayers % transmission mode 
        case 1  % single antenna transmission, TS 36.211, Section 6.3.3.1
            layer_xHI = HI_d_al;
        case 2  % transmit diversity, TS 36.211, Section 6.3.3.3
            layer_xHI(1,:) = HI_d_al(1:2:end);
            layer_xHI(2,:) = HI_d_al(2:2:end);
        case 4
            if (mod(length(HI_d),4) ~= 0)
                HI_d_al = [HI_d_al,0,0];
            end  
            layer_xHI(1,:) = HI_d_al(1:4:end);
            layer_xHI(2,:) = HI_d_al(2:4:end);
            layer_xHI(3,:) = HI_d_al(3:4:end);
            layer_xHI(4,:) = HI_d_al(4:4:end);
        otherwise
            error('Number of layers not supported (LTE-tx-PHICH)');
    end
     
    %% Precoding
    % TS 36.211 V8.9.0 Section 6.9.2 (6.3.4.1 and 6.3.4.3)
    AtPorts = LTE_params.BS_config.nAtPort;
    nLayers = 1; % 1 layer for PHICH 
    tx_mode = UE.mode; 
    codebook_index = 0; % no meaning for PHICH - spatial multiplexing is not used for PHICH
    RI = 0;             % no meaning for PHICH
    CDD = 0;            % no meaning for PHICH 
    rb_numbers = 0;     % no meaning for PHICH
    indices = 0;        % no meaning for PHICH
    slots = 0;          % no meaning for PHICH
 
    if UE.mode > 2; % spatial multiplexing is not used for PCFICH
        error('Spatial multiplexing is not used for PHICH (LTE-tx-PHICH)')
    end
     
    if BS.nAtPort == 4 % special precoding is used in this case 
        Z_phich1 = [1,  0, 0,  0,  1j,  0,   0,  0; 
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0, -1, 0,  0,   0, 1j,   0,  0;
                    0,  0, 0,  0,   0,  0,   0,  0; 
                    0,  1, 0,  0,   0, 1j,   0,  0;
                    0,  0, 0,  0,   0,  0,   0,  0; 
                    1,  0, 0,  0, -1j,  0,   0,  0;
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0,  0, 1,  0,   0,  0,  1j,  0;
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0,  0, 0, -1,   0,  0,   0, 1j;
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0,  0, 0,  1,   0,  0,   0, 1j;
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0,  0, 1,  0,   0,  0, -1j,  0;
                    0,  0, 0,  0,   0,  0,   0,  0];

        Z_phich2 = [0,  0, 0,  0,   0,  0,   0,  0; 
                    1,  0, 0,  0,  1j,  0,   0,  0;
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0, -1, 0,  0,   0, 1j,   0,  0; 
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0,  1, 0,  0,   0, 1j,   0,  0; 
                    0,  0, 0,  0,   0,  0,   0,  0;
                    1,  0, 0,  0, -1j,  0,   0,  0;
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0,  0, 1,  0,   0,  0,  1j,  0;
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0,  0, 0, -1,   0,  0,   0, 1j;
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0,  0, 0,  1,   0,  0,   0, 1j;
                    0,  0, 0,  0,   0,  0,   0,  0;
                    0,  0, 1,  0,   0,  0, -1j,  0];
                    
        c = length(layer_xHI(1,:));
        precode_yHI = zeros(4,4*c); % NOTE: cleanup
        for ii = 1:3 % index i - section 6.9.2
            if (mod((ii-1)+(n_g_PH-1),2)==0 && strcmp(LTE_params.CyclicPrefix,'normal')) || (mod((ii-1)+floor((n_g_PH-1)/2),2)==0 && strcmp(LTE_params.CyclicPrefix,'extended'))
                Z = Z_phich1;
            else
                Z = Z_phich2;
            end
            precode_y_te = zeros(4,4);
            X = 1/sqrt(2)*Z*[real(layer_xHI(1,ii));real(layer_xHI(2,ii));real(layer_xHI(3,ii));real(layer_xHI(4,ii));...
                             imag(layer_xHI(1,ii));imag(layer_xHI(2,ii));imag(layer_xHI(3,ii));imag(layer_xHI(4,ii))];
            precode_y_te(1,1) = X(1,:); 
            precode_y_te(2,1) = X(2,:);
            precode_y_te(3,1) = X(3,:); 
            precode_y_te(4,1) = X(4,:);
            precode_y_te(1,2) = X(5,:); 
            precode_y_te(2,2) = X(6,:);
            precode_y_te(3,2) = X(7,:); 
            precode_y_te(4,2) = X(8,:);
            precode_y_te(1,3) = X(9,:); 
            precode_y_te(2,3) = X(10,:);
            precode_y_te(3,3) = X(11,:); 
            precode_y_te(4,3) = X(12,:);
            precode_y_te(1,4) = X(13,:); 
            precode_y_te(2,4) = X(14,:);
            precode_y_te(3,4) = X(15,:); 
            precode_y_te(4,4) = X(16,:);
            
            precode_yHI(:,4*(ii-1)+1:4*(ii-1)+4) = precode_y_te; 
        end
   else % standard precoding is used TS 36.211 Section 6.3.4.1 or 6.3.4.3
        [precode_yHI PRE] = LTE_precoding(tx_mode,layer_xHI,BS.nAtPort,codebook_index,nLayers,indices,LTE_params.Ntot,CDD,LTE_params,slots);
        [precode_yHI2,UE,D,W,U] = LTE_tx_precoding(LTE_params,layer_xHI,UE,AtPorts,codebook_index,RI,CDD,rb_numbers);
        if sum(precode_yHI2 ~= precode_yHI) > 0
            error('Wrong PHICH scrambling or spreading (LTE-tx-PHICH)' )
        end
    end
    HI_y_tx_el(:,:,n_g_PH) = precode_yHI; % individual PHICH processes in one PHICH group
%     
  end
% TS 36.211 V8.9.0 Section 6.9.3
HI_y_tx(:,:,i0) = sum(HI_y_tx_el,3);

end

LTE_params.PHICH(b_,subframe_corr).HI_y_tx = HI_y_tx;
