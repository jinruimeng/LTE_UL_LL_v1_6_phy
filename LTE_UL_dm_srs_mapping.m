function [rs] = LTE_UL_dm_srs_mapping(LTE_params, BS, BS_output, UE_output, UE, bb, uu, type, layer)

% generation of basic sequences for DM-RS and SRS signal
% DM-RS signal generation 
% SRS signal generation 
% mapping into subframe
% according to TS 136 211 V11.1.0 section 5.5
% 2010 Jan Prokopec
%[bb, uu]
local_uu = utils.globalToLocalUser(LTE_params.connection_table, bb, uu);

UE_mapping = BS_output.UE_signaling_UL(bb,local_uu).MCS_and_scheduling_UL.UE_mapping;
m = sum(UE_mapping,1);
m = m(1);

% Require user mapping, number of resource blocks

% m = 1;  % 1=<m=<NrbMAXul==110
Nsc = LTE_params.Nsc;
Nrb = LTE_params.Nrb;

% n = 1;                    % 0 =< n =< MscRS-1
%% higher layer configuration for DMRS
cyclicShift = LTE_params.cyclicShift;           % {0 to 7}
% CSFiDCIFormat = LTE_params.CSFiDCIFormat;       % Cyclic Shif Field in DCI format 0 [3] (page 26)
deltass = LTE_params.deltass;                   % {0..29}
groupHopping = LTE_params.groupHopping;         % enable sequence group hopping
sequenceHopping = LTE_params.sequenceHopping;   % enable sequence hopping
nID_PUSCH = LTE_params.nID_PUSCH;               % can be provided, otherwise NAN        
nID_PUCCH = LTE_params.nID_PUCCH;               % can be provided, otherwise NAN 

%% higher layer configuration for SRS
Csrs = LTE_params.Csrs;                     % srs-BandwidthConfig {0..7}
srsMaxUpPTS = LTE_params.srsMaxUpPTS;       % {true, false}
frameStructure = LTE_params.frameStructure; % frame structure type {1,2} = {FDD,TDD}
Nsp = LTE_params.Nsp;                       % number of uplink to downlink switching points {1,2}
% CSFiDCIFormat   = LTE_params.CSFiDCIFormat;

% UE specific
CSFiDCIFormat   = UE.CSFiDCIFormat;
Tsrs            = UE.Tsrs;                  % srs configuration index (0..1023)
Toffset         = UE.Toffset;
nCS_SRS         = UE.nCS_SRS;               % cyclicShift {0..7} 
Bsrs            = UE.Bsrs;                  % srs-Bandwidth {0..3}
b_hop           = UE.b_hop;                 % srs-HoppingBandwidth {0..3}
nRRC            = UE.nRRC;                  % freqDomainPosition {0..31}
kTC             = UE.kTC;                   % transmissionComb {0..1}

transmission = 'PUSCH';

subframe = mod(BS.clock.current_TTI-1,LTE_params.Nsfr);
nf = floor((BS.clock.current_TTI-1)/LTE_params.Nsfr);
ns = LTE_params.Nslot*subframe:LTE_params.Nslot*(subframe+1)-1;

%% Generating the base sequence

fi1=[ -1  1  3 -3  3  3  1  1  3  1 -3  3;
    1  1  3  3  3 -1  1 -3 -3  1 -3  3;
    1  1 -3 -3 -3 -1 -3 -3  1 -3  1 -1;
    -1  1  1  1  1 -1 -3 -3  1 -3  3 -1;
    -1  3  1 -1  1 -1 -3 -1  1 -1  1  3;
    1 -3  3 -1 -1  1  1 -1 -1  3 -3  1;
    -1  3 -3 -3 -3  3  1 -1  3  3 -3  1;
    -3 -1 -1 -1 -1  1 -3  3 -1  1 -3  3;
    1 -3  3  1 -1 -1 -1  1  1  3 -1  1;
    1 -3 -1  3  3 -1 -3  1  1  1  1  1;
    -1  3 -1  1  1 -3 -3 -1 -3 -3  3 -1;
    3  1 -1 -1  3  3 -3  1  3  1  3  3;
    1 -3  1  1 -3  1  1  1 -3 -3 -3  1;
    3  3 -3  3 -3  1  1  3 -1 -3  3  3;
    -3  1 -1 -3 -1  3  1  3  3  3 -1  1;
    3 -1  1 -3 -1 -1  1  1  3  1 -1 -3;
    1  3  1 -1  1  3  3  3 -1 -1  3 -1;
    -3  1  1  3 -3  3 -3 -3  3  1  3 -1;
    -3  3  1  1 -3  1 -3 -3 -1 -1  1 -3;
    -1  3  1  3  1 -1 -1  3 -3 -1 -3 -1;
    -1 -3  1  1  1  1  3  1 -1  1 -3 -1;
    -1  3 -1  1 -3 -3 -3 -3 -3  1 -1 -3;
    1  1 -3 -3 -3 -3 -1  3 -3  1 -3  3;
    1  1 -1 -3 -1 -3  1 -1  1  3 -1  1;
    1  1  3  1  3  3 -1  1 -1 -3 -3  1;
    1 -3  3  3  1  3  3  1 -3 -1 -1  3;
    1  3 -3 -3  3 -3  1 -1 -1  3 -1 -3;
    -3 -1 -3 -1 -3  3  1 -1  1  3 -3 -3;
    -1  3 -3  3 -1  3  3 -3  3  3 -1 -1;
    3 -3 -3 -1 -1 -3 -1  3 -3  3  1 -1];

fi2=[ -1  3  1 -3  3 -1  1  3 -3  3  1  3 -3  3  1  1 -1  1  3 -3  3 -3 -1 -3;
    -3  3 -3 -3 -3  1 -3 -3  3 -1  1  1  1  3  1 -1  3 -3 -3  1  3  1  1 -3;
    3 -1  3  3  1  1 -3  3  3  3  3  1 -1  3 -1  1  1 -1 -3 -1 -1  1  3  3;
    -1 -3  1  1  3 -3  1  1 -3 -1 -1  1  3  1  3  1 -1  3  1  1 -3 -1 -3 -1;
    -1 -1 -1 -3 -3 -1  1  1  3  3 -1  3 -1  1 -1 -3  1 -1 -3 -3  1 -3 -1 -1;
    -3  1  1  3 -1  1  3  1 -3  1 -3  1  1 -1 -1  3 -1 -3  3 -3 -3 -3  1  1;
    1  1 -1 -1  3 -3 -3  3 -3  1 -1 -1  1 -1  1  1 -1 -3 -1  1 -1  3 -1 -3;
    -3  3  3 -1 -1 -3 -1  3  1  3  1  3  1  1 -1  3  1 -1  1  3 -3 -1 -1  1;
    -3  1  3 -3  1 -1 -3  3 -3  3 -1 -1 -1 -1  1 -3 -3 -3  1 -3 -3 -3  1 -3;
    1  1 -3  3  3 -1 -3 -1  3 -3  3  3  3 -1  1  1 -3  1 -1  1  1 -3  1  1;
    -1  1 -3 -3  3 -1  3 -1 -1 -3 -3 -3 -1 -3 -3  1 -1  1  3  3 -1  1 -1  3;
    1  3  3 -3 -3  1  3  1 -1 -3 -3 -3  3  3 -3  3  3 -1 -3  3 -1  1 -3  1;
    1  3  3  1  1  1 -1 -1  1 -3  3 -1  1  1 -3  3  3 -1 -3  3 -3 -1 -3 -1;
    3 -1 -1 -1 -1 -3 -1  3  3  1 -1  1  3  3  3 -1  1  1 -3  1  3 -1 -3  3;
    -3 -3  3  1  3  1 -3  3  1  3  1  1  3  3 -1 -1 -3  1 -3 -1  3  1  1  3;
    -1 -1  1 -3  1  3 -3  1 -1 -3 -1  3  1  3  1 -1 -3 -3 -1 -1 -3 -3 -3 -1;
    -1 -3  3 -1 -1 -1 -1  1  1 -3  3  1  3  3  1 -1  1 -3  1 -3  1  1 -3 -1;
    1  3 -1  3  3 -1 -3  1 -1 -3  3  3  3 -1  1  1  3 -1 -3 -1  3 -1 -1 -1;
    1  1  1  1  1 -1  3 -1 -3  1  1  3 -3  1 -3 -1  1  1 -3 -3  3  1  1 -3;
    1  3  3  1 -1 -3  3 -1  3  3  3 -3  1 -1  1 -1 -3 -1  1  3 -1  3 -3 -3;
    -1 -3  3 -3 -3 -3 -1 -1 -3 -1 -3  3  1  3 -3 -1  3 -1  1 -1  3 -3  1 -1;
    -3 -3  1  1 -1  1 -1  1 -1  3  1 -3 -1  1 -1  1 -1 -1  3  3 -3 -1  1 -3;
    -3 -1 -3  3  1 -1 -3 -1 -3 -3  3 -3  3 -3 -1  1  3  1 -3  1  3  3 -1 -3;
    -1 -1 -1 -1  3  3  3  1  3  3 -3  1  3 -1  3 -1  3  3 -3  3  1 -1  3  3;
    1 -1  3  3 -1 -3  3 -3 -1 -1  3 -1  3 -1 -1  1  1  1  1 -1 -1 -3 -1  3;
    1 -1  1 -1  3 -1  3  1  1 -1 -1 -3  1  1 -3  1  3 -3  1  1 -3 -3 -1 -1;
    -3 -1  1  3  1  1 -3 -1 -1 -3  3 -3  3  1 -3  3 -3  1 -1  1 -3  1  1  1;
    -1 -3  3  3  1  1  3 -1 -3 -1 -1 -1  3  1 -3 -3 -1  3 -3 -1 -3 -1 -3 -1;
    -1 -3 -1 -1  1 -3 -1 -1  1 -1 -3  1  1 -3  1 -3 -3  3  1  1 -1  3 -1 -1;
    1  1 -1 -1 -3 -1  3 -1  3 -1  1  3  1 -1  3  1  3 -3 -3  1 -1 -1  1  3 ];

if strcmp(LTE_params.CyclicPrefix,'extended')
    NsymbUL = 6;               % number of SC-FDMA symbols in uplink slot, (page 13 tab. 5.2.3-1)
else
    NsymbUL = 7;               % number of SC-FDMA symbols in uplink slot, (page 13 tab. 5.2.3-1)
end

Nc = 1600;
Npn = 8*NsymbUL*(subframe+1)*LTE_params.Nslot+8;
NumBit = Npn+Nc;
c_init = zeros(1,31);
c_init(31) = 1;

pn_gen = commsrc.pn('GenPoly', [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1], ...
    'InitialStates',    c_init,   ...
    'Shift',         0,               ...
    'NumBitsOut',    NumBit);
pn_seq1 = generate(pn_gen);

%% selection of sequence group (section 5.5.1.3)
nIDRS = BS.NIDcell;       % values for SRS
fss = mod(nIDRS,30);
switch transmission
    case 'PUCCH'
        if not(isnan(nID_PUCCH))
            nIDRS = nID_PUCCH;
        end
        fss = mod(nIDRS,30);
    case 'PUSCH'
        if not(isnan(nID_PUSCH))
            nIDRS = nID_PUSCH;
            fss = mod(nIDRS,30);
        else
            fss = mod(BS.NIDcell + deltass,30);
        end
end

fgh = zeros(size(ns));
if groupHopping
    ii = 0:7;
    
   
    
    c_init = de2bi(floor(nIDRS/30),31,'left-msb');
    pn_gen = commsrc.pn('GenPoly', [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1], ...
        'InitialStates',   c_init,   ...
        'Shift',         0,               ...
        'NumBitsOut',    NumBit);
    pn_seq2 = generate(pn_gen);
    
    for jj=ns
        pn_seq_c = xor(pn_seq1(8*jj+ii+Nc+1),pn_seq2(8*jj+ii+Nc+1));
        fgh(mod(jj-1,LTE_params.Nslot)+1) = bi2de(pn_seq_c');
    end  
end
u = mod(fgh + fss,30);

%% selection of sequence (section 5.5.1.4)
v = zeros(size(ns));

if strcmp(type,'srs')
    [mSRS Nb mSRS0] = LTE_UL_mSRS(Csrs,Bsrs,Nrb); 
    Msc = mSRS(Bsrs + 1)*Nsc/2;
    m = Msc / Nsc;
else
    Msc = m * Nsc;
end

if m >= 6 && not(groupHopping) && sequenceHopping
    if STRCMP(type,'srs')
        c_init = de2bi(floor(nIDRS/30)*2^5 + mod(nIDRS + deltass,30),31,'left-msb');
    else
        c_init = de2bi(floor(nIDRS/30)*2^5 + fss,31,'left-msb');        
    end
    pn_gen = commsrc.pn('GenPoly', [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1], ...
        'InitialStates',   c_init,   ...
        'Shift',         0,               ...
        'NumBitsOut',    NumBit);
    pn_seq2 = generate(pn_gen);

    for jj=ns
        v(mod(jj-1,LTE_params.Nslot)+1) = xor(pn_seq1(jj+Nc+1),pn_seq2(jj+Nc+1));
    end  
end

r = zeros(length(u),Msc);
switch m
    case 1
        for ii=1:size(r,1)
            r(ii,:) = exp( 1i * fi1(u(ii)+1,:) * pi/4 );
        end
    case 2
        for ii=1:size(r,1)
            r(ii,:) = exp( 1i * fi2(u(ii)+1,:) * pi/4 );
        end
    otherwise
        NzcRS = max(primes(Msc));
        m = 0 : NzcRS - 1;
        q_ = NzcRS*(u+1)/31;
        q = floor(q_+1/2)+v.*(-1).^(floor(2*q_));

        for ii=1:size(r,1)
            x_q = exp( -1i * pi * q(ii) * m .*( m + 1 )/NzcRS );
            temp = repmat(x_q,1,ceil(Msc/NzcRS));
            r(ii,:) = temp(1:Msc);
        end
end

%% Generating DM RS
% select base sequence 
% 1 base sequence for m < 1
% 2 seguences v = 0,1 for m > 5
% compute alfa
switch type
    case 'dmrs'
        %NumBit = NmaxRB*4;
        nDRMS1 = [0 2 3 4 6 8 9 10];    % tab. 5.5.2.1.1-2

        [w, nDMRS2] = LTE_UL_cyclicShift(CSFiDCIFormat);
        
        % spratsch: MUMIMO modification
%         if (LTE_params.UE_config.mode == 5)
%             nDMRS2 = mod(nDMRS2 + 12/LTE_params.nUE * (uu-1), 12);
%         end
        
        if(LTE_params.activate_DMRS_with_OCC == false)  % 5.5.2.1.1 if OCC is activated or not 
            w = ones(size(w));
        end
        
        % Pseudo-random sequence sec:7.2                 
        c_init = de2bi(floor(BS.NIDcell/30)*2^5 + mod(BS.NIDcell + deltass,30),31,'left-msb');
        pn_gen = commsrc.pn('GenPoly', [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1], ...
            'InitialStates',   c_init,   ...
            'Shift',         0,               ...
            'NumBitsOut',    NumBit);
        pn_seq2 = generate(pn_gen);
        
        ii = 0:7;                       % 36.211 sec. 5.5.2.1.1
        nprs = zeros(size(ns));
        
        for jj=ns
            pn_seq_c = xor(pn_seq1(8*NsymbUL*jj+ii+Nc+1),pn_seq2(8*NsymbUL*jj+ii+Nc+1));
            nprs(mod(jj-1,LTE_params.Nslot)+1) = bi2de(pn_seq_c');
        end

        n_cs = mod((nDRMS1(cyclicShift+1)+ nDMRS2(layer)+nprs),12);
        alpha = pi*n_cs / 6;

        n=0:Msc-1;
        rs = transpose([w(layer,1)*exp(1i*alpha(1)*n).*r(1,:); w(layer,2)*exp(1i*alpha(2)*n).*r(2,:)]);
        % addiditonaly store the base sequence
        rs_base = transpose([w(layer,1).*r(1,:); w(layer,2).*r(2,:)]);
        
        % write to genie
        if(layer==1)
            % if thie is the first layer to generate the DMRS
            % clear the DMRS in the genie before you overwrite it
            UE_output.UE_genie.rs = [];             % the DMRS
            UE_output.UE_genie.rs_base_seq = [];    % and the base sequence
        end
        UE_output.UE_genie.rs(:,:,layer) = rs;                  % user allocation of RE
        UE_output.UE_genie.rs_base_seq(:,:,layer) = rs_base;    % store base sequence for CE
        UE_output.UE_genie.CSFiDCIFormat = CSFiDCIFormat;
        
    case 'srs'                                  % see TS 36.211 section 5.5.3
        if frameStructure == 2 && (subframe == 1 || (subframe == 6 && Nsp == 2))
            UpPTS = true;
        else
            UpPTS = false;
        end
        
        Nap = LTE_params.UE_config.nTX;
        p = 0:Nap - 1;
        nCS_SRS = mod(nCS_SRS + 8*p/Nap,8);
        alpha = pi*nCS_SRS/4;
        
        if Nap == 4
            kTC = kTC + (1 - 2*kTC)*mod(p,2).*(nCS_SRS >= 4);
        else
            kTC = kTC*ones(Nap,1);
        end
                        
        if b_hop >= Bsrs && Tsrs == 0
            nb = mod(floor(4*n_RRC./mSRS),Nb);
        else
            if Tsrs == 2 && frameStructure == 2
                nSRS = 2*Nsp*nf + 2*(Nsp - 1)*floor(ns/10) + floor(Toffset/Toffset_max);
            else
                nSRS = floor(nf*10 + subframe)/Tsrs;
            end
            
            nb = mod(floor(4*nRRC./mSRS),Nb);
            prod1 = 1;
            prod2 = 1;
            
            Nb(b_hop + 1) = 1;
            for b = b_hop + 2:Bsrs + 1
                prod1 = prod1*Nb(b);
                
                if mod(Nb(b),2) == 0
                    Fb = Nb(b)/2*(floor(mod(nSRS,prod1)/prod2) + floor(mod(nSRS,prod1)/(2*prod2)));
                else
                    Fb = floor(Nb(b)/2)*floor(nSRS/prod2);
                end
                
                prod2 = prod1;
                nb(b) = mod(Fb + floor(4*nRRC/mSRS(b)),Nb(b));
            end
        end
        
        if UpPTS            
            n_hf = (subframe > 5);
            
            if mod(mod(nf,2)*(2-Nsp) + n_hf,2) == 0 
                if srsMaxUpPTS
                    Nra = 0;   
                    mSRS_max = max(mSRS0.*(mSRS0 <= Nrb - 6*Nra));
                else
                    mSRS_max = mSRS(1);
                end
                k0 = (Nrb - mSRS_max)*Nsc + kTC;
            else
                k0 = kTC;
            end
        else
            k0 = (floor(Nrb/2) - mSRS(1)/2)*Nsc + kTC;
        end
        k0 = k0 + Nsc*sum(mSRS.*nb);
                
        n = (0:Msc-1)';
        r = transpose(r(1,:));
        rs = zeros(LTE_params.Ntot,Nap);
        
        for ii = 1:Nap
            rs(k0(ii)+(1:2:2*Msc),ii) = exp(1i*n*alpha(ii)).*r/sqrt(Nap); 
        end
end        


 