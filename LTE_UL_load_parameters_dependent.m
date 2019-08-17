% Calculate configuration parameters which are dependent on the main
% configuration file.
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

LTE_params.N_subframes = N_subframes;


%% SNR vector
LTE_params.SNR_vec = SNR_vec;


%% connection_table generate a standardised one if non specified
if ~exist('connection_table','var')  % automatically generate connection_table form nBS and nUE
    LTE_params.connection_table = utils.generate_connection_table(LTE_params.nBS, LTE_params.nUE);    
else
    if LTE_params.connection_table == 0
       LTE_params.connection_table = utils.generate_connection_table(LTE_params.nBS, LTE_params.nUE);
    else
        LTE_params.connection_table = LTE_params.connection_table; %emphasise that the connection table is not altered
    end
end

if ~exist('delay_table','var')
    LTE_params.delay_table = zeros(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
else
    LTE_params.delay_table = zeros(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);  % asynchronous network case, change here.
end

%% Calculate number of available Resource Blocks
% DEFINED IN STANDARD 3GPP TS 36.211 V8.2.0 (2008-03)
LTE_params.ResourceBlock = 180e3; % fixed badwidth of resource block in Hz, page 33
LTE_params.Nsc = LTE_params.ResourceBlock/LTE_params.SubcarrierSpacing; % number of subcarriers in one resource block, fixed length of resource block in Hz, page 33
if(LTE_params.Bandwidth == 1.4e6)
    LTE_params.Nrb = 6; % number of resource blocks
else
    LTE_params.Nrb = (LTE_params.Bandwidth*0.9) / LTE_params.ResourceBlock; % number of resource blocks, transmission BW is 90% of the total BW
end
LTE_params.Ntot = LTE_params.Nsc*LTE_params.Nrb; % Total number of subcarriers not NULL

% %% Calculate CQI feedback granularity
% if isnumeric(LTE_params.UE_config.CQI_fb_granularity)                   % insert number of clustersize for CQI feedback
%     if((LTE_params.UE_config.CQI_fb_granularity <=0) || (LTE_params.UE_config.CQI_fb_granularity >LTE_params.Nrb) || (mod(LTE_params.Nrb,LTE_params.UE_config.CQI_fb_granularity)))
%         %  granularity has to be positive and must divide total number of RBs
%         error('CQI feedback granularity (cluster size) not valid');
%     end
% else
%     switch LTE_params.UE_config.CQI_fb_granularity                      % choose coarse or fine
%         case 'fine'
%             LTE_params.UE_config.CQI_fb_granularity = 1;                % coarse feeback, only one cluster for feedback
%         case 'coarse'
%             LTE_params.UE_config.CQI_fb_granularity = LTE_params.Nrb;   % fine feeback, each RB is a feeback cluster
%         otherwise
%             error('CQI feedback granularity (cluster size) not valid');
%     end
% end

%% Calculate CP duration, number of OFDM symbols
% DEFINED IN STANDARD 3GPP TS 36.211 V8.2.0 (2008-03)
LTE_params.FrameDur = 10e-3; % fixed frame duration 10 ms, page 9

%Table 6.12-1: OFDM parameters + page 9
if (strcmp(LTE_params.CyclicPrefix,'normal') && LTE_params.SubcarrierSpacing == 15e3)
    LTE_params.Tg = zeros(1,2);
    LTE_params.Tg(1) = 160/(15000*2048); % normal CP time - 1. OFDM symbol in slot
    LTE_params.Tg(2) = 144/(15000*2048); % normal CP time - remaining OFDM symbols in slot
    %Table 6.2.3-1: Physical resource blocks parameters.
    LTE_params.Ns = 7; % number of OFDM symbols in one slot
    LTE_params.Nsub = 2*7; % number of OFDM symbols in one subframe = 2*in one slot
elseif(strcmp(LTE_params.CyclicPrefix,'extended') && LTE_params.SubcarrierSpacing == 15e3)
    LTE_params.Tg = 512/(15000*2048); % extended CP time for all OFDM symbols in slot
    %Table 6.2.3-1: Physical resource blocks parameters.
    LTE_params.Ns = 6; % number of OFDM symbols in one slot
    LTE_params.Nsub = 2*6; % number of OFDM symbols in one subframe = 2*in one slot
else
    error('Wrong combination of value of subcarrier spacing and type of cyclic prefix')
end

%% SRS parameters (see TS 36.211 section 5.5.3.3)
if LTE_params.frameStructure == 1       % FDD
    switch LTE_params.srsSubframeConfig
        case 0
            LTE_params.SRS_period = 1;
            LTE_params.SRS_offset = 0;
        case 1
            LTE_params.SRS_period = 2;
            LTE_params.SRS_offset = 0;
        case 2
            LTE_params.SRS_period = 2;
            LTE_params.SRS_offset = 1;
        case 3
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = 0;
        case 4
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = 1;
        case 5
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = 2;
        case 6
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = 3;
        case 7
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = [0 1];
        case 8
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = [2 3];
        case 9
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = 0;
        case 10
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = 1;
        case 11
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = 2;
        case 12
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = 3;
        case 13
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = [0 1 2 3 4 6 8];
        case 14
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = [0 1 2 3 4 5 6 8];
        case 15
            LTE_params.SRS_period = 0;
            LTE_params.SRS_offset = 0;
        otherwise
            error('SRS SubframeConfig not supported');
    end
elseif LTE_params.frameStructure == 2      % TDD
    switch LTE_params.srsSubframeConfig
        case 0
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = 1;
        case 1
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = [1 2];
        case 2
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = [1 3];
        case 3
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = [1 4];
        case 4
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = [1 2 3];
        case 5
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = [1 2 4];
        case 6
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = [1 3 4];
        case 7
            LTE_params.SRS_period = 5;
            LTE_params.SRS_offset = [1 2 3 4];
        case 8
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = [1 2 6];
        case 9
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = [1 3 6];
        case 10
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = [1 6 7];
        case 11
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = [1 2 6 8];
        case 12
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = [1 3 6 9];
        case 13
            LTE_params.SRS_period = 10;
            LTE_params.SRS_offset = [1 4 6 7];
        case 15
            LTE_params.SRS_period = 0;
            LTE_params.SRS_offset = 0;
        otherwise
            error('SRS SubframeConfig not supported');
    end
else
    error('framestructure not supported');
end

if strcmp(LTE_params.CyclicPrefix,'normal')
    LTE_params.UpPTS_length = LTE_params.specialSubframeConf;
else
    LTE_params.UpPTS_length = LTE_params.specialSubframeConf + 1;
end
LTE_params.UpPTS_length = floor(LTE_params.UpPTS_length / 5) + 1;

%% Create ChanMod object (channel model)
ChanMod.filtering            = LTE_params.ChanMod_config.filtering;
ChanMod.interpolation_method = LTE_params.ChanMod_config.interpolation_method;
ChanMod.type                 = LTE_params.ChanMod_config.type;
ChanMod.nTX                  = LTE_params.UE_config.nTX;
ChanMod.nRX                  = LTE_params.BS_config.nRX;
ChanMod.corr_coefRX          = LTE_params.ChanMod_config.corr_coefRX;
ChanMod.corr_coefTX          = LTE_params.ChanMod_config.corr_coefTX;
ChanMod.sin_num              = LTE_params.ChanMod_config.sin_num;
ChanMod.time_correlation     = LTE_params.ChanMod_config.time_correlation;

switch ChanMod.type
    case {'PedA'}
        ChanMod.PDP_dB = [0 -9.7 -19.2 -22.8;  % Average power [dB]
            0 110*10^-9 190*10^-9 410*10^-9]; % delay (s)
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));
    case {'PedB', 'PedBcorr'}
        ChanMod.PDP_dB = [0   -0.9  -4.9  -8    -7.8  -23.9; % Average power [dB]
            0 200*10^-9 800*10^-9 1200*10^-9 2300*10^-9 3700*10^-9]; % delay (s)
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));
    case {'VehA'}
        ChanMod.PDP_dB = [0   -1  -9  -10    -15  -20; % Average power [dB]
            0 310*10^-9 710*10^-9 1090*10^-9 1730*10^-9 2510*10^-9]; % delay (s)
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));
    case {'VehB'}       
        ChanMod.PDP_dB = [-2.5   0  -12.8  -10    -25.2  -16; % Average power [dB]
            0 300*10^-9 8900*10^-9 12900*10^-9 17100*10^-9 20000*10^-9]; % delay (s)
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));
    case {'TU'}
        ChanMod.PDP_dB = [-5.7000 -7.6000 -10.1000 -10.2000 -10.2000 -11.5000 -13.4000 -16.3000 -16.9000 -17.1000 -17.4000,...
            -19.0000 -19.0000 -19.8000 -21.5000 -21.6000 -22.1000 -22.6000 -23.5000 -24.3000; % Average power [dB]
            0 0.2170 0.5120 0.5140 0.5170 0.6740 0.8820 1.2300 1.2870 1.3110 1.3490 1.5330 1.5350,...
            1.6220 1.8180 1.8360 1.8840 1.9430 2.0480 2.1400];% delay (us)
        ChanMod.PDP_dB(2,:) = ChanMod.PDP_dB(2,:)*10^-6;
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));
    case {'RA'}
        ChanMod.PDP_dB = [-5.2000 -6.4000 -8.4000 -9.3000 -10.0000 -13.1000 -15.3000 -18.5000 -20.4000 -22.4000; % Average power [dB]
            0 0.0420 0.1010 0.1290 0.1490 0.2450 0.3120 0.4100 0.4690 0.5280]; % delay (us)
        ChanMod.PDP_dB(2,:) = ChanMod.PDP_dB(2,:)*10^-6;
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));
    case {'HT'}
        ChanMod.PDP_dB = [-3.6000 -8.9000 -10.2000 -11.5000 -11.8000 -12.7000 -13.0000 -16.2000 -17.3000 -17.700 -17.6000 -22.7000,...
            -24.1000 -25.8000 -25.8000 -26.2000 -29.0000 -29.9000 -30.0000 -30.7000; % Average power [dB]
            0 0.3560 0.4410 0.5280 0.5460 0.6090 0.6250 0.8420 0.9160 0.9410 15.0000 16.1720 16.4920 16.8760 16.8820,...
            16.9780 17.6150 17.827 17.8490 18.0160]; % delay (us)
        ChanMod.PDP_dB(2,:) = ChanMod.PDP_dB(2,:)*10^-6;
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));
    case {'EPedA'}  % extended multipath channel model in 36.101 added, 02092009
        ChanMod.PDP_dB = [0 -1.0,    -2.0     -3.0     -8.0      -17.2    -20.8;  % Average power [dB]
            0 30*10^-9 70*10^-9 90*10^-9 110*10^-9 190*10^-9 410*10^-9]; % delay (s)
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));
    case {'EVehA'}  % extended multipath channel model in 36.101 added, 02092009
        ChanMod.PDP_dB = [0   -1.5   -1.4  -3.6 -0.6 -9.1 -7.0 -12.0 -16.9; % Average power [dB]
            0 30*10^-9 150*10^-9 310*10^-9 370*10^-9 710*10^-9 1090*10^-9 1730*10^-9 2510*10^-9]; % delay (s)
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));
    case {'ETU'}    % extended multipath channel model in 36.101 added, 02092009
        ChanMod.PDP_dB = [-1.0 -1.0 -1.0 0.0 0.0 0.0 -3.0 -5.0 -7.0; % Average power [dB]
            0 50*10^-9 120*10^-9 200*10^-9 230*10^-9 500*10^-9 1600*10^-9 2300*10^-9 5000*10^-9];% delay (s)
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));
    case {'ePDP'} % exponential PDP with pre-specified rms-delay spread
        sample_interv = 1/LTE_params.SubcarrierSpacing*1/2^ceil(log2(LTE_params.Ntot))*5;
        tau_rms_samples = LTE_params.ChanMod_config.tau_rms/sample_interv;
        syms x
        tau = solve('exp(-1/x)*(1+exp(-1/x))/(1-exp(-1/x))^2 - (exp(-1/x)/(1-exp(-1/x)))^2 = tau_rms_samples^2','IgnoreAnalyticConstraints',true);
        tau = eval(tau(eval(tau) > 0));
        k = 0:10*tau_rms_samples;
        Ps = (1-exp(-1/tau))*exp(-k/tau);
        Ps = Ps(Ps > Ps*0.001);
        ChanMod.PDP_dB = [10*log10(Ps);(0:length(Ps)-1)*sample_interv];
        ChanMod.normH = sqrt(sum(10.^(ChanMod.PDP_dB(1,:)/10)));      
    case  {'AWGN','flat Rayleigh','winner_II','Rayleigh','TR 36.873'}
        % Do nothing (these are not channels defined by a PDP)
    otherwise
        error('Channel not supported');
        
end

%% Store interleaver sequences and the sort for the turbo coding and rate matching
% The sub-block interleaver permutation pattern (column permutation)
LTE_params.sub_block_interleaver_permutation_pattern = [0,16,8,24,4,20,12,28,2,18,10,26,6,22,14,30,1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31];
LTE_params.sub_block_interleaver_permutation_pattern_plus_one = LTE_params.sub_block_interleaver_permutation_pattern+1;
% Turbo encoder/decoder interleaving table
LTE_params.turbo_interleaver_table = [40,3,10;48,7,12;56,19,42;64,7,16;72,7,18;80,11,20;88,5,22;96,11,24;104,7,26;112,41,84;120,103,90;128,15,32;136,9,34;144,17,108;152,9,38;160,21,120;168,101,84;176,21,44;184,57,46;192,23,48;200,13,50;208,27,52;216,11,36;224,27,56;232,85,58;240,29,60;248,33,62;256,15,32;264,17,198;272,33,68;280,103,210;288,19,36;296,19,74;304,37,76;312,19,78;320,21,120;328,21,82;336,115,84;344,193,86;352,21,44;360,133,90;368,81,46;376,45,94;384,23,48;392,243,98;400,151,40;408,155,102;416,25,52;424,51,106;432,47,72;440,91,110;448,29,168;456,29,114;464,247,58;472,29,118;480,89,180;488,91,122;496,157,62;504,55,84;512,31,64;528,17,66;544,35,68;560,227,420;576,65,96;592,19,74;608,37,76;624,41,234;640,39,80;656,185,82;672,43,252;688,21,86;704,155,44;720,79,120;736,139,92;752,23,94;768,217,48;784,25,98;800,17,80;816,127,102;832,25,52;848,239,106;864,17,48;880,137,110;896,215,112;912,29,114;928,15,58;944,147,118;960,29,60;976,59,122;992,65,124;1008,55,84;1024,31,64;1056,17,66;1088,171,204;1120,67,140;1152,35,72;1184,19,74;1216,39,76;1248,19,78;1280,199,240;1312,21,82;1344,211,252;1376,21,86;1408,43,88;1440,149,60;1472,45,92;1504,49,846;1536,71,48;1568,13,28;1600,17,80;1632,25,102;1664,183,104;1696,55,954;1728,127,96;1760,27,110;1792,29,112;1824,29,114;1856,57,116;1888,45,354;1920,31,120;1952,59,610;1984,185,124;2016,113,420;2048,31,64;2112,17,66;2176,171,136;2240,209,420;2304,253,216;2368,367,444;2432,265,456;2496,181,468;2560,39,80;2624,27,164;2688,127,504;2752,143,172;2816,43,88;2880,29,300;2944,45,92;3008,157,188;3072,47,96;3136,13,28;3200,111,240;3264,443,204;3328,51,104;3392,51,212;3456,451,192;3520,257,220;3584,57,336;3648,313,228;3712,271,232;3776,179,236;3840,331,120;3904,363,244;3968,375,248;4032,127,168;4096,31,64;4160,33,130;4224,43,264;4288,33,134;4352,477,408;4416,35,138;4480,233,280;4544,357,142;4608,337,480;4672,37,146;4736,71,444;4800,71,120;4864,37,152;4928,39,462;4992,127,234;5056,39,158;5120,39,80;5184,31,96;5248,113,902;5312,41,166;5376,251,336;5440,43,170;5504,21,86;5568,43,174;5632,45,176;5696,45,178;5760,161,120;5824,89,182;5888,323,184;5952,47,186;6016,23,94;6080,47,190;6144,263,480];

%% Define symbol alphabet
% DEFINED IN STANDARD 3GPP TS 36.211 V8.2.0 (2008-03), page 61-63
LTE_params.SymbolAlphabet{1} = [ 1+1j, -1-1j].'/sqrt(2);
LTE_params.SymbolAlphabet{2} = [ 1+1j, 1-1j, -1+1j, -1-1j].'/sqrt(2);
LTE_params.SymbolAlphabet{4} = [
    complex( 1,  1)
    complex( 1,  3)
    complex( 3,  1)
    complex( 3,  3)
    complex( 1, -1)
    complex( 1, -3)
    complex( 3, -1)
    complex( 3, -3)
    complex(-1,  1)
    complex(-1,  3)
    complex(-3,  1)
    complex(-3,  3)
    complex(-1, -1)
    complex(-1, -3)
    complex(-3, -1)
    complex(-3, -3)
    ] / sqrt(10);
LTE_params.SymbolAlphabet{6} = [
    complex( 3,  3)
    complex( 3,  1)
    complex( 1,  3)
    complex( 1,  1)
    complex( 3,  5)
    complex( 3,  7)
    complex( 1,  5)
    complex( 1,  7)
    complex( 5,  3)
    complex( 5,  1)
    complex( 7,  3)
    complex( 7,  1)
    complex( 5,  5)
    complex( 5,  7)
    complex( 7,  5)
    complex( 7,  7) % symbol 0-15
    complex( 3, -3)
    complex( 3, -1)
    complex( 1, -3)
    complex( 1, -1)
    complex( 3, -5)
    complex( 3, -7)
    complex( 1, -5)
    complex( 1, -7)
    complex( 5, -3)
    complex( 5, -1)
    complex( 7, -3)
    complex( 7, -1)
    complex( 5, -5)
    complex( 5, -7)
    complex( 7, -5)
    complex( 7, -7) % symbol 16-31
    complex(-3,  3)
    complex(-3,  1)
    complex(-1,  3)
    complex(-1,  1)
    complex(-3,  5)
    complex(-3,  7)
    complex(-1,  5)
    complex(-1,  7)
    complex(-5,  3)
    complex(-5,  1)
    complex(-7,  3)
    complex(-7,  1)
    complex(-5,  5)
    complex(-5,  7)
    complex(-7,  5)
    complex(-7,  7) % symbol 32-47
    complex(-3, -3)
    complex(-3, -1)
    complex(-1, -3)
    complex(-1, -1)
    complex(-3, -5)
    complex(-3, -7)
    complex(-1, -5)
    complex(-1, -7)
    complex(-5, -3)
    complex(-5, -1)
    complex(-7, -3)
    complex(-7, -1)
    complex(-5, -5)
    complex(-5, -7)
    complex(-7, -5)
    complex(-7, -7) ] / sqrt(42); % symbol 48-63

LTE_params.bittable{1} = logical([0,1]);
LTE_params.bittable{2} = logical(   [0,1,0,1;
    0,0,1,1]);
LTE_params.bittable{4} = logical(   [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1;
    0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
    0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
    0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1]);
LTE_params.bittable{6} = logical(   [   0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1;
    0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
    0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
    0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]);

%% Calculate FFT lengths, CP length and indexes, OFDM symbol, slot and subframe duration, # transmit symbols per subframe, sampling time
if(LTE_params.Bandwidth == 15e6 && LTE_params.SubcarrierSpacing == 15e3)
    LTE_params.Nfft = 1536;                             % number of FFT points
else
    LTE_params.Nfft =  2^ceil(log2(LTE_params.Ntot));   % number of FFT points
end

LTE_params.Fs = LTE_params.SubcarrierSpacing*LTE_params.Nfft; % sampling frequency
LTE_params.Tb = 1/LTE_params.SubcarrierSpacing; % useful Symbol Time
if(length(LTE_params.Tg)==2)    % normal CP
    LTE_params.Ng = zeros(1,2);
    LTE_params.Ng(1) = LTE_params.Tg(1)*LTE_params.Fs; % number of CP points of normal CP for 1. first OFDM symbol in slot
    LTE_params.Ng(2) = round(LTE_params.Tg(2)*LTE_params.Fs); % number of CP points of normal CP for remaining OFDM symbols in slot
    LTE_params.Index_TxCyclicPrefix{1} = [LTE_params.Nfft-LTE_params.Ng(1)+1:LTE_params.Nfft 1:LTE_params.Nfft];
    LTE_params.Index_TxCyclicPrefix{2} = [LTE_params.Nfft-LTE_params.Ng(2)+1:LTE_params.Nfft 1:LTE_params.Nfft];
    
    LTE_params.Ts = zeros(1,2);
    LTE_params.Ts(1) = LTE_params.Tb + LTE_params.Tg(1); % 1. OFDM symbol time
    LTE_params.Ts(2) = LTE_params.Tb + LTE_params.Tg(2); % remaining OFDM symbol time
    
    LTE_params.Tslot = LTE_params.Tg(1) + (LTE_params.Ns-1)*LTE_params.Tg(2) + LTE_params.Ns*LTE_params.Tb; % fixed duration of the slot = 0.5ms, page 9
    LTE_params.Tsubframe = 2*LTE_params.Tslot; % fixed duration of the subframe = 1ms, page 9
    LTE_params.TxSymbols = 2*(LTE_params.Ng(1) + (LTE_params.Ns-1)*LTE_params.Ng(2) + LTE_params.Nfft*LTE_params.Ns);  % number of symbols in one subframe
    
    %DERIVED PARAMETERS FOR RECEIVER
    LTE_params.NfftCP{1} = length(LTE_params.Index_TxCyclicPrefix{1});
    LTE_params.NfftCP{2} = length(LTE_params.Index_TxCyclicPrefix{2});
    LTE_params.Index_RxCyclicPrefix{1} = LTE_params.Ng(1)+1: LTE_params.NfftCP{1};
    LTE_params.Index_RxCyclicPrefix{2} = LTE_params.Ng(2)+1: LTE_params.NfftCP{2};
else                            % extended CP
    LTE_params.Ng = LTE_params.Tg*LTE_params.Fs;    % number of CP points of extended CP for OFDM symbols in slot
    LTE_params.Index_TxCyclicPrefix{1} = [LTE_params.Nfft-LTE_params.Ng+1:LTE_params.Nfft 1:LTE_params.Nfft];
    LTE_params.Index_TxCyclicPrefix{2} = LTE_params.Index_TxCyclicPrefix{1};    % both the same for extended CP
    
    LTE_params.Ts = LTE_params.Tb + LTE_params.Tg;  % OFDM symbol time
    
    LTE_params.Tslot = LTE_params.Ns*LTE_params.Tg + LTE_params.Ns*LTE_params.Tb; % fixed duration of the slot = 0.5ms, page 9
    LTE_params.Tsubframe = 2*LTE_params.Tslot;      % fixed duration of the slot = 1ms, page 9
    LTE_params.TxSymbols = 2*(LTE_params.Ns*LTE_params.Ng + LTE_params.Nfft*LTE_params.Ns);  % number of symbols in one subframe
    
    %DERIVED PARAMETERS FOR RECEIVER
    LTE_params.NfftCP{1} = length(LTE_params.Index_TxCyclicPrefix{1});
    LTE_params.NfftCP{2} = LTE_params.NfftCP{1};    % both the same for extended CP
    LTE_params.Index_RxCyclicPrefix{1} = LTE_params.Ng+1: LTE_params.NfftCP{1};
    LTE_params.Index_RxCyclicPrefix{2} = LTE_params.Index_RxCyclicPrefix{1}; % both the same for extended CP
end

LTE_params.SamplingTime = LTE_params.Tb/LTE_params.Nfft;
% Number of OFDM Symbols transmitted in one frame
LTE_params.OFDMNsym = round(LTE_params.FrameDur/LTE_params.Ts(1)/4)*4; % number of OFDM symbols in one frame

%% ICI emulation parameters
LTE_params.ICI_power = zeros(1,LTE_params.Ntot);
if (~LTE_params.UE_config.ignore_ISI_ICI) || (LTE_params.UE_config.user_speed < eps)
    LTE_params.ICI_power = LTE_UL_estimate_ICI_power(LTE_params);
end

%% Store CQI parameters
% CQI 1 is index 1, CQI 2 is index 2, etc...
LTE_params.CQI_params(1).CQI = 1;
LTE_params.CQI_params(1).modulation = 'QPSK';
LTE_params.CQI_params(1).modulation_order = 2;
LTE_params.CQI_params(1).coding_rate_x_1024 = 78;
LTE_params.CQI_params(1).efficiency = 0.1523;

LTE_params.CQI_params(2).CQI = 2;
LTE_params.CQI_params(2).modulation = 'QPSK';
LTE_params.CQI_params(2).modulation_order = 2;
LTE_params.CQI_params(2).coding_rate_x_1024 = 120;
LTE_params.CQI_params(2).efficiency = 0.2344;

LTE_params.CQI_params(3).CQI = 3;
LTE_params.CQI_params(3).modulation = 'QPSK';
LTE_params.CQI_params(3).modulation_order = 2;
LTE_params.CQI_params(3).coding_rate_x_1024 = 193;
LTE_params.CQI_params(3).efficiency = 0.3770;

LTE_params.CQI_params(4).CQI = 4;
LTE_params.CQI_params(4).modulation = 'QPSK';
LTE_params.CQI_params(4).modulation_order = 2;
LTE_params.CQI_params(4).coding_rate_x_1024 = 308;
LTE_params.CQI_params(4).efficiency = 0.6016;

LTE_params.CQI_params(5).CQI = 5;
LTE_params.CQI_params(5).modulation = 'QPSK';
LTE_params.CQI_params(5).modulation_order = 2;
LTE_params.CQI_params(5).coding_rate_x_1024 = 449;
LTE_params.CQI_params(5).efficiency = 0.8770;

LTE_params.CQI_params(6).CQI = 6;
LTE_params.CQI_params(6).modulation = 'QPSK';
LTE_params.CQI_params(6).modulation_order = 2;
LTE_params.CQI_params(6).coding_rate_x_1024 = 602;
LTE_params.CQI_params(6).efficiency = 1.1758;

LTE_params.CQI_params(7).CQI = 7;
LTE_params.CQI_params(7).modulation = '16QAM';
LTE_params.CQI_params(7).modulation_order = 4;
LTE_params.CQI_params(7).coding_rate_x_1024 = 378;
LTE_params.CQI_params(7).efficiency = 1.4766;

LTE_params.CQI_params(8).CQI = 8;
LTE_params.CQI_params(8).modulation = '16QAM';
LTE_params.CQI_params(8).modulation_order = 4;
LTE_params.CQI_params(8).coding_rate_x_1024 = 490;
LTE_params.CQI_params(8).efficiency = 1.9141;

LTE_params.CQI_params(9).CQI = 9;
LTE_params.CQI_params(9).modulation = '16QAM';
LTE_params.CQI_params(9).modulation_order = 4;
LTE_params.CQI_params(9).coding_rate_x_1024 = 616;
LTE_params.CQI_params(9).efficiency = 2.4063;

LTE_params.CQI_params(10).CQI = 10;
LTE_params.CQI_params(10).modulation = '64QAM';
LTE_params.CQI_params(10).modulation_order = 6;
LTE_params.CQI_params(10).coding_rate_x_1024 = 466;
LTE_params.CQI_params(10).efficiency = 2.7305;

LTE_params.CQI_params(11).CQI = 11;
LTE_params.CQI_params(11).modulation = '64QAM';
LTE_params.CQI_params(11).modulation_order = 6;
LTE_params.CQI_params(11).coding_rate_x_1024 = 567;
LTE_params.CQI_params(11).efficiency = 3.3223;

LTE_params.CQI_params(12).CQI = 12;
LTE_params.CQI_params(12).modulation = '64QAM';
LTE_params.CQI_params(12).modulation_order = 6;
LTE_params.CQI_params(12).coding_rate_x_1024 = 666;
LTE_params.CQI_params(12).efficiency = 3.9023;

LTE_params.CQI_params(13).CQI = 13;
LTE_params.CQI_params(13).modulation = '64QAM';
LTE_params.CQI_params(13).modulation_order = 6;
LTE_params.CQI_params(13).coding_rate_x_1024 = 772;
LTE_params.CQI_params(13).efficiency = 4.5234;

LTE_params.CQI_params(14).CQI = 14;
LTE_params.CQI_params(14).modulation = '64QAM';
LTE_params.CQI_params(14).modulation_order = 6;
LTE_params.CQI_params(14).coding_rate_x_1024 = 873;
LTE_params.CQI_params(14).efficiency = 5.1152;

LTE_params.CQI_params(15).CQI = 15;
LTE_params.CQI_params(15).modulation = '64QAM';
LTE_params.CQI_params(15).modulation_order = 6;
LTE_params.CQI_params(15).coding_rate_x_1024 = 930; %938; %948
LTE_params.CQI_params(15).efficiency = 5.5547;

LTE_params.CQI_params(16).CQI = 16;
LTE_params.CQI_params(16).modulation = 'QPSK';
LTE_params.CQI_params(16).modulation_order = 2;
LTE_params.CQI_params(16).coding_rate_x_1024 = 1024/3;
LTE_params.CQI_params(16).efficiency = 2/3;

LTE_params.CQI_params(17).CQI = 17;
LTE_params.CQI_params(17).modulation = '16QAM';
LTE_params.CQI_params(17).modulation_order = 4;
LTE_params.CQI_params(17).coding_rate_x_1024 = 1024/2;
LTE_params.CQI_params(17).efficiency = 4/2;

% CQI 20 is used as a zero rate CQI (= CQI 0)

LTE_params.CQI_params(20).CQI = 20;
LTE_params.CQI_params(20).modulation = '0QAM';
LTE_params.CQI_params(20).modulation_order = 2;
LTE_params.CQI_params(20).coding_rate_x_1024 = 0;
LTE_params.CQI_params(20).efficiency = 0;

% CQIs used for HARQ testing
LTE_params.CQI_params(50).CQI = 50;
LTE_params.CQI_params(50).modulation = 'QPSK';
LTE_params.CQI_params(50).modulation_order = 2;
LTE_params.CQI_params(50).coding_rate_x_1024 = 1024/3;
LTE_params.CQI_params(50).efficiency = 2/3;

LTE_params.CQI_params(51).CQI = 51;
LTE_params.CQI_params(51).modulation = '16QAM';
LTE_params.CQI_params(51).modulation_order = 4;
LTE_params.CQI_params(51).coding_rate_x_1024 = 1024/3;
LTE_params.CQI_params(51).efficiency = 4/3;

LTE_params.CQI_params(52).CQI = 52;
LTE_params.CQI_params(52).modulation = '64QAM';
LTE_params.CQI_params(52).modulation_order = 6;
LTE_params.CQI_params(52).coding_rate_x_1024 = 1024/3;
LTE_params.CQI_params(52).efficiency = 6/3;

% Add RAN R1-07196 CQIs (for testing)
LTE_params.CQI_params(101:127) = LTE_common_R1_071667_CQI_parameters;

% Add HARQ test CQIs
LTE_params.CQI_params(55:55+25) = LTE_common_HARQ_CQI_parameters;

%% Create the correct resampler for the chanel model's PDP
switch ChanMod.type
    case {'AWGN','flat Rayleigh','winner_II','Rayleigh','TR 36.873'}
        % Do nothing
    otherwise
        switch ChanMod.interpolation_method
            case 'shift_to_nearest_neighbor'
                % create an appropriate object that does an equivalent thing
                ChanMod.interpolator = channel_resampling.nearestNeighborInterpolator(...
                    LTE_params.Fs,...                  % Sampling frequency (Hz)
                    size(ChanMod.PDP_dB,2),...         % Number of Taps
                    ChanMod.PDP_dB(2,:,:),...          % delay (s)
                    ChanMod.PDP_dB(1,:,:),...          % Average power [dB]
                    ChanMod.nTX,ChanMod.nRX);
                ChanMod.tap_delays = round(ChanMod.PDP_dB(2,:) * LTE_params.Fs);
            case 'sinc_interpolation'
                ChanMod.interpolator = channel_resampling.sincInterpolator(...
                    LTE_params.Fs,...                    % Sampling frequency (Hz)
                    ChanMod.PDP_dB(2,:,:),...            % delay (s)
                    ChanMod.PDP_dB(1,:,:),...            % Average power [dB]
                    false,...                            % Do not include the acausal part in the sinc interpolation
                    min(LTE_params.Tg)*LTE_params.Fs,... % Length of the impulse response. Set to the minimum length of the CP
                    ChanMod.nTX,ChanMod.nRX);
                
                ChanMod.tap_delays = min(LTE_params.Tg)*LTE_params.Fs;
        end
end

%% Set up random number generator for channel parameters
if LTE_params.random_channel_param_seeding
    LTE_params.channel_param_RandStream = RandStream('mt19937ar','Seed',LTE_params.channel_param_seed);
else
    LTE_params.channel_param_RandStream = RandStream('mt19937ar','Seed',rand*intmax);
end

%% Set up random number generators
if LTE_params.random_noise_seeding && strcmp(LTE_params.channel_matrix_source,'generated')
    LTE_params.noise_RandStream = RandStream('mt19937ar','Seed',LTE_params.noise_seed);
else
    LTE_params.noise_RandStream = RandStream('mt19937ar','Seed',rand*intmax);
end

%% Set up random number generator for transmit data
if LTE_params.random_data_seeding
    LTE_params.data_RandStream = RandStream('mt19937ar','Seed',LTE_params.data_seed);
else
    LTE_params.data_RandStream = RandStream('mt19937ar','Seed',rand*intmax);
end

%% Set up channel trace storage
if strcmp(LTE_params.channel_matrix_source,'generated') && LTE_params.store_channel_trace
    LTE_params.channel_matrix_trace = channels.channelMatrixTrace;
    LTE_params.channel_matrix_trace.noise_seed = LTE_params.noise_RandStream.Seed;
    
    LTE_params.channel_matrix_trace.bandwidth  = LTE_params.Bandwidth;
    LTE_params.channel_matrix_trace.type       = LTE_params.ChanMod_config.type;
    LTE_params.channel_matrix_trace.nTx        = LTE_params.BS_config.nRX;
    LTE_params.channel_matrix_trace.nRx        = LTE_params.UE_config.nTX;
    LTE_params.channel_matrix_trace.TTI_length = N_subframes;
    LTE_params.channel_matrix_trace.Nfft       = LTE_params.Nfft;
else
    LTE_params.store_channel_trace = false;
end

% Generate filename automatically
the_date = clock;
output_filename_UL = sprintf('%3.1fMHz_%s_%dx%d_%dTTI_%04d%02d%02d_%02d%02d%02d',...
    LTE_params.Bandwidth/1e6,...
    LTE_params.ChanMod_config.type,...
    LTE_params.UE_config.nTX,...
    LTE_params.BS_config.nRX,...
    N_subframes,...
    the_date(1),...      % Date: year
    the_date(2),...      % Date: month
    the_date(3),...      % Date: day
    the_date(4),...      % Date: hour
    the_date(5),...      % Date: minutes
    floor(the_date(6))); % Date: seconds)
    
if strcmp(LTE_params.channel_matrix_source,'generated') && strcmp(LTE_params.channel_matrix_tracefile,'auto')
    LTE_params.channel_matrix_tracefile = fullfile('./results/',['H_trace_' output_filename_UL '.mat']);
end

%% Control that channel tracing and noise seeding is only done under non-parallel mode and only once (length(SNR_vec)==1)
if LTE_params.random_noise_seeding && strcmp(LTE_params.simulation_method,'parallel')
    error('Using noise seeding with parallel simulations is not allowed. Results may not be repeatable');
end
if LTE_params.store_channel_trace && strcmp(LTE_params.simulation_method,'parallel')
    error('Storing the channel trace and (specially the) noise seed is not allowed in parallel mode. Results may not be repeatable');
end
if LTE_params.store_channel_trace && length(SNR_vec)~=1
    error('Simulations that store the channel trace can consist only of an SNR vector of length one (1)');
end

% If in trace mode, load the channel matrix trace
if strcmp(LTE_params.channel_matrix_source,'trace')
    load(LTE_params.channel_matrix_tracefile);
    
    % Check trace consistency
    if channel_matrix_tracefile.bandwidth~=LTE_params.Bandwidth
        error('Loaded channel trace is at a different sampling frequency. %3.2f MHz bandwidth found, %3.2f MHz required',channel_matrix_tracefile.bandwidth/1e6,LTE_params.Bandwidth/1e6);
    elseif ~strcmp(channel_matrix_tracefile.type,LTE_params.ChanMod_config.type)
        error('Config file specifies %s channel. Trace is %s',LTE_params.ChanMod_config.type,channel_matrix_tracefile.type);
    elseif (channel_matrix_tracefile.nTx~=LTE_params.UE_config.nTX) || (channel_matrix_tracefile.nRx~=LTE_params.BS_config.nRX)
        error('Config file specifies %dx%d system. %dx%d found in trace',LTE_params.UE_config.nTX,LTE_params.BS_config.nRX,channel_matrix_tracefile.nTx,channel_matrix_tracefile.nRx);
    elseif channel_matrix_tracefile.TTI_length < N_subframes
        error('Trace (%d TTIs) is shorter than simulation (%d TTIs)',channel_matrix_tracefile.TTI_length,N_subframes);
    elseif ~strcmp(LTE_params.ChanMod_config.filtering,'BlockFading')
        error('Tracing only supported for Block Fading due to memory consumption issues');
    end
    
    LTE_params.channel_matrix_trace = channel_matrix_tracefile;
    
    % Seed the noise RandStream with the seed specified in the trace
    LTE_params.noise_RandStream = RandStream('mt19937ar','Seed',channel_matrix_tracefile.noise_seed);
    LTE_params.read_channel_from_trace = true;
else
    LTE_params.read_channel_from_trace = false;
end

%% Generate channel realization using Winner II Channel Model
switch ChanMod.type
    case 'winner_II'
        LTE_params.Arrays = LTE_winner_channel_model_antenna_init();
        [channel, delays, out] = LTE_UL_winner_channel_model(N_subframes,LTE_params.Arrays);
    case 'TR 36.873'
        cd('TR36873_3D_Model_standalone');
        % NOTE: Ugly workaround for multisim and multi-CQI simulations
        if exist('channel', 'var')
            save('temp_3D');
            clear('classes'); %#ok<CLCLS>
            load('temp_3D');
            channel_existed = true;
        else
            channel_existed = false;
        end
        glob_rand_stream = RandStream.setGlobalStream(LTE_params.channel_param_RandStream);
        [channel, pathlosses] = channel_matrix_UL_LL(N_subframes, LTE_params);
        RandStream.setGlobalStream(glob_rand_stream);
        cd('..');
        % NOTE: Ugly workaround for multisim and multi-CQI simulations
        if channel_existed
            save(fullfile('TR36873_3D_Model_standalone', 'temp_3D'));
            clear('classes'); %#ok<CLCLS>
            load(fullfile('TR36873_3D_Model_standalone', 'temp_3D'));
            delete(fullfile('TR36873_3D_Model_standalone', 'temp_3D.mat'));
        end
        out = [];  
    otherwise
        channel = [];
        out = [];
end
