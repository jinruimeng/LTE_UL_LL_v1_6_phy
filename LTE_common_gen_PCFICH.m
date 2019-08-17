function [PCFICHmapping,UsedElements,NElements,k_sort,UsedREgroups] = LTE_common_gen_PCFICH(Nrb,Nsc,Nsub,Ns,NoData,NIDcell,BS,RefSig)
% Reserves space for the PCFICH 
% TS 36.211 V8.9.0, Section 6.7.4  
% Author: Petr Kej�k, xkejik00@stud.feec.vutbr.cz
% 2010
% 
% input :   Nrb       ... [1 x 1]double number of resource blocks
%           Nsc       ... [1 x 1]double number - resource block size
%           Nsub      ... [1 x 1]double number of OFDM symbols in subframe
%           Ns        ... [1 x 1]double number of OFDM symbols in slot
%           NoData    ... reserved space for synchronization signal
% output:   k_sort    ... [4 x 1] double - defines order of used resource
%                         element groups

%% Warnings, comments, possible changes
% The PCFICH shall be transmitted on the same antenna ports as PBCH.
% The PCFICH shall be transmitted only when the number of OFDM symbols for PDCCH is greater than zero.
% help variable and sorting ?

%% Mapping to resource elements
PCFICHmapping = zeros(Nrb*Nsc,Nsub,BS.nAtPort);

% NIDcell=3; % test case
k_add = (Nsc/2)*mod(NIDcell,2*Nrb); % addition
k = zeros(4,1);
k(1,1) = mod(k_add,Nrb*Nsc) + 1;                      
k(2,1) = mod(k_add + floor(Nrb/2)*(Nsc/2),Nrb*Nsc) + 1;    
k(3,1) = mod(k_add + floor(2*Nrb/2)*(Nsc/2),Nrb*Nsc) + 1;  
k(4,1) = mod(k_add + floor(3*Nrb/2)*(Nsc/2),Nrb*Nsc) + 1; 
% + 1 -- corrected to matlab indices (zero excluded) 

% help variable - it enables to map n-th quadruplet to the right (n-th) 
% resource element group (example: k = [43; 61; 7; 25] -> k_sort = [3; 4; 1; 2])
k_help = sort(k);
k_sort = zeros(4,1);
for i0 = 1:4
    k_sort(i0,1) = find(k_help(i0,1) == k(:,1)); 
end

%% Own mapping
PCFICHmapping = zeros(Nrb*Nsc,Nsub,BS.nAtPort);
% REs occupied by reference signal have to be omitted
for i2 = 1:BS.nAtPort
  PCFICHmap_temp = PCFICHmapping(:,:,i2); 
  for i3 = 1:4 % 4 resource element groups
    kp = k(i3,1);
    RE_assign = 0;
    while RE_assign < 4 % 4 resource elements in grup
        if NoData(kp,1) == 1; % the resource element is occupied by reference signal
            kp = kp + 1;
        else
           PCFICHmap_temp(kp,1) = 1;
           RE_assign = RE_assign+1;
           kp = kp + 1;
        end
    end
  end
  PCFICHmapping(:,:,i2) = PCFICHmap_temp;
end

%% Image and control

% check collisions between PCFICH and reference signals
if sum(sum((sum(PCFICHmapping,3)).*NoData)) ~= 0;
    error('Wrong PCFICH mapping (LTE-common-gen-PCFICH)')
end
if sum(PCFICHmapping(:,1,1)) ~= 16;
    error('Wrong PCFICH mapping (LTE-common-gen-PCFICH)')
end

%% UsedElements
UsedElements = zeros(Nrb*2,BS.nAtPort);
for i5 = 1:BS.nAtPort
    for i6 = 1:2
        for i7 = 1:Nrb
            PCFICHmapping_temp = PCFICHmapping(:,:,i5);
            UsedElements_temp(i7+(i6-1)*Nrb) = sum(sum(PCFICHmapping_temp((i7-1)*Nsc+1:i7*Nsc,(i6-1)*Ns+1:Ns*i6)));
        end
    end
    UsedElements(:,i5) = UsedElements_temp;
    NElements(:,i5) = sum(sum(UsedElements(:,i5)));
end

UsedREgroups = zeros(Nrb*2,BS.nAtPort); % used RE groups in the first OFDM symbol
for i8 = 1:BS.nAtPort
    for i9 = 1:2*Nrb % 2 RE groups per 12 subcarriers
            PCFICHmapping_temp = PCFICHmapping(:,:,i8);
            UsedREgroups_temp(i9,1) = sum(sum(PCFICHmapping_temp((i9-1)*6+1:(i9-1)*6+6)));
        
    end
    UsedREgroups(:,i8) = UsedREgroups_temp;
end
 