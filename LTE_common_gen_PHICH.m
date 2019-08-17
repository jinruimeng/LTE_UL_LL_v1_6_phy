function [PHICHmapping,UsedElements,NElements,k_sort] = LTE_common_gen_PHICH(Nrb,Nsc,Nsub,Ns,NoData,NIDcell,BS,RefSig,LTE_params,b_,subframe_num)
% Reserves space for the PHICH 
% TS 36.211 V8.9.0, Section 6.9.3  
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
% The PHICH mapping assumes only 1 PHICH group and non MBSFN !!!!!!!!!!!!!!

%% Mapping to resource elements
% TS 36.211 V8.9.0 Section 6.9.3
% 1)
for l = 1:1 % only 1st OFDM symbol is assumed to be used 
    % 2)
    % number of resource elemenet groups not assigned to PCFICH
    UnUsedElementsPCFICH = find(LTE_params.PCFICH(b_,subframe_num).Used_REgroups(:,1)==0);
    if l == 1
        n_l0 = length(UnUsedElementsPCFICH);
        n_l = n_l0;
    else
        disp('error')
    end
    % 3)
    % 4)
    m = 0; % PHICH group number
    % 5)
    for i = 1:3
        % 7)
        li = 0;
        % NOTE there are other 3 cases - MBSFN, frame structure type 2 and
        % other which are not implemented
        % 8)
        switch i
            case 1
                ni = mod((floor(BS.NIDcell*n_l/n_l0)+m),n_l);
            case 2
                ni = mod((floor(BS.NIDcell*n_l/n_l0)+m+floor(n_l/3)),n_l);
            case 3
                ni = mod((floor(BS.NIDcell*n_l/n_l0)+m+floor(2*n_l/3)),n_l);
            otherwise
                error('Error (LTE_PHICH)')
        end
        m = m + 1;
        if l ==1
            UsedElements_PHICH(i,1) = UnUsedElementsPCFICH(ni+1,1);
        else
            error('This combination is not implmented yet (LTE-common-gen-PHICH)')
        end
    end
end

% Used RE groups -> used frequences (example: [5; 9; 2] -> [31; 55; 13;])
k = zeros(3,1);
for i1=1:3;
k(i1,1) = (UsedElements_PHICH(i1,1)-1)*6+1; 
end

% help variable - it enables to map n-th quadruplet to the right (n-th) 
% resource element group (example: k = [31; 55; 13;] -> k_sort = [3; 1; 2;])
k_help = sort(k);
k_sort = zeros(3,1);
for i0 = 1:3
    k_sort(i0,1) = find(k_help(i0,1) == k(:,1)); 
end

%% Own mapping
PHICHmapping = zeros(Nrb*Nsc,Nsub,BS.nAtPort);
% REs occupied by reference signal have to be omitted
for i2 = 1:BS.nAtPort
  PHICHmap_temp = PHICHmapping(:,:,i2); 
  for i3 = 1:3 % 3 resource element groups
    kp = k(i3,1);
    RE_assign = 0;
    while RE_assign < 4 % 4 resource elements in group
        if NoData(kp,1) == 1; % the resource element is occupied by reference signal
            kp = kp + 1;
        else
           PHICHmap_temp(kp,1) = 1;
           RE_assign = RE_assign+1;
           kp = kp + 1;
        end
    end
  end
  PHICHmapping(:,:,i2) = PHICHmap_temp;
end

%% Image and control

% check collisions between PHICH, PCFICH and reference signals
if sum(sum((sum(PHICHmapping,3)).*NoData)) ~= 0;
    error('Wrong PHICH mapping 1 (LTE-common-gen-PHICH)')
end
if sum(sum((sum(PHICHmapping(:,:,1),3)).*LTE_params.PCFICH(1,1).Mapping(:,:,1))) ~= 0;
    error('Wrong PHICH mapping 2 (LTE-common-gen-PHICH)')
end
if sum(PHICHmapping(:,1,1)) ~= 12;
    error('Wrong PHICH mapping 3 (LTE-common-gen-PHICH)')
end

%% Used elements, NEelements
UsedElements = zeros(Nrb*2,BS.nAtPort);
for i5 = 1:BS.nAtPort
    for i6 = 1:2
        for i7 = 1:Nrb
            PHICHmapping_temp = PHICHmapping(:,:,i5);
            UsedElements_temp(i7+(i6-1)*Nrb) = sum(sum(PHICHmapping_temp((i7-1)*Nsc+1:i7*Nsc,(i6-1)*Ns+1:Ns*i6)));
        end
    end
    UsedElements(:,i5) = UsedElements_temp;
    NElements(:,i5) = sum(sum(UsedElements(:,i5)));
end


