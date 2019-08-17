function [ allocation ] = LTE_UL_allocate_RS(LTE_params, UE_MCS_and_scheduling_info, UE, subframe_i, subframe_corr, srs_subframe)
% RS allocation in a subframe for one UE
% 
% Stefan Pratschner, stefan.pratschner@nt.tuwien.ac.at
% www.nt.tuwien.ac.at

allocation = false(1,3);
[~,indy] = find(UE_MCS_and_scheduling_info.UE_mapping);

SRS = false;
if strcmp(UE.SRS_duration,'periodic') && srs_subframe
    switch LTE_params.frameStructure
        case 1  % FDD
            SRS = (mod(subframe_i - 1 - UE.Toffset,UE.Tsrs) == 0);
        case 2  % TDD
            if mod(subframe_corr-1,5) ~= 0
                if UE.Tsrs == 2
                    SRS = (mod(subframe_i - 1 - UE.Toffset,5) == 0);
                else
                    SRS = (mod(subframe_i - 1 - UE.Toffset,UE.Tsrs) == 0);
                end

                if LTE_params.UpPTS_length == 2 && mod(subframe_corr-1,5) == 1
                    SRS = [SRS (mod(subframe_i - 2 - UE.Toffset,UE.Tsrs) == 0)];
                end

                SRS = max(SRS);
            end
        otherwise
            error('frame structure not supported');
    end
end

% SRS allocation
if SRS;
    allocation([indy; 3]) = true; % SRS and DMRS allocated
else
    allocation(indy) = true;      % only DMRS (in every scheduled slot)               
end

end

