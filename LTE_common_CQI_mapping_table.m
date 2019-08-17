function CQI = LTE_common_CQI_mapping_table(CQI_mapping_params,SINRs,CQIs)
% The function takes as input argument a number of SINRs corresponding
% to the input set of CQIs and delivers as output the highest possible CQI 
% with BLER <= 0.1

temp = zeros(size(CQIs));
temp(CQI_mapping_params.table(CQIs) <= SINRs) = 1;
CQI = find(temp,1,'last')-1;


CQI(CQI == 0) = 20;  