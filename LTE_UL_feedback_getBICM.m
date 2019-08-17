function BICM = LTE_UL_feedback_getBICM(LTE_params,SNR)
% this function delivers the BICM in bits per channel use for 4/16/64 QAM
% and the SNR values handed over
% author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at

BICM_temp = 0;

for i2 = 1:3
    BICM_temp = 0;
    for i = 1:length(SNR)
        [~,SNR_temp] = min(abs(LTE_params.MI_data(i2).SNR-SNR(i)));
        BICM_temp = BICM_temp+LTE_params.MI_data(i2).BICM(SNR_temp);
    end
    BICM_temp_temp(i2) = BICM_temp;
end
BICM = max(BICM_temp_temp);