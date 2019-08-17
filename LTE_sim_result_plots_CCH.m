function LTE_sim_result_plots_CCH(simulation_results)
% Plots BER for PCFICH channel
% Author: Petr Kejík, xkejik00@stud.feec.vutbr.cz
% 2010

%% Figure settings
SNR_vector = simulation_results.SNR_vector;
minimum_for_semilogy = 10^-4;

%% BER for PCFICH and PCHICH
figure();
semilogy(SNR_vector,sum(simulation_results.UE_specific.BER_PCFICH_CFI(:,:),1)/(length(simulation_results.UE_specific.BER_PCFICH_CFI(:,:))*32),'b');
hold on
semilogy(SNR_vector,sum(simulation_results.UE_specific.BER_PCFICH(:,:),1)/(length(simulation_results.UE_specific.BER_PCFICH_CFI(:,:))),'r');
for i=1:1
semilogy(SNR_vector,sum(simulation_results.UE_specific.BER_PHICH(:,:,i),1)/(length(simulation_results.UE_specific.BER_PHICH(:,:,i))),'g');
end
hold off
legend('uncoded BER (before descrambling)','BER (final CFI value)','PHICH');
xlabel('SNR [dB]');
ylabel('BER');
title('PCFICH BER');
axis([min(SNR_vector) max(SNR_vector),minimum_for_semilogy 10^0])
grid on
for i=1:1
sum(simulation_results.UE_specific.BER_PHICH(:,:,i),1)/(length(simulation_results.UE_specific.BER_PHICH(:,:,i)))
end

% %% BER for PHICH
% figure();
% semilogy(SNR_vector,sum(simulation_results.UE_specific.BER_PCFICH_CFI(:,:),1)/(length(simulation_results.UE_specific.BER_PCFICH_CFI(:,:))*32),'b');
% hold on
% semilogy(SNR_vector,sum(simulation_results.UE_specific.BER_PCFICH(:,:),1)/(length(simulation_results.UE_specific.BER_PCFICH_CFI(:,:))),'r');
% hold off
% legend('uncoded BER (before descrambling)','BER (final CFI value)');
% xlabel('SNR [dB]');
% ylabel('BER');
% title('PCFICH BER');
% axis([min(SNR_vector) max(SNR_vector),minimum_for_semilogy 10^0])
% grid on
    