function [recovered_Sig] = Recover_Delay_Phase_Noise(delay_phase_distorted_rxSig,rxSig_org)
% Works with 2SpS
avgstd = zeros(500,2);
for k=1:500
    phasedif = angle(circshift(delay_phase_distorted_rxSig(1:2:2*length(rxSig_org)), -k)./rxSig_org);
    avgstd(k,1) = mean(phasedif);
    avgstd(k,2) = std(phasedif);
end
% figure;
% plot(1:500, avgstd(:,1));
% title('AVG');
% figure;
% plot(1:500, avgstd(:,2));
% title('STD');
[M I] = min(avgstd(:,2));
recovered_Sig = delay_phase_distorted_rxSig(2*I+1:end).*exp(-1i*avgstd(I,1));