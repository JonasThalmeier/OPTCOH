function [recovered_Sig] = Recover_Delay_Phase_Noise(delay_phase_distorted_rxSig,rxSig_org)
% Works with 2SpS
N = floor(length(rxSig_org)/160);
stdarray = zeros(N,1);
for k=1:N
    phasedif = angle(circshift(delay_phase_distorted_rxSig(1:2:2*length(rxSig_org)), -k)./rxSig_org);
    avgstd(k) = std(phasedif);
end
[M I] = min(avgstd(:,2));
avg = mean(angle(circshift(delay_phase_distorted_rxSig(1:2:2*length(rxSig_org)), -I)./rxSig_org);)
recovered_Sig = delay_phase_distorted_rxSig(2*I+1:end).*exp(-1i*avg);