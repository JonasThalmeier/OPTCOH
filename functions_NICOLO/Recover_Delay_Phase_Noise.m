function [recovered_TXSig_Apol] = Recover_Delay_Phase_Noise(rxSig_Apol,delay_phase_distorted_RX_Apol)

% ----- Recover from delay and phase at 8 SpS, this way no problems with downsampling

% Delay
[corr1, lag_1] = xcorr(rxSig_Apol,delay_phase_distorted_RX_Apol);
figure(), stem(lag_1,abs(corr1));
[max_corr, max_index] = max(abs(corr1));
fprintf('The Xpol tracked delay is of %d samples.\n', abs(lag_1(max_index)));
recovered_TXSig_Apol = delay_phase_distorted_RX_Apol(abs(lag_1(max_index))+1:end);
[corr1, lag_1] = xcorr(rxSig_Apol,recovered_TXSig_Apol);
figure(), stem(lag_1,abs(corr1));

% Phase 

[theta,rho] = cart2pol(real(rxSig_Apol),imag(rxSig_Apol));
[theta2,rho] = cart2pol(real(recovered_TXSig_Apol),imag(recovered_TXSig_Apol));
angle = theta2-theta;
fprintf('The Xpol tracked phase is of %.0f.\n', (mean(angle)*180/pi));

recovered_TXSig_Apol = recovered_TXSig_Apol.*exp(-1i .* angle);

% 
% Ki = 0.01; % learning rate
% N = length(recovered_TXSig_Apol);
% phase_estimate = zeros(1, N);
% phase_estimate1 = zeros(1, N);
% phase_estimate(1) = angle(rxSig_Apol(1));
% phase_estimate1(1) = angle(recovered_TXSig_Apol(1));
% %phase_corrected_signal = zeros(size(rxSig_Apol));
% 
% for k = 2:length(recovered_TXSig_Apol)
%     % Phase detector
%     phase_estimate(k) = phase_estimate(k-1) + Ki.*imag(rxSig_Apol(k).*exp(-1i*phase_estimate(k-1)));
%     phase_estimate1(k) = phase_estimate1(k-1) + Ki.*imag(recovered_TXSig_Apol(k).*exp(-1i*phase_estimate1(k-1)));
%     
% end
% 
% phase_compensation = phase_estimate1 - phase_estimate; 
% recovered_TXSig_Apol = recovered_TXSig_Apol.*exp(-1i.*phase_compensation'); 

fprintf('The Xpol sequence is correctly recovered [1 yes/ 0 no]: %d\n', isequal(round(rxSig_Apol,8),round(recovered_TXSig_Apol,8)));
% we have to add this round because when it comes to compare to the distorted version matlab doesn't set them equal due to to high precision

end