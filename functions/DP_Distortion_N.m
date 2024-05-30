function [delay_phase_distorted_RX_Xpol,delay_phase_distorted_RX_Ypol] = DP_Distortion_N(TX_Xpol,TX_Ypol)
%Performs delay, phase interferences and convolution

delay = randi(floor(length(TX_Xpol)/160),1); % maximum delay of half a period
fprintf('The random delay introduced is (x8): %d\n', delay);

phase_init = randi(361,1);
fprintf('The random phase introduced is (degrees): %d\n', (phase_init-1));
phase_init = (phase_init-1) *pi / 180; %radians

delay_phase_distorted_RX_Xpol = [zeros(delay, 1, 'like', TX_Xpol)', TX_Xpol']' .* exp(1i * phase_init); %add zeros at beginning to simulate delay
delay_phase_distorted_RX_Ypol = [zeros(delay, 1, 'like', TX_Ypol)', TX_Ypol']' .* exp(1i * phase_init);

end