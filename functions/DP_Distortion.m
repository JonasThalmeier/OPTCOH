function [delay_phase_distorted_RX_Xpol,delay_phase_distorted_RX_Ypol] = DP_Distortion(TX_Xpol,TX_Ypol)
%Performs delay, phase interferences and convolution

%delay = randi(floor(length(TX_Xpol)/160),1); % maximum delay of half a period
%delay = randi(400,1);
delay = 8 * randi(8*50,1);

fprintf('The random delay introduced is (x8): %d\n', delay);

phase = randi(361,1);
fprintf('The random phase introduced is (degrees): %d\n', (phase-1));
phase = (phase-1) *pi / 180; %radians

delay_phase_distorted_RX_Xpol = [zeros(delay, 1, 'like', TX_Xpol)', TX_Xpol']' .* exp(1i * phase); %add zeros at beginning to simulate delay
delay_phase_distorted_RX_Ypol = [zeros(delay, 1, 'like', TX_Ypol)', TX_Ypol']' .* exp(1i * phase);

end