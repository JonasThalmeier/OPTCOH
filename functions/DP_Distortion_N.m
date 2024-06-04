function [delay_phase_distorted_RX_Xpol,delay_phase_distorted_RX_Ypol] = DP_Distortion_N(TX_Xpol,TX_Ypol)
%Performs delay, phase interferences and convolution

delay = randi(floor(length(TX_Xpol)/160),1); % maximum delay of half a period
fprintf('The random delay introduced is (x8): %d\n', delay);

phase_init = randi(361,1);
fprintf('The random phase introduced is (degrees): %d\n', (phase_init-1));
phase_init = (phase_init-1) *pi / 180; %radians

%-------------Wiener Process/Random Walk-----------------------------------
delta_nu = 50e3;  % Laser linewidth. 50kHz seems to be realisitic
var = 2*pi*delta_nu;
samp_rate = 58*64e9; % deltaW of a Wiener process have a variance equal to the stime between steps
n = 1000;
lenX = length(TX_Xpol)+delay;
lenY = length(TX_Ypol)+delay;

% Generating the random steps
stepsX = sqrt(var/samp_rate)*randn(1,lenX);
phaseX = cumsum(stepsX)+phase_init;
stepsY = sqrt(var/samp_rate)*randn(1,lenY);
phaseY = cumsum(stepsY)+phase_init;



% Rotate constellation
delay_phase_distorted_RX_Xpol = [zeros(delay, 1, 'like', TX_Xpol)', TX_Xpol']' .* exp(1i * phaseX)'; %add zeros at beginning to simulate delay
delay_phase_distorted_RX_Ypol = [zeros(delay, 1, 'like', TX_Ypol)', TX_Ypol']' .* exp(1i * phaseY)';

% delay_phase_distorted_RX_Xpol = [zeros(delay, 1, 'like', TX_Xpol)', TX_Xpol']' .* exp(1i * phase_init); %add zeros at beginning to simulate delay
% delay_phase_distorted_RX_Ypol = [zeros(delay, 1, 'like', TX_Ypol)', TX_Ypol']' .* exp(1i * phase_init);

end