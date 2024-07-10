function [delay_phase_distorted_RX_Xpol,delay_phase_distorted_RX_Ypol] = DP_Distortion(TX_Xpol,TX_Ypol, delta_nu, rad_sec, SIG_symbolRate)
%Performs delay, phase interferences and convolution

delay = randi(floor(length(TX_Xpol)/160),1); % maximum delay of half a period

fprintf('The random delay introduced is (x8): %d\n', delay);

phase_init = randi(361,1);
fprintf('The random phase introduced is (degrees): %d\n', (phase_init-1));
phase_init = (phase_init-1) *pi / 180; %radians

% delay_phase_distorted_RX_Xpol = [zeros(delay, 1, 'like', TX_Xpol)', TX_Xpol']' .* exp(1i * phase_init); %add zeros at beginning to simulate delay
% delay_phase_distorted_RX_Ypol = [zeros(delay, 1, 'like', TX_Ypol)', TX_Ypol']' .* exp(1i * phase_init);

%-------------Wiener Process/Random Walk-----------------------------------
% delta_nu = 50e3;  % Laser linewidth. 50kHz seems to be realisitic
var = 2*pi*delta_nu;
%var = var;
samp_rate = 8*SIG_symbolRate; % deltaW of a Wiener process have a variance equal to the stime between steps
% n = 1000;
lenX = length(TX_Xpol)+delay;
% lenY = length(TX_Ypol)+delay;

% Generating the random steps
steps = sqrt(var/samp_rate)*randn(1,lenX);
phase = cumsum(steps)+phase_init;


% Rotate constellation
delay_phase_distorted_RX_Xpol_1 = [zeros(delay, 1, 'like', TX_Xpol)', TX_Xpol']' .* exp(1i * phase)'; %add zeros at beginning to simulate delay
delay_phase_distorted_RX_Ypol_1 = [zeros(delay, 1, 'like', TX_Ypol)', TX_Ypol']' .* exp(1i * phase)';

%--------------------frequency impairment-----------------
N = length(delay_phase_distorted_RX_Xpol_1);            
%after 6e9 f_offset stops work
f_offset = 1e9;      % Frequency offset (1 GHz)
T = 1/SIG_symbolRate; 
 
t = (1:N)* T;
delay_phase_distorted_RX_Xpol_1 = delay_phase_distorted_RX_Xpol_1 .* exp(1j*2*pi*f_offset*t)';
delay_phase_distorted_RX_Ypol_1 = delay_phase_distorted_RX_Ypol_1 .* exp(1j*2*pi*f_offset*t)';

% ------------------Jones Matrix (Pol.rotation)-----------------------------
% How to handle lenX~=lenY????


% kappa = 0;
% rad_sec = 1e6;
Theta = randi([0,0]) *pi / 180;
stepsPol = sqrt(rad_sec/samp_rate)*randn(1,lenX);
phasePol = cumsum(stepsPol)+Theta;


amplitude = 1;
phi = pi/3;

pol_xx = amplitude * cos(phasePol);
pol_xy = exp(-1i*phi)*amplitude * sin(phasePol);
pol_yx = -exp(1i*phi)*amplitude * sin(phasePol);
pol_yy = amplitude * cos(phasePol);


delay_phase_distorted_RX_Xpol = pol_xx.'.*delay_phase_distorted_RX_Xpol_1 + pol_xy.'.*delay_phase_distorted_RX_Ypol_1;
delay_phase_distorted_RX_Ypol = pol_yx.'.*delay_phase_distorted_RX_Xpol_1 + pol_yy.'.*delay_phase_distorted_RX_Ypol_1;

% TX_Xpol = [zeros(delay, 1, 'like', TX_Xpol)', TX_Xpol']';
% TX_Ypol = [zeros(delay, 1, 'like', TX_Ypol)', TX_Ypol']';
% J = zeros(2, 2, lenX);
% R = zeros(2, 2, lenX);
% JR = zeros(2, 2, lenX);
% for idx=1:lenX
%     J(:,:,idx) = [exp(1i * phase(idx))', kappa;kappa,exp(1i * phase(idx))];
%     R(:,:,idx) = [cos(phasePol(idx)),-sin(phasePol(idx));sin(phasePol(idx)),cos(phasePol(idx))];
%     JR(:,:,idx) = R(:,:,idx) * J(:,:,idx);
%     dist = [TX_Xpol(idx),TX_Ypol(idx)]*JR(:,:,idx);
%     delay_phase_distorted_RX_Xpol(idx) = dist(1);
%     delay_phase_distorted_RX_Ypol(idx) = dist(2);
% end
% end