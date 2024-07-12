function [TX_Xpol, TX_Ypol] = DP_Distortion(TX_Xpol, TX_Ypol, delta_nu, rad_sec, SIG_symbolRate, f_offset)
% DP_Distortion Applies delay, phase interferences, and convolution to the input signals.
% 
% Inputs:
%   TX_Xpol - Input signal for X polarization
%   TX_Ypol - Input signal for Y polarization
%   delta_nu - Linewidth of the laser (related to phase noise)
%   rad_sec - Variance of the polarization rotation process
%   SIG_symbolRate - Symbol rate of the signal
%   f_offset - Frequency offset
%
% Outputs:
%   TX_Xpol - Distorted signal for X polarization
%   TX_Ypol - Distorted signal for Y polarization

% Generate a random delay up to half a period
delay = randi(floor(length(TX_Xpol) / 160), 1);

% -------------------- Phase Rotation -------------------------------------
Theta_0 = randi([0,360]) *pi / 180; % Initial phase offset
var_Theta = 2 * pi * delta_nu; % Variance of phase noise
samp_rate = 8 * SIG_symbolRate; % Sampling rate

lenX = length(TX_Xpol) + delay; % Length of the signal after delay

% Generate random steps for the Wiener process (phase noise)
Delta_Theta = sqrt(var_Theta / samp_rate) * randn(1, lenX);
Theta_n = cumsum(Delta_Theta) + Theta_0; % Cumulative sum to get phase noise

% Rotate the constellation by applying the phase noise
TX_Xpol = [zeros(delay, 1, 'like', TX_Xpol)', TX_Xpol']' .* exp(1i * Theta_n)'; % Add zeros at beginning to simulate delay
TX_Ypol = [zeros(delay, 1, 'like', TX_Ypol)', TX_Ypol']' .* exp(1i * Theta_n)';

% -------------------- Frequency Impairment ------------------------------
N = length(TX_Xpol); % Length of the signal
T = 1 / samp_rate; % Sampling period
t = (1:N) * T; % Time vector

% Apply frequency offset to the signals
TX_Xpol = TX_Xpol .* exp(1j * 2 * pi * f_offset * t)';
TX_Ypol = TX_Ypol .* exp(1j * 2 * pi * f_offset * t)';

% ------------------ Jones Matrix (Polarization Rotation) -----------------
Phi_0 = 0; % Initial phase for polarization rotation
Delta_Phi = sqrt(rad_sec / samp_rate) * randn(1, lenX); % Random steps for polarization rotation
Phi_n = cumsum(Delta_Phi) + Phi_0; % Cumulative sum to get polarization rotation

amplitude = 1; % Amplitude of the Jones matrix components
phi = pi / 3; % Phase shift for the off-diagonal components

% Define the Jones matrix components
pol_xx = amplitude * cos(Phi_n);
pol_xy = exp(-1i * phi) * amplitude * sin(Phi_n);
pol_yx = -exp(1i * phi) * amplitude * sin(Phi_n);
pol_yy = amplitude * cos(Phi_n);

% Apply the Jones matrix to the signals
TX_Xpol = pol_xx.' .* TX_Xpol + pol_xy.' .* TX_Ypol;
TX_Ypol = pol_yx.' .* TX_Xpol + pol_yy.' .* TX_Ypol;
