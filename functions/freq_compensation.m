function [X_CD_rec, Y_CD_rec] = freq_compensation(TX_Xpol, TX_Ypol, SIG_SpS, SIG_symbolRate)
% freq_compensation Compensates for frequency offset in received optical signals.
%
% Inputs:
%   TX_Xpol - Received signal for X polarization
%   TX_Ypol - Received signal for Y polarization
%   SIG_SpS - Samples per symbol
%   SIG_symbolRate - Symbol rate of the signal
%
% Outputs:
%   X_CD_rec - Frequency-compensated signal for X polarization
%   Y_CD_rec - Frequency-compensated signal for Y polarization

% Calculate the bandwidth and frequency axis
Bs = SIG_SpS * SIG_symbolRate; % System bandwidth
f = (-Bs/2:Bs/length(TX_Xpol):Bs/2 - Bs/length(TX_Xpol)); % Frequency axis centered at 0

% Calculate the time axis
T = 1 / SIG_symbolRate; % Sampling period
N = length(TX_Xpol); % Number of samples
t = (1:N) * T; % Time vector

% Frequency offset estimation and compensation for X polarization
y = TX_Xpol .^ 4; % Fourth power nonlinearity for frequency estimation
Y = fftshift(abs(fft(y))); % FFT and shift to center
f_est = abs((1/4) * f(Y == max(Y))) / SIG_SpS; % Estimate the frequency offset
if isnan(max(Y))
    f_est=0;
end

X_CD_rec = TX_Xpol .* exp(-1j * 2 * pi * f_est * t)'; % Compensate for the frequency offset


% Frequency offset estimation and compensation for Y polarization
y = TX_Ypol .^ 4; % Fourth power nonlinearity for frequency estimation
Y = fftshift(abs(fft(y))); % FFT and shift to center
f_est = abs((1/4) * f(Y == max(Y))) / SIG_SpS; % Estimate the frequency offset
if isnan(max(Y))
    f_est=0;
end

Y_CD_rec = TX_Ypol .* exp(-1j * 2 * pi * f_est * t)'; % Compensate for the frequency offset


end

