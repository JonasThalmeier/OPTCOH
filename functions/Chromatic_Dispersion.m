function [Xpol_CD, Ypol_CD] = Chromatic_Dispersion(Xpol, Ypol, SpS, flag)
% Chromatic_Dispersion Simulates the effect of chromatic dispersion on optical signals.
%
% Inputs:
%   Xpol - Input signal for X polarization
%   Ypol - Input signal for Y polarization
%   SpS - Samples per symbol
%   flag - Determines the direction of dispersion compensation
%          1: Apply chromatic dispersion
%          2: Compensate for chromatic dispersion
%
% Outputs:
%   Xpol_CD - Output signal for X polarization after dispersion
%   Ypol_CD - Output signal for Y polarization after dispersion

% Constants
f0 = 0; % Center frequency offset
fc = 191; % Center frequency for C-band [THz] (approximately 1570 nm)
c = 3e8; % Speed of light in vacuum (m/s)
D = 17; % Dispersion parameter in ps/(nm*km)
L = 100; % Length of the fiber in km

% Calculate the dispersion coefficient beta_2 [ps^2/km]
beta_2 = -(D * c * 1e-3) / (2 * pi * (fc)^2); % Convert D to beta_2

% Bandwidth and frequency axis
Bs = SpS * 64e9; % System bandwidth
f = (-Bs/2:Bs/length(Xpol):Bs/2 - Bs/length(Xpol)).'; % Frequency axis
f = f / 1e11 * 1e-1; % Normalize the frequency axis

% Apply or compensate chromatic dispersion based on the flag
if flag == 1
    % Apply chromatic dispersion
    H_cd = exp(-1i * 2 * pi^2 * beta_2 * (f - f0).^2 * L); % Transfer function for dispersion
    Xpol_CD_fft = fft(Xpol) .* H_cd; % Apply dispersion in the frequency domain
    Ypol_CD_fft = fft(Ypol) .* H_cd;
    Xpol_CD = ifft(Xpol_CD_fft); % Transform back to the time domain
    Ypol_CD = ifft(Ypol_CD_fft);
    
elseif flag == 2
    % Compensate for chromatic dispersion
    H_cd = exp(1i * 2 * pi^2 * beta_2 * (f - f0).^2 * L); % Inverse transfer function for dispersion
    Xpol_CD_fft = fft(Xpol) .* H_cd; % Apply compensation in the frequency domain
    Ypol_CD_fft = fft(Ypol) .* H_cd;
    Xpol_CD = ifft(Xpol_CD_fft); % Transform back to the time domain
    Ypol_CD = ifft(Ypol_CD_fft);
end
end
