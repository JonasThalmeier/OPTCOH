function [TX_noise, Noise] = WGN_Noise_Generation(TX, SpS, M, EbN0)
% WGN_Noise_Generation Adds white Gaussian noise (WGN) to the input signal.
%
% Inputs:
%   delay_phase_distorted_RX_Apol - Input signal with delay and phase distortion
%   SpS - Samples per symbol
%   M - Modulation order (e.g., 4 for QPSK, 16 for 16QAM)
%   EbN0 - Energy per bit to noise power spectral density ratio (dB)
%   Rs - Symbol rate (unused in this function)
%
% Outputs:
%   delay_phase_noise_distorted_RX_Apol - Output signal with added WGN
%   Noise - The generated noise that was added to the signal

% Number of bits per symbol
Nbit = log2(M);

% Convert Eb/N0 from dB to linear scale
EbN0_lin = 10^(EbN0/10);

% Number of samples in the input signal
Nsamples = length(TX);

% Evaluate the signal power
SignalPower = mean(abs(TX).^2);

% Evaluate the noise power
% E_s = SignalPower * SpS (Energy per symbol)
% E_b = E_s / Nbit (Energy per bit)
% N_0 = E_b / EbN0_lin (Noise power spectral density)
% Noise variance = N_0 / 2 (for complex noise)
NoisePower = SignalPower * SpS / (EbN0_lin);

% Generate normalized complex WGN (mean 0, variance 1)
NoiseNormalized = (1/sqrt(2)) * (randn(1, Nsamples) + 1i * randn(1, Nsamples));

% Scale the noise to the desired noise power
Noise = NoiseNormalized * sqrt(NoisePower/4);

% Add the noise to the input signal
TX_noise = TX + Noise.';

end
