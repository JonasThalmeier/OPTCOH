function [wML] = MLFilterViterbi(M, Delta_nu, Rs, OSNRdB, Es, NPol, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MLFilterViterbi calculates the maximum likelihood (ML) filter for the
% Viterbi & Viterbi algorithm.
%
% The ML filter is calculated as:
%   wML = (1^T * C^-1)^T
% where (.)^T indicates transpose and C is a covariance matrix, which
% depends on the signal-to-noise ratio and on the phase noise magnitude.
%
% Inputs:
% - M: Modulation order of the M-PSK modulation format
% - Delta_nu: Sum of transmitter and local oscillator laser linewidths in Hz
% - Rs: Symbol rate in symbols/second
% - OSNRdB: Channel OSNR in dB
% - Es: Symbol energy (per polarization orientation) in W
% - NPol: Number of polarization orientations used
% - N: Number of past and future symbols used in the Viterbi & Viterbi 
%      algorithm for phase noise estimates. The block length is L = 2*N + 1
%
% Output:
% - wML: Maximum likelihood filter to be used in the Viterbi & Viterbi 
%        algorithm for phase noise estimation. wML is a column vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate block length, symbol period, and phase noise variance
L = 2 * N + 1;
Ts = 1 / Rs;
Sigma_DeltaTheta2 = 2 * pi * Delta_nu * Ts;

% Calculate additive noise variance
SNRLin = 10^(OSNRdB / 10) * (2 * 12.5e9) / (NPol * Rs);
Sigma_eta2 = Es / (2 * SNRLin);

% Initialize K matrix
KAux = zeros(N + 1);
K = zeros(L);
for i = 0:N
    for ii = 0:N
        KAux(i + 1, ii + 1) = min(i, ii);
    end
end

% Construct K matrix
K(1:N + 1, 1:N + 1) = rot90(KAux(1:N + 1, 1:N + 1), 2);
K(N + 1:L, N + 1:L) = KAux(1:N + 1, 1:N + 1);

% Identity matrix
I = eye(L);

% Obtain the covariance matrix
C = Es^M * M^2 * Sigma_DeltaTheta2 * K + Es^(M - 1) * M^2 * Sigma_eta2 * I;

% Calculate filter coefficients
wML = (ones(L, 1)' / C).';
wML = wML / max(wML);
end

function [v] = vit_n_vit(z, Delta_nu, Rs, OSNRdB, Es, NPol, M, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ViterbiCPR performs phase recovery of signals with M-PSK modulation
% formats using the Viterbi & Viterbi algorithm. The signal 'z' is first
% raised to the M-th power, and then ML filtering is performed before 
% phase estimation. The phase estimates are applied to a phase unwrapper,
% and finally, the phase noise is compensated.
%
% Inputs:
% - z: Signal in which phase recovery will be performed. For transmission
%      in single polarization orientation, 'z' must be a column vector.
%      For transmission with polarization multiplexing, 'z' must be a 
%      matrix with two column-oriented vectors, where each column vector
%      corresponds to the signal of one polarization orientation. Signal
%      'z' must be obtained at 1 sample per symbol and must be normalized
%      to unitary power.
% - Delta_nu: Sum of the transmitter and local oscillator laser linewidths
%             in Hz
% - Rs: Symbol rate in symbols/second
% - OSNRdB: Channel OSNR in dB
% - Es: Symbol energy (per polarization orientation) in W
% - NPol: Number of polarization orientations used
% - M: Modulation order of the M-PSK modulation format
% - N: Number of past and future symbols used in the Viterbi & Viterbi 
%      algorithm for phase noise estimation. The block length is then
%      'L = 2*N + 1'.
%
% Outputs:
% - v: Signal produced after compensating for the phase noise present in 'z'
% - varargout: When the flag 'ParamViterbi.PEstimate' is 'true', the estimated
%              phase noise is also an output of the function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate block length for phase estimation
L = 2 * N + 1;

% Calculate ML filter
wML = MLFilterViterbi(M, Delta_nu, Rs, OSNRdB, Es, NPol, N);

% Initialize the vector of phase estimates
ThetaML = zeros(size(z, 1), NPol);

% Phase noise estimation
for Pol = 1:NPol
    % Prepare input blocks
    zBlocks = [zeros(floor(L / 2), 1); z(:, Pol); zeros(floor(L / 2), 1)];
    zBlocks = convmtx(zBlocks.', L);
    zBlocks = flipud(zBlocks(:, L:end-L+1));
    
    % Generate phase estimates
    ThetaML(:, Pol) = (1/M) * angle(wML.' * (zBlocks.^M)) - pi/M;
end
clearvars zBlocks;

% Vector of phase estimates after phase unwrapping
ThetaPU = zeros(size(ThetaML, 1), NPol);

% Initial 'previous phase' for unwrapping operation
ThetaPrev = zeros(1, NPol);

% Phase unwrapping
for i = 1:size(ThetaML, 1)
    % Phase unwrapper
    n = floor(1/2 + (ThetaPrev - ThetaML(i, :)) / (2 * pi / M));
    ThetaPU(i, :) = ThetaML(i, :) + n * (2 * pi / M);
    ThetaPrev = ThetaPU(i, :);
end

% Phase noise compensation
v = z .* exp(-1i * ThetaPU);

end
