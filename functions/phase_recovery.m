function outputSig = phase_recovery(inputSig, revSig)
%PHASE_RECOVERY Recovers the phase shift between a reference and a received signal.
% This function estimates the phase shift introduced to a signal during transmission
% by comparing the received signal with a reference signal. It adjusts the phase
% of the input signal to align it with the received signal.
%
% Inputs:
%   inputSig - The original transmitted signal (reference).
%   revSig   - The received signal which may have undergone a phase shift.
%
% Output:
%   outputSig - The phase-adjusted version of the input signal, aligned as closely
%               as possible with the received signal.
%
% The function first downsamples the input signal to match the length of the received
% signal. It then computes the cross-correlation between the received signal and
% multiple phase-shifted versions of the input signal. The phase shift yielding
% the highest correlation indicates the estimated phase shift. The output signal
% is the input signal adjusted by this estimated phase shift.

% Validate lengths and preprocess signals
len = length(revSig);                     % Length of the received signal
inputSig_1Sps = inputSig(2:2:len);             % Downsample the input signal to match the received signal length

% Initialize parameters
N = 100;                                  % Number of phase shifts to test
corr = zeros(N+1,1);                      % Preallocate correlation array

% Compute correlation for various phase shifts
for k = 0:N
    [c, lags] = xcorr(angle(inputSig_1Sps .* exp(1i*pi*k/50)), angle(revSig)); % Cross-correlation between phase-shifted input and received signal
    corr(k+1) = max(abs(c));                   % Store the maximum correlation for each phase shift
end
figure;
plot(corr);
% Find the optimal phase shift
[M, I] = max(corr);                       % Find the index of the highest correlation
outputSig = inputSig .* exp(1i*pi*(I-1)/50);  % Apply the optimal phase shift to the input signal
end
