clear;
close all;
clc;

% Load the .mat file
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "16-QAM"];
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
r = 2;
fprintf('The transmitted moduluation is: %s\n', modulation(r));
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_64GBaud.mat'));

% Parameters

SpS_down = 4;
SpS_up = SIG.Sps/SpS_down;

% txSig is now loaded along with other structures like SIG and PulseShaping

% I put here constellation estimation in order to generate properly noise
% since it depends on Nbit too. BUT HOW TO DO IT WITHOUT THE DOWNSAMPLED
% VERSION?

if r == 1
    M = 4;
else
    M = 16;
end

% Apply matched filtering
% Assuming b_coeff is your filter coefficients from the PulseShaping structure

rxSig_Xpol = conv(PulseShaping.b_coeff, SIG.Xpol.txSig);
rxSig_Ypol = conv(PulseShaping.b_coeff, SIG.Ypol.txSig);

% Create delay and phase convolved signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, PulseShaping.b_coeff);

% Adding the noise
[X_distorted_AGWN, NoiseX] = WGN_Noise_Generation(X_distorted,SIG.Sps, M, 15);
[Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_distorted,SIG.Sps, M, 15);

%----------------Downsample the signal---------------------------
X_2Sps = downsample(X_distorted_AGWN, SpS_down);
Y_2Sps = downsample(Y_distorted_AWGN, SpS_down);

%-----------------Recover the right sampling time-------------------------
[X_2Sps,Y_2Sps] = samp_time_recovery(X_2Sps,Y_2Sps,SIG.Sps);
% Plot constellation
figure;
scatter(real(X_2Sps(1:2:end)), imag(X_2Sps(1:2:end)), ".", "k");
title('Xpol constellation, samp time recovered');
grid on;

% ----- Recover from delay and phase -------------------------------------
X_2Sps = Recover_Delay_Phase_Noise(X_2Sps,SIG.Xpol.txSymb);
Y_2Sps = Recover_Delay_Phase_Noise(Y_2Sps,SIG.Ypol.txSymb);


% Plot constellation
figure;
scatter(real(X_2Sps(1:2:end)), imag(X_2Sps(1:2:end)), ".", "k");
title('Xpol constellation');
grid on;

figure;
scatter(real(Y_2Sps(1:2:end)), imag(Y_2Sps(1:2:end)), ".", "k");
title('Ypol constellation');
grid on;
%%
[counts, binEdges] = histcounts(angle(X_2Sps(1:2:end)), 12, 'Normalization', 'probability');
%histogram(angle(downsampledSig_Xpol(2:2:end)), 12, 'Normalization', 'probability');
if max(counts) > .17
    fprintf('The tracked moduluation is: QPSK\n');
    M = 2;
    [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QPSK_demapping(X_2Sps,Y_2Sps);
else
    fprintf('The tracked moduluation is: 16-QAM\n');
    M = 4;
    [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QAM_16_demapping(X_2Sps,Y_2Sps);
end

%----------------Majority voting over the repeated symbols-----------------
% Let's assume demappedBits_Xpol is your bit outcomes with size [N * Npp, 2]
% Where N is the number of unique symbols and Npp is the number of repetitions

% N = length(demappedBits_Xpol) / SIG.Npp;  % Calculate the number of unique symbols
N = length(SIG.Xpol.txSymb);  % Calculate the number of unique symbols

consolidatedBits_Xpol = zeros(N, size(SIG.Xpol.bits,2));
consolidatedBits_Ypol = zeros(N, size(SIG.Xpol.bits,2));

for i = 1:N
    consolidatedBits_Xpol(i, :) = mode(demappedBits_Xpol(i:N:end,:),1);
    consolidatedBits_Ypol(i, :) = mode(demappedBits_Ypol(i:N:end,:),1);
end
% Now, consolidatedBits contains the 'averaged' bit decisions



%-------------------------BER calculation----------------------------------

if size(consolidatedBits_Xpol) ~= size(SIG.Xpol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Xpol = sum(sum(consolidatedBits_Xpol ~= SIG.Xpol.bits));

    % Calculate the total number of bits
    totalBits_Xpol = numel(consolidatedBits_Xpol);

    % Calculate the Bit Error Rate (BER)
    BER_Xpol = bitErrors_Xpol / totalBits_Xpol;

    % Display the results
    fprintf('Total number of bits for X polarization: %d\n', totalBits_Xpol);
    fprintf('Number of bit errors for X polarization: %d\n', bitErrors_Xpol);
    fprintf('Bit Error Rate (BER) for X polariztion: %f\n', BER_Xpol);
end

if size(consolidatedBits_Ypol) ~= size(SIG.Ypol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Ypol = sum(sum(consolidatedBits_Ypol ~= SIG.Ypol.bits));

    % Calculate the total number of bits
    totalBits_Ypol = numel(consolidatedBits_Ypol);

    % Calculate the Bit Error Rate (BER)
    BER_Ypol = bitErrors_Ypol / totalBits_Ypol;

    % Display the results
    fprintf('Total number of bits for Y polarization: %d\n', totalBits_Ypol);
    fprintf('Number of bit errors for Y polarization: %d\n', bitErrors_Ypol);
    fprintf('Bit Error Rate (BER) for Y polariztion: %f\n', BER_Ypol);
end
