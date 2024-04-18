clear;
close all;
clc;

% Load the .mat file
QPSK = load('TXsequences/TXsequence_QPSK_64GBaud.mat');
QAM16 = load('TXsequences/TXsequence_16QAM_64GBaud.mat');

MODULATIONS = [QPSK,QAM16];

modulation = ["QPSK" "16-QAM"];
r = randi([1, 2], 1); % Get a 1 or 2 randomly.
fprintf('The transmitted moduluation is: %s\n', modulation(r));

% Parameters

SpS_down = 4;
SpS_up = MODULATIONS(r).SIG.Sps/SpS_down;

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

rxSig_Xpol = conv(MODULATIONS(r).PulseShaping.b_coeff, MODULATIONS(r).SIG.Xpol.txSig);
rxSig_Ypol = conv(MODULATIONS(r).PulseShaping.b_coeff, MODULATIONS(r).SIG.Ypol.txSig);

% Create delay and phase convolved signals
[delay_phase_distorted_RX_Xpol, delay_phase_distorted_RX_Ypol] = DP_Distortion(MODULATIONS(r).SIG.Xpol.txSig, MODULATIONS(r).SIG.Ypol.txSig, MODULATIONS(r).PulseShaping.b_coeff);

% Add of the noise
[delay_phase_noise_distorted_RX_Xpol, Noise] = WGN_Noise_Generation(delay_phase_distorted_RX_Xpol,MODULATIONS(r).SIG.Sps, M);
[delay_phase_noise_distorted_RX_Ypol, Noise2] = WGN_Noise_Generation(delay_phase_distorted_RX_Ypol,MODULATIONS(r).SIG.Sps, M);

%%
% ----- Recover from delay and phase at 8 SpS, this way no problems with downsampling

recovered_TXSig_Xpol = Recover_Delay_Phase_Noise(rxSig_Xpol,delay_phase_noise_distorted_RX_Xpol);
recovered_TXSig_Ypol = Recover_Delay_Phase_Noise(rxSig_Ypol,delay_phase_noise_distorted_RX_Ypol);



% %-------Plotting Phase and Amplitude of the filtered signal----------------
% % Select the first 30 elements
% n1 = 1;
% n2 = 200;
% selectedElements = rxSig_Xpol(n1:n2);
% 
% % Calculate amplitude and phase
% amplitude = abs(selectedElements);
% phase = angle(selectedElements)*180/pi;
% 
% % Create a new figure
% figure;
% 
% % Plot amplitude
% subplot(2, 1, 1);  % This creates a subplot with 2 rows, 1 column, and selects the 1st subplot.
% plot(n1:n2, amplitude, '-o');
% title('Amplitude of the First 30 Elements of rxSig_Xpol');
% xlabel('Element Index');
% ylabel('Amplitude');
% 
% % Plot phase
% subplot(2, 1, 2);  % This selects the 2nd subplot in the same figure.
% plot(n1:n2, phase, '-x');
% title('Phase of the First 30 Elements of rxSig_Xpol');
% xlabel('Element Index');
% ylabel('Phase (deg)');
% yline(-180);
% yline(-90);
% yline(90);
% yline(180);

%----------------Downsample/Demapping the signal---------------------------

downsampledSig_Xpol = downsample(recovered_TXSig_Xpol, SpS_down);
downsampledSig_Ypol = downsample(recovered_TXSig_Ypol, SpS_down);

downsampledSig_Xpol = downsampledSig_Xpol/abs(real(median(downsampledSig_Xpol))); % normalize over the median value since gaussian shape, take oly real part because it represents the unit in the non-normalized case
downsampledSig_Ypol = downsampledSig_Ypol/abs(real(median(downsampledSig_Ypol)));


% Apply 2 SpS decoding

upsampledSig_Xpol_txSymb = upsample(MODULATIONS(r).SIG.Xpol.txSymb, SpS_up);
upsampledSig_Xpol_txSymb(2:2:end) = upsampledSig_Xpol_txSymb(1:2:end);

upsampledSig_Ypol_txSymb = upsample(MODULATIONS(r).SIG.Ypol.txSymb, SpS_up);
upsampledSig_Ypol_txSymb(2:2:end) = upsampledSig_Ypol_txSymb(1:2:end);

% Correlation sample 1
downsampledSig_Xpol_1 = downsampledSig_Xpol(1:2:length(downsampledSig_Xpol));
downsampledSig_Ypol_1 = downsampledSig_Ypol(1:2:length(downsampledSig_Ypol));

figure();
[c1,lags1] = xcorr(downsampledSig_Xpol_1(1:length(upsampledSig_Xpol_txSymb(1:2:end))), upsampledSig_Xpol_txSymb(1:2:end));
stem(lags1,real(c1))
[M1,I1] = max(c1);
downsampledSig_Xpol_1 = downsampledSig_Xpol_1(lags1(I1)+1:end-lags1(I1), :);

[c1,lags1] = xcorr(downsampledSig_Ypol_1(1:length(upsampledSig_Ypol_txSymb(1:2:end))), upsampledSig_Ypol_txSymb(1:2:end));
[M1,I1] = max(c1);
downsampledSig_Ypol_1 = downsampledSig_Ypol_1(lags1(I1)+1:end-lags1(I1), :);

figure();
[c1,lags1] = xcorr( downsampledSig_Xpol_1(1:length(upsampledSig_Xpol_txSymb(1:2:end-1))), upsampledSig_Xpol_txSymb(1:2:end-1));
stem(lags1,real(c1))

downsampledSig_Xpol_1 = downsampledSig_Xpol_1/abs(real(median(downsampledSig_Xpol_1))); % normalize over the median value since gaussian shape, take oly real part because it represents the unit in the non-normalized case
downsampledSig_Ypol_1 = downsampledSig_Ypol_1/abs(real(median(downsampledSig_Ypol_1)));

% Plot constellation
figure;
scatter(real(downsampledSig_Xpol_1), imag(downsampledSig_Xpol_1), ".", "k");
grid on;



% Correlation sample 2
downsampledSig_Xpol_2 = downsampledSig_Xpol(2:2:length(downsampledSig_Xpol));
downsampledSig_Ypol_2 = downsampledSig_Ypol(2:2:length(downsampledSig_Ypol));

figure();
[c2,lags2] = xcorr(downsampledSig_Xpol_2(1:length(upsampledSig_Xpol_txSymb(2:2:end))), upsampledSig_Xpol_txSymb(2:2:end));
stem(lags2,real(c2))
[M2,I2] = max(c2);
downsampledSig_Xpol_2 = downsampledSig_Xpol_2(lags2(I2)+1:end-lags2(I2), :);

[c2,lags2] = xcorr(downsampledSig_Ypol_2(1:length(upsampledSig_Ypol_txSymb(2:2:end))), upsampledSig_Ypol_txSymb(2:2:end));
[M2,I2] = max(c2);
downsampledSig_Ypol_2 = downsampledSig_Ypol_2(lags2(I2)+1:end-lags2(I2), :);

figure();
[c2,lags2] = xcorr( downsampledSig_Xpol_2(1:length(upsampledSig_Xpol_txSymb(2:2:end))), upsampledSig_Xpol_txSymb(2:2:end));
stem(lags2,real(c2))

downsampledSig_Xpol_2 = downsampledSig_Xpol_2/abs(real(median(downsampledSig_Xpol_2))); % normalize over the median value since gaussian shape, take oly real part because it represents the unit in the non-normalized case
downsampledSig_Ypol_2 = downsampledSig_Ypol_2/abs(real(median(downsampledSig_Ypol_2)));

% Plot constellation
figure;
scatter(real(downsampledSig_Xpol_2), imag(downsampledSig_Xpol_2), ".", "k");
grid on;


% Decide which phase is better using the variance of the energies

downsampledSig_Xpol = downsampledSig_Xpol_1;
downsampledSig_Ypol = downsampledSig_Ypol_1;
% if mean(histcounts(abs(downsampledSig_Xpol_1))) > mean(histcounts(abs(downsampledSig_Xpol_2)))
%     downsampledSig_Xpol = downsampledSig_Xpol_1;
% else
%     downsampledSig_Xpol = downsampledSig_Xpol_2;
% end
% 
% if mean(histcounts(abs(downsampledSig_Ypol_1))) > mean(histcounts(abs(downsampledSig_Ypol_2)))
%     downsampledSig_Ypol = downsampledSig_Ypol_1;
% else
%     downsampledSig_Ypol = downsampledSig_Ypol_2;
% end

%clear downsampledSig_Xpol_1 downsampledSig_Xpol_2 downsampledSig_Ypol_2 downsampledSig_Ypol_1;

% Plot constellation
figure;
scatter(real(downsampledSig_Xpol), imag(downsampledSig_Xpol), ".", "k");
title('Xpol constellation');
grid on;

figure;
scatter(real(downsampledSig_Ypol), imag(downsampledSig_Ypol), ".", "k");
title('Ypol constellation');
grid on;
%%
mean_energy = mean(abs(downsampledSig_Xpol));

if mean_energy < 2
    fprintf('The tracked moduluation is: QPSK\n');
    [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QPSK_demapping(downsampledSig_Xpol,downsampledSig_Ypol);
else
    fprintf('The tracked moduluation is: 16-QAM\n');
    [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QAM_16_demapping(downsampledSig_Xpol,downsampledSig_Ypol);
end

%----------------Majority voting over the repeated symbols-----------------
% Let's assume demappedBits_Xpol is your bit outcomes with size [N * Npp, 2]
% Where N is the number of unique symbols and Npp is the number of repetitions

% N = length(demappedBits_Xpol) / SIG.Npp;  % Calculate the number of unique symbols
N = length(MODULATIONS(r).SIG.Xpol.txSymb);  % Calculate the number of unique symbols

consolidatedBits_Xpol = zeros(N, size(MODULATIONS(r).SIG.Xpol.bits,2));
consolidatedBits_Ypol = zeros(N, size(MODULATIONS(r).SIG.Xpol.bits,2));

for i = 1:N
    consolidatedBits_Xpol(i, :) = mode(demappedBits_Xpol(i:N:end,:),1);
    consolidatedBits_Ypol(i, :) = mode(demappedBits_Ypol(i:N:end,:),1);
end
% Now, consolidatedBits contains the 'averaged' bit decisions



%-------------------------BER calculation----------------------------------

if size(consolidatedBits_Xpol) ~= size(MODULATIONS(r).SIG.Xpol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Xpol = sum(sum(consolidatedBits_Xpol ~= MODULATIONS(r).SIG.Xpol.bits));

    % Calculate the total number of bits
    totalBits_Xpol = numel(consolidatedBits_Xpol);

    % Calculate the Bit Error Rate (BER)
    BER_Xpol = bitErrors_Xpol / totalBits_Xpol;

    % Display the results
    fprintf('Total number of bits for X polarization: %d\n', totalBits_Xpol);
    fprintf('Number of bit errors for X polarization: %d\n', bitErrors_Xpol);
    fprintf('Bit Error Rate (BER) for X polariztion: %f\n', BER_Xpol);
end

if size(consolidatedBits_Ypol) ~= size(MODULATIONS(r).SIG.Ypol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Ypol = sum(sum(consolidatedBits_Ypol ~= MODULATIONS(r).SIG.Ypol.bits));

    % Calculate the total number of bits
    totalBits_Ypol = numel(consolidatedBits_Ypol);

    % Calculate the Bit Error Rate (BER)
    BER_Ypol = bitErrors_Ypol / totalBits_Ypol;

    % Display the results
    fprintf('Total number of bits for Y polarization: %d\n', totalBits_Ypol);
    fprintf('Number of bit errors for Y polarization: %d\n', bitErrors_Ypol);
    fprintf('Bit Error Rate (BER) for Y polariztion: %f\n', BER_Ypol);
end
