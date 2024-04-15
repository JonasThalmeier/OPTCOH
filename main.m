clear;
close all;
clc;

% Load the .mat file
% load('TXsequences/TXsequence_QPSK_64GBaud.mat');
load('TXsequences/TXsequence_16QAM_64GBaud.mat');

if size(SIG.Xpol.bits, 2) == 2
    MODULATIONS = 'QPSK';
else
    MODULATIONS = 'QAM16';
end

% Parameters

% txSig is now loaded along with other structures like SIG and PulseShaping

% Apply matched filtering
% Assuming b_coeff is your filter coefficients from the PulseShaping structure

% modulation = ["QPSK" "16-QAM"];
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
% fprintf('The transmitted moduluation is: %s\n', modulation(r));

%------------------Fiber Simulation----------------------------------------
% [Xpol_dc, Ypol_dc] = cd_sim(SIG.Xpol.txSig, SIG.Ypol.txSig, SIG.Sps, SIG.symbolRate, 17, .001, 1550);
Xpol = SIG.Xpol.txSig;
Ypol = SIG.Ypol.txSig;

r = randi(8,1);
Xpol = Xpol(r:end);
Ypol = Ypol(r:end);

rxSig_Xpol = conv(PulseShaping.b_coeff, Xpol);
rxSig_Ypol = conv(PulseShaping.b_coeff, Ypol);

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

downsampledSig_Xpol = downsample(rxSig_Xpol, SIG.Sps/2);
downsampledSig_Ypol = downsample(rxSig_Ypol, SIG.Sps/2);

[downsampledSig_Xpol,downsampledSig_Ypol] = samp_phase_recovery(downsampledSig_Xpol,downsampledSig_Ypol,SIG.Sps);

[c,lags] = xcorr(downsampledSig_Xpol(2:2:length(SIG.Xpol.txSymb)), SIG.Xpol.txSymb);
stem(lags,real(c));

[M,I] = max(c);

downsampledSig_Xpol = downsampledSig_Xpol(2*(lags(I)+1):end-2*(lags(I)), :);
downsampledSig_Ypol = downsampledSig_Ypol(2*(lags(I)+1):end-2*(lags(I)), :);

downsampledSig_Xpol = downsampledSig_Xpol/abs(real(median(downsampledSig_Xpol(2:2:end)))); % normalize over the median value since gaussian shape, take oly real part because it represents the unit in the non-normalized case
downsampledSig_Ypol = downsampledSig_Ypol/abs(real(median(downsampledSig_Ypol(2:2:end))));

% Plot constellation
figure;
scatter(real(downsampledSig_Xpol(2:2:end)), imag(downsampledSig_Xpol(2:2:end)), ".", "k");
grid on;

% max_energy = max(abs(downsampledSig_Xpol(2:2:end)));
% 
% if max_energy < 2
%     fprintf('The tracked moduluation is: QPSK\n');
%     M = 2;
%     [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QPSK_demapping(downsampledSig_Xpol,downsampledSig_Ypol);
% else
%     fprintf('The tracked moduluation is: 16-QAM\n');
%     M = 4;
%     [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QAM_16_demapping(downsampledSig_Xpol,downsampledSig_Ypol);
% end

[counts, binEdges] = histcounts(angle(downsampledSig_Xpol(2:2:end)), 12, 'Normalization', 'probability');
%histogram(angle(downsampledSig_Xpol(2:2:end)), 12, 'Normalization', 'probability');
if max(counts) > .17
    fprintf('The tracked moduluation is: QPSK\n');
    M = 2;
    [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QPSK_demapping(downsampledSig_Xpol,downsampledSig_Ypol);
else
    fprintf('The tracked moduluation is: 16-QAM\n');
    M = 4;
    [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QAM_16_demapping(downsampledSig_Xpol,downsampledSig_Ypol);
end


%----------------Majority voting over the repeated symbols-----------------
% Let's assume demappedBits_Xpol is your bit outcomes with size [N * Npp, 2]
% Where N is the number of unique symbols and Npp is the number of repetitions
N = length(SIG.Xpol.txSymb);  % Calculate the number of unique symbols
[consolidatedBits_Xpol,consolidatedBits_Ypol] = voting(N, M, demappedBits_Xpol, demappedBits_Ypol);



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