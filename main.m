clear;
close all;
clc;

% Load the .mat file
load('TXsequences/TXsequence_16QAM_64GBaud.mat');

% Variables
modulation = 'QAM16';
SpS = 8;

% txSig is now loaded along with other structures like SIG and PulseShaping

% Apply matched filtering
% Assuming b_coeff is your filter coefficients from the PulseShaping structure

% modulationArray = {'QPSK', 'QAM16'};
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
% fprintf('The chosen modluation is: ');
% modulation = modulationArray(r)  % Extract the modulation for this random number.


rxSig_Xpol = conv(PulseShaping.b_coeff, SIG.Xpol.txSig);
rxSig_Ypol = conv(PulseShaping.b_coeff, SIG.Ypol.txSig);

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

downsampledSig_Xpol = own_downsample(rxSig_Xpol, SIG.Sps, SIG.Npp);
downsampledSig_Ypol = own_downsample(rxSig_Ypol, SIG.Sps, SIG.Npp);

[c,lags] = xcorr(downsampledSig_Xpol(1:length(SIG.Xpol.txSymb)),SIG.Xpol.txSymb);
stem(lags,c)

[M,I] = max(c);

downsampledSig_Xpol = downsampledSig_Xpol(lags(I)+1:end-lags(I), :);
downsampledSig_Ypol = downsampledSig_Ypol(lags(I)+1:end-lags(I), :);

% Plot constellation
figure;
scatter(real(downsampledSig_Xpol), imag(downsampledSig_Xpol), ".", "k");
grid on;

%[demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QPSK_demapping(downsampledSig_Xpol,downsampledSig_Ypol);
[demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QAM_16_demapping(downsampledSig_Xpol,downsampledSig_Ypol);
%%

%----------------Majority voting over the repeated symbols-----------------
% Let's assume demappedBits_Xpol is your bit outcomes with size [N * Npp, 2]
% Where N is the number of unique symbols and Npp is the number of repetitions

% N = length(demappedBits_Xpol) / SIG.Npp;  % Calculate the number of unique symbols
N = length(SIG.Xpol.txSymb);  % Calculate the number of unique symbols

consolidatedBits_Xpol = zeros(N, 2);
consolidatedBits_Ypol = zeros(N, 2);
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
