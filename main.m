clear all;
close all;
clc;

% Load the .mat file
load('TXsequences/TXsequence_QPSK_64GBaud.mat');
% txSig is now loaded along with other structures like SIG and PulseShaping

% Apply matched filtering
% Assuming b_coeff is your filter coefficients from the PulseShaping structure
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

downsampledSig_Xpol = downsample(rxSig_Xpol, SIG.Sps);
downsampledSig_Ypol = downsample(rxSig_Ypol, SIG.Sps);
% downsampledSig_Xpol = downsampledSig_Xpol(16:end, :);
% downsampledSig_Ypol = downsampledSig_Ypol(16:end, :);
downsampledSig_Xpol = downsampledSig_Xpol(16:end-15, :);
downsampledSig_Ypol = downsampledSig_Ypol(16:end-15, :);

% Plot constalation
figure;
scatter(real(downsampledSig_Xpol), imag(downsampledSig_Xpol), ".", "k");

% Symbol Demapping
demappedBits_Xpol = zeros(length(downsampledSig_Xpol),2); % Adjust size accordingly for bit pairs, etc.
demappedSymb_Xpol = zeros(length(downsampledSig_Xpol),1);
demappedBits_Ypol = zeros(length(downsampledSig_Ypol),2); % Adjust size accordingly for bit pairs, etc.
demappedSymb_Ypol = zeros(length(downsampledSig_Ypol),1);

% Assuming QPSK and not accounting for noise, just a direct mapping
for i = 1:length(downsampledSig_Xpol)
    % This is a simplistic approach; real demapping would consider noise, etc.
    if real(downsampledSig_Xpol(i)) > 0
        if imag(downsampledSig_Xpol(i)) > 0
            demappedBits_Xpol(i, :) = [1 0];
            demappedSymb_Xpol(i) = 2;
        else
            demappedBits_Xpol(i, :) = [1 1];
            demappedSymb_Xpol(i) = 3;
        end
    else
        if imag(downsampledSig_Xpol(i)) > 0
            demappedBits_Xpol(i, :) = [0 0];
            demappedSymb_Xpol(i) = 0;
        else
            demappedBits_Xpol(i, :) = [0 1];
            demappedSymb_Xpol(i) = 1;
        end
    end
end

for i = 1:length(downsampledSig_Ypol)
    % This is a simplistic approach; real demapping would consider noise, etc.
    if real(downsampledSig_Ypol(i)) > 0
        if imag(downsampledSig_Ypol(i)) > 0
            demappedBits_Ypol(i, :) = [1 0];
            demappedSymb_Ypol(i) = 2;
        else
            demappedBits_Ypol(i, :) = [1 1];
            demappedSymb_Ypol(i) = 3;
        end
    else
        if imag(downsampledSig_Ypol(i)) > 0
            demappedBits_Ypol(i, :) = [0 0];
            demappedSymb_Ypol(i) = 0;
        else
            demappedBits_Ypol(i, :) = [0 1];
            demappedSymb_Ypol(i) = 1;
        end
    end
end



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
