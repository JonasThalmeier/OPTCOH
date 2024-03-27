clear all;
close all;
clc;

% Load the .mat file
load('TXsequences/TXsequence_QPSK_64GBaud.mat');

% txSig is now loaded along with other structures like SIG and PulseShaping

% Apply matched filtering
% Assuming b_coeff is your filter coefficients from the PulseShaping structure

% Create the SRRC filter (same as used for TX or initial filtering)
srrcFilter = rcosdesign(PulseShaping.rho, PulseShaping.Span, SIG.Sps, 'sqrt');

% Assume rxSig is your received SRRC filtered signal
% Apply the matched filter (which is the same SRRC filter here)
%rxSig_Xpol = filter(srrcFilter, 1, SIG.Xpol.txSig);
%rxSig_Ypol = filter(srrcFilter, 1, SIG.Ypol.txSig);

delay = PulseShaping.Span * SIG.Sps / 2;
rxSig_Xpol = filter(PulseShaping.b_coeff, 1, SIG.Xpol.txSig);
%rxSig_Xpol = rxSig_Xpol(delay+1:end);
rxSig_Ypol = filter(PulseShaping.b_coeff, 1, SIG.Ypol.txSig);
%rxSig_Ypol = rxSig_Ypol(delay+1:end);

%-------------------------------------------------------------------------
% Select the first 30 elements
n1 = 1;
n2 = 200;
selectedElements = rxSig_Xpol(n1:n2);

% Calculate amplitude and phase
amplitude = abs(selectedElements);
phase = angle(selectedElements)*180/pi;

% Create a new figure
figure;

% Plot amplitude
subplot(2, 1, 1);  % This creates a subplot with 2 rows, 1 column, and selects the 1st subplot.
plot(n1:n2, amplitude, '-o');
title('Amplitude of the First 30 Elements of rxSig_Xpol');
xlabel('Element Index');
ylabel('Amplitude');

% Plot phase
subplot(2, 1, 2);  % This selects the 2nd subplot in the same figure.
plot(n1:n2, phase, '-x');
title('Phase of the First 30 Elements of rxSig_Xpol');
xlabel('Element Index');
ylabel('Phase (deg)');
yline(-180);
yline(-90);
yline(90);
yline(180);



%--------------------------------------------------------------------------

% Downsample the signal
% Assuming Sps is your number of samples per symbol
downsampledSig_Xpol = downsample(rxSig_Xpol, SIG.Sps);
downsampledSig_Ypol = downsample(rxSig_Ypol, SIG.Sps);

% Symbol Demapping
% The specifics of this process depend on your modulation scheme
% For QPSK as an example, you might have a demapping function like this:

% This is a placeholder; actual implementation will vary based on modulation and needs
demappedBits_Xpol = zeros(length(downsampledSig_Xpol),2); % Adjust size accordingly for bit pairs, etc.
demappedBits_Ypol = zeros(length(downsampledSig_Ypol),2); % Adjust size accordingly for bit pairs, etc.

% Assuming QPSK and not accounting for noise, just a direct mapping
for i = 1:length(downsampledSig_Xpol)
    % This is a simplistic approach; real demapping would consider noise, etc.
    if real(downsampledSig_Xpol(i)) > 0
        if imag(downsampledSig_Xpol(i)) > 0
            demappedBits_Xpol(i, :) = [0 0];
        else
            demappedBits_Xpol(i, :) = [0 1];
        end
    else
        if imag(downsampledSig_Xpol(i)) > 0
            demappedBits_Xpol(i, :) = [1 0];
        else
            demappedBits_Xpol(i, :) = [1 1];
        end
    end
end

for i = 1:length(downsampledSig_Ypol)
    % This is a simplistic approach; real demapping would consider noise, etc.
    if real(downsampledSig_Ypol(i)) > 0
        if imag(downsampledSig_Ypol(i)) > 0
            demappedBits_Ypol(i, :) = [0 0];
        else
            demappedBits_Ypol(i, :) = [0 1];
        end
    else
        if imag(downsampledSig_Ypol(i)) > 0
            demappedBits_Ypol(i, :) = [1 0];
        else
            demappedBits_Ypol(i, :) = [1 1];
        end
    end
end

% Now demappedBits contains your recovered bitstream.

% Let's assume demappedBits_Xpol is your bit outcomes with size [N * Npp, 2]
% Where N is the number of unique symbols and Npp is the number of repetitions

N = length(demappedBits_Xpol) / SIG.Npp;  % Calculate the number of unique symbols
consolidatedBits_Xpol = zeros(N, 2);  % Initialize the array for consolidated bits
for i = 1:N
    % Extract bits corresponding to the same symbol across all repetitions
    symbolBits = reshape(demappedBits_Xpol((i-1)*SIG.Npp + 1:i*SIG.Npp, :), [], 2, SIG.Npp);
    
    % For binary data, decide each bit based on the majority across repetitions
    % This example uses a simple majority, which is straightforward for binary outcomes
    for j = 1:2  % Assuming binary phase shift keying for simplification
        % Count the occurrences of '1's and '0's and decide based on majority
        bitMajority = mode(symbolBits(:, j, :), 3);
        consolidatedBits_Xpol(i, j) = bitMajority;
    end
end
consolidatedBits_Ypol = zeros(N, 2);  % Initialize the array for consolidated bits
for i = 1:N
    % Extract bits corresponding to the same symbol across all repetitions
    symbolBits = reshape(demappedBits_Ypol((i-1)*SIG.Npp + 1:i*SIG.Npp, :), [], 2, SIG.Npp);
    
    % For binary data, decide each bit based on the majority across repetitions
    % This example uses a simple majority, which is straightforward for binary outcomes
    for j = 1:2  % Assuming binary phase shift keying for simplification
        % Count the occurrences of '1's and '0's and decide based on majority
        bitMajority = mode(symbolBits(:, j, :), 3);
        consolidatedBits_Ypol(i, j) = bitMajority;
    end
end

% Now, consolidatedBits contains the 'averaged' bit decisions



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
