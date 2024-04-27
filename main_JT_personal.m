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
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig);

% Adding the noise
[X_distorted_AGWN, NoiseX] = WGN_Noise_Generation(X_distorted,SIG.Sps, M, 20);
%[Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_distorted,SIG.Sps, M, 15);

%----------------Pulse shaping and Downsample the signal-------------------
X_distorted_AGWN = conv(PulseShaping.b_coeff, X_distorted_AGWN);
%Y_distorted_AWGN = conv(PulseShaping.b_coeff, Y_distorted_AWGN);
% X_2Sps = downsample(X_distorted_AGWN, SpS_down);
% Y_2Sps = downsample(Y_distorted_AWGN, SpS_down);
% 
% The code below is done by the built-in function FINDDELAY
% [corr1, lag_1] = xcorr(rxSig_Xpol(1:65536), X_distorted_AGWN(1:65536));
% figure(), stem(lag_1,abs(corr1));
% [max_corr, max_index] = max(abs(corr1));
% fprintf('The Xpol tracked delay is of %d samples.\n', abs(lag_1(max_index)));

% X_distorted_AGWN = downsample(X_distorted_AGWN,SpS_down);
% rxSig_Xpol = downsample(rxSig_Xpol,SpS_down);
rx_Xpol_Delay = finddelay(X_distorted_AGWN, rxSig_Xpol);
fprintf('The Xpol tracked delay is of %d samples.\n', abs(rx_Xpol_Delay));

% X_distorted_AGWN = X_distorted_AGWN(abs(lag_1(max_index))+1:end);
% [corr1, lag_1] = xcorr(rxSig_Xpol(1:65536),X_distorted_AGWN);
% figure(), stem(lag_1,abs(corr1));

%---------------------------EQ---------------------------------------------
%From 15 numtaps are the best for qpsk, also for qam, but 60 and 120 seem
%better. The reference tap for qpsk are 1,3,65,121.
if r == 1
    constellation = pskmod(0:3, 4, pi/4);
    stepsize = 1e-3; %best result
    numtaps = 9;
    referencetap = (numtaps-1)/2;
    modulation = 'QPSK';
else
    constellation = qammod(0:15, 16);
    stepsize = 5e-6;
    numtaps = 9; % bigger filter works better here
    referencetap = (numtaps-1)/2;
    modulation = 'QAM';
end

aw = true;
X_2_dist = X_distorted_AGWN(1:4:end);

% Define range for number of taps
numtaps_values = linspace(5, 33, 15); % 10 values from 5 to 150

% Preallocate array for error storage
errors = zeros(1, length(numtaps_values));

% % Create figure for constellation plots
% figure;
% % Loop through each tap value
% for idx = 1:length(numtaps_values)
%     numtaps = numtaps_values(idx);
%     EQ = comm.LinearEqualizer('Algorithm', 'CMA', 'StepSize', stepsize, ...
%                               'NumTaps', numtaps, 'InputSamplesPerSymbol', 2, ...
%                               'Constellation', constellation, ...
%                               'ReferenceTap', floor((numtaps-1)/2), 'InputDelay', abs(rx_Xpol_Delay)/4);
% 
%     % Equalize the signal
%     [X_eq, err] = EQ(X_2_dist(1:end-mod(length(X_2_dist),2)));
%     errors(idx) = norm(err);
% 
%     % Plot constellation diagram for current number of taps
%     subplot(3, 5, idx);
%     plot(real(X_eq), imag(X_eq), '.');
%     title(['Constellation with ', num2str(numtaps), ' Taps']);
%     xlabel('In-phase');
%     ylabel('Quadrature');
%     axis square; % Keeping the aspect ratio square
%     grid on;
% end
% 
% % Create figure for error plot
% figure;
% plot(numtaps_values, errors, '-o');
% title('Error vs. Number of Taps');
% xlabel('Number of Taps');
% ylabel('Error (norm)');
% legend('Error', 'Location', 'best');
% grid on;


% % Define ranges for number of taps and step sizes
% numtaps_values = linspace(5, 33, 15); % 10 values from 5 to 150
% stepsize_values = logspace(-8.5, -7, 15); % Example step sizes
% 
% % Initialize matrix to store errors for heatmap
% error_matrix = zeros(length(numtaps_values), length(stepsize_values));
% 
% % Initialize variables to track the best configuration
% min_error = inf; % Start with a very large number
% best_numtaps = 0;
% best_stepsize = 0;
% 
% % Grid search over the number of taps and step sizes
% for i = 1:length(numtaps_values)
%     for j = 1:length(stepsize_values)
%         numtaps = numtaps_values(i);
%         stepsize = stepsize_values(j);
%         EQ = comm.LinearEqualizer('Algorithm', 'CMA', 'StepSize', stepsize, ...
%                                   'NumTaps', numtaps, 'InputSamplesPerSymbol', 2, ...
%                                   'Constellation', constellation, ...
%                                   'ReferenceTap', (numtaps-1)/2, 'InputDelay', abs(rx_Xpol_Delay)/4);
% 
%         % Equalize the signal
%         [~, err] = EQ(X_2_dist(1:end-mod(length(X_2_dist),2)));
%         current_error = norm(err); % Compute the norm of the error vector
% 
%         % Store error in matrix
%         error_matrix(i, j) = current_error;
% 
%         % Check if the current configuration is better
%         if current_error < min_error
%             min_error = current_error;
%             best_numtaps = numtaps;
%             best_stepsize = stepsize;
%         end
%     end
% end
% 
% % Display the best configuration and the minimum error
% fprintf('Best Configuration: NumTaps = %d, StepSize = %f\n', best_numtaps, best_stepsize);
% fprintf('Minimum Error: %f\n', min_error);
% 
% % Plotting the heatmap of errors
% figure;
% heatmap(stepsize_values, numtaps_values, error_matrix);
% xlabel('Step Size');
% ylabel('Number of Taps');
% title('Heatmap of Errors for Different Configurations');
% colorbar; % Add a colorbar to indicate error magnitude

EQ = comm.LinearEqualizer('Algorithm', 'CMA', 'StepSize', stepsize,'NumTaps', numtaps, 'InputSamplesPerSymbol', SIG.Sps, 'Constellation', constellation, 'ReferenceTap', referencetap, 'InputDelay', abs(rx_Xpol_Delay));
[X_eq,err] = EQ(X_distorted_AGWN(1:end-mod(length(X_distorted_AGWN),8)));
figure; plot(abs(err));
% constell = comm.ConstellationDiagram('NumInputPorts', 1, 'SamplesPerSymbol', SpS_up, 'ReferenceConstellation', constellation, 'Title', 'Before phase correction');
% constell(X_eq);
scatterplot(X_eq);
%%

carrSynch = comm.CarrierSynchronizer("Modulation", modulation,"SamplesPerSymbol", SIG.Sps);

[X_eq, phEst] = carrSynch(X_eq(1:end-mod(length(X_eq),8)));
fprintf('The random phase recovered is (degrees): %d\n', (mean(phEst) *180 /pi));
% constell2 = comm.ConstellationDiagram('NumInputPorts', 1, 'SamplesPerSymbol', SIG.Sps, 'ReferenceConstellation', constellation, 'Title', 'After phase correction');
% constell2(X_eq);

transient_Xpol = abs(finddelay(X_eq(1:65536), SIG.Xpol.txSymb));
X_eq = X_eq(transient_Xpol+1:end-transient_Xpol);
X_eq = X_eq(1:end-mod(length(X_eq),8));
% constell3 = comm.ConstellationDiagram('NumInputPorts', 1, 'SamplesPerSymbol', SIG.Sps, 'ReferenceConstellation', constellation, 'Title', 'After phase correction');
% constell3(X_eq);

X_2Sps = downsample(X_eq, SpS_down);
% constell4 = comm.ConstellationDiagram('NumInputPorts', 1, 'SamplesPerSymbol', SpS_up, 'ReferenceConstellation', constellation, 'Title', 'After phase correction');
% constell4(X_2Sps);

Y_2Sps =X_2Sps; %Just to give it a value, for the moment i test only the X_pol

[counts, binEdges] = histcounts(angle(X_2Sps(1:2:end)), 12, 'Normalization', 'probability');
if max(counts) > .17
    fprintf('The tracked moduluation is: QPSK\n');
    M = 4;
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_2Sps,Y_2Sps);
%     X_demappedBits_2 = pskdemod(X_2Sps(1:2:end),M, pi/4); it doesn't demodulate in the same way as our function
else
    fprintf('The tracked moduluation is: 16-QAM\n');
    M = 16;
%     X_demappedBits = qamdemod(X_2Sps(1:2:end),M);
 [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_2Sps,Y_2Sps);
end


X_BER = sum(~isequal(X_demappedSymb(1:length(SIG.Xpol.decSymbols)), SIG.Xpol.decSymbols))/length(SIG.Xpol.txSymb);
fprintf('The SER on Xpol is: %.9f\n', X_BER);
% scatterplot(X_eq);
% eyediagram(X_eq,2*SpS_down);

% plot(abs(err))
% xlabel('Symbols')
% ylabel('Error Magnitude')
% grid on
% title('Time-Varying Channel Without Retraining')
%%

%-------------Remove Transient at the end of transmission-----------------
X_2Sps = X_2Sps(1:end-(length(PulseShaping.b_coeff)/(SIG.Sps)));
Y_2Sps = Y_2Sps(1:end-(length(PulseShaping.b_coeff)/(SIG.Sps)));


% Plot constellation
figure;
scatter(real(X_2Sps(1:2:end)), imag(X_2Sps(1:2:end)), ".", "k");
title('Xpol constellation');
grid on;

figure;
scatter(real(Y_2Sps(1:2:end)), imag(Y_2Sps(1:2:end)), ".", "k");
title('Ypol constellation');
grid on;

[counts, binEdges] = histcounts(angle(X_2Sps(1:2:end)), 12, 'Normalization', 'probability');
if max(counts) > .17
    fprintf('The tracked moduluation is: QPSK\n');
    M = 2;
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_2Sps,Y_2Sps);
else
    fprintf('The tracked moduluation is: 16-QAM\n');
    M = 4;
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_2Sps,Y_2Sps);
end

%----------------Majority voting over the repeated symbols-----------------
% Let's assume demappedBits_Xpol is your bit outcomes with size [N * Npp, 2]
% Where N is the number of unique symbols and Npp is the number of repetitions

N = length(SIG.Xpol.txSymb);  % Calculate the number of unique symbols

X_consolidated = zeros(N, size(SIG.Xpol.bits,2));
Y_consolidatedBits = zeros(N, size(SIG.Xpol.bits,2));

for i = 1:N
    X_consolidated(i, :) = mode(X_demappedBits(i:N:end,:),1);
    Y_consolidatedBits(i, :) = mode(Y_demappedBits(i:N:end,:),1);
end
% Now, consolidatedBits contains the 'averaged' bit decisions




% % ----- Recover from delay and phase -------------------------------------
% X_2Sps = Recover_Delay_Phase_Noise(X_2Sps,SIG.Xpol.txSymb);
% Y_2Sps = Recover_Delay_Phase_Noise(Y_2Sps,SIG.Ypol.txSymb);
%-------------------------BER calculation----------------------------------

if size(X_consolidated) ~= size(SIG.Xpol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Xpol = sum(sum(X_consolidated ~= SIG.Xpol.bits));

    % Calculate the total number of bits
    totalBits_Xpol = numel(X_consolidated);

    % Calculate the Bit Error Rate (BER)
    BER_Xpol = bitErrors_Xpol / totalBits_Xpol;

    % Display the results
    fprintf('Total number of bits for X polarization: %d\n', totalBits_Xpol);
    fprintf('Number of bit errors for X polarization: %d\n', bitErrors_Xpol);
    fprintf('Bit Error Rate (BER) for X polariztion: %f\n', BER_Xpol);
end

if size(Y_consolidatedBits) ~= size(SIG.Ypol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Ypol = sum(sum(Y_consolidatedBits ~= SIG.Ypol.bits));

    % Calculate the total number of bits
    totalBits_Ypol = numel(Y_consolidatedBits);

    % Calculate the Bit Error Rate (BER)
    BER_Ypol = bitErrors_Ypol / totalBits_Ypol;

    % Display the results
    fprintf('Total number of bits for Y polarization: %d\n', totalBits_Ypol);
    fprintf('Number of bit errors for Y polarization: %d\n', bitErrors_Ypol);
    fprintf('Bit Error Rate (BER) for Y polariztion: %f\n', BER_Ypol);
end

function y = phaseCorrection(y)
a = angle(y((real(y) > 0) & (imag(y) > 0)));
a(a < 0.1) = a(a < 0.1) + pi/2;
theta = mean(a) - pi/4;
y = y * exp(-1i*theta);
fprintf('The random phase recovered is (degrees): %d\n', (theta *180 /pi));
end
