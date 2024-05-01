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
if r == 1
    M = 4;
else
    M = 16;
end

% Create delay and phase convolved signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig);

% Adding the noise
[X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_distorted,SIG.Sps, M, 20);
[Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_distorted,SIG.Sps, M, 20);

%----------------Pulse shaping and Downsample the signal-------------------
X_distorted_AWGN = conv(PulseShaping.b_coeff, X_distorted_AWGN);
Y_distorted_AWGN = conv(PulseShaping.b_coeff, Y_distorted_AWGN);

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
X_dist_2SpS = X_distorted_AWGN(1:4:end);
Y_dist_2SpS = Y_distorted_AWGN(1:4:end);

% Define range for number of taps
% numtaps_values = linspace(5, 33, 15); % 10 values from 5 to 150

% Preallocate array for error storage
% errors = zeros(1, length(numtaps_values));

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

EQ = comm.LinearEqualizer('Algorithm', 'CMA', 'StepSize', stepsize,'NumTaps', numtaps, 'InputSamplesPerSymbol', 2, 'Constellation', constellation, 'ReferenceTap', referencetap, 'InputDelay', 0);
[X_eq,errX] = EQ(X_dist_2SpS(1:end-mod(length(X_dist_2SpS),2)));
[Y_eq,errY] = EQ(Y_dist_2SpS(1:end-mod(length(Y_dist_2SpS),2)));
% figure; plot(abs(err));
% constell = comm.ConstellationDiagram('NumInputPorts', 1, 'SamplesPerSymbol', SpS_up, 'ReferenceConstellation', constellation, 'Title', 'Before phase correction');
% constell(X_eq);
% scatterplot(X_eq);

carrSynch = comm.CarrierSynchronizer("Modulation", modulation,"SamplesPerSymbol", 1);
[X_eq, phEstX] = carrSynch(X_eq);
[Y_eq, phEstY] = carrSynch(Y_eq);

%-----------------Consolidation--------------------------------------------
N = length(SIG.Xpol.txSymb);  % Calculate the number of unique symbols
X_consolidated = zeros(N, 1);
Y_consolidated = zeros(N, 1);
for idx = 1:N
    X_consolidated(idx) = mean(X_eq(idx:N:end,:));
    Y_consolidated(idx) = mean(Y_eq(idx:N:end,:));
end

%-----------find right rotation and delay----------------------------------
X_avg_amp = mean(abs(X_consolidated))/sqrt(2);
X_rotation = zeros(4,1);
X_delay = zeros(4,1);
Y_avg_amp = mean(abs(X_consolidated))/sqrt(2);
Y_rotation = zeros(4,1);
Y_delay = zeros(4,1);
for idx=1:4
    Xrotated = X_consolidated.*exp(1i*idx*pi/2);
    Xdelay(idx) = finddelay(Xrotated, SIG.Xpol.txSymb);
    Xrotated = circshift(Xrotated,Xdelay(idx));
    X_rotation(idx)=std(Xrotated-X_avg_amp.*SIG.Xpol.txSymb);
    Yrotated = Y_consolidated.*exp(1i*idx*pi/2);
    Ydelay(idx) = finddelay(Yrotated, SIG.Ypol.txSymb);
    Yrotated = circshift(Yrotated,Ydelay(idx));
    Y_rotation(idx)=std(Yrotated-Y_avg_amp.*SIG.Ypol.txSymb);
end
[M I] = min(X_rotation);
X_consolidated_recovered = circshift(X_consolidated.*exp(1i*I*pi/2), Xdelay(I));
[M I] = min(X_rotation);
Y_consolidated_recovered = circshift(Y_consolidated.*exp(1i*I*pi/2), Ydelay(I));
scatterplot(X_consolidated_recovered)

if r == 1
    M = 2;
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_consolidated_recovered,Y_consolidated_recovered);
else
    M = 4;
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_consolidated_recovered,Y_consolidated_recovered);
end

%-------------------------BER calculation----------------------------------

if size(X_demappedBits) ~= size(SIG.Xpol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Xpol = sum(sum(X_demappedBits ~= SIG.Xpol.bits));

    % Calculate the total number of bits
    totalBits_Xpol = numel(X_demappedBits);

    % Calculate the Bit Error Rate (BER)
    BER_Xpol = bitErrors_Xpol / totalBits_Xpol;

    % Display the results
    fprintf('Total number of bits for X polarization: %d\n', totalBits_Xpol);
    fprintf('Number of bit errors for X polarization: %d\n', bitErrors_Xpol);
    fprintf('Bit Error Rate (BER) for X polariztion: %e\n', BER_Xpol);
end

if size(Y_demappedBits) ~= size(SIG.Ypol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Ypol = sum(sum(Y_demappedBits ~= SIG.Ypol.bits));

    % Calculate the total number of bits
    totalBits_Ypol = numel(Y_demappedBits);

    % Calculate the Bit Error Rate (BER)
    BER_Ypol = bitErrors_Ypol / totalBits_Ypol;

    % Display the results
    fprintf('Total number of bits for Y polarization: %d\n', totalBits_Ypol);
    fprintf('Number of bit errors for Y polarization: %d\n', bitErrors_Ypol);
    fprintf('Bit Error Rate (BER) for Y polariztion: %e\n', BER_Ypol);
end
