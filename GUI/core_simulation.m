function [Ber_Tot] = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac, scatplot)
% core_simulation Simulates the optical communication system and calculates the BER.
%
% Inputs:
%   X_CD - Received signal for X polarization with chromatic dispersion
%   Y_CD - Received signal for Y polarization with chromatic dispersion
%   r - Modulation index (1 for QPSK, 2 for 16QAM, 3 for 64QAM)
%   Rs - Symbol rate in GBaud
%   OSNR_dB - Optical Signal-to-Noise Ratio in dB
%   EQ_mode - Equalization mode (e.g., 'LMS')
%   EQ_N_tap - Number of taps for the equalizer
%   EQ_mu - Learning rate for the training phase of the equalizer
%   EQ_mu2 - Learning rate for the tracking phase of the equalizer
%   EQ_N1 - Number of symbols used for training
%   CarSync_DampFac - Damping factor for carrier synchronization (unused in this code)
%   scatplot - Flag to generate scatter plots (0 or 1)
%
% Outputs:
%   Ber_Tot - Total Bit Error Rate (BER) for the system

% Define modulation schemes and parameters
MODULATIONS = ["QPSK", "16QAM", "64QAM"];
modulation = ["QPSK", "QAM", "QAM"];
Baud_rate = num2str(Rs);

% Load the transmitted sequence based on modulation and baud rate
load(strcat('TXsequences/TXsequence_', MODULATIONS(r), '_', Baud_rate, 'GBaud.mat'));

% Repeat the transmitted bits and symbols 10 times to simulate the original transmission
TX_BITS_Xpol = repmat(SIG.Xpol.bits, 10, 1);
TX_BITS_Ypol = repmat(SIG.Ypol.bits, 10, 1);
TX_SYMB = [repmat(SIG.Xpol.txSymb, 10, 1), repmat(SIG.Ypol.txSymb, 10, 1)];

% Set modulation order and power normalization factor
if r == 1
    M = 4;
    power_norm = 2;
elseif r == 2
    M = 16;
    power_norm = 10;
else
    M = 64;
    power_norm = 42;
end

% Add AWGN noise to the signals
[X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB, SIG.symbolRate);
[Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB, SIG.symbolRate);

% Compensate for chromatic dispersion
[X_CD_rec, Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, SIG.Sps, 2);

% Perform frequency compensation
[X_freq_rec, Y_freq_rec] = freq_compensation(X_CD_rec, Y_CD_rec, SIG.Sps, SIG.symbolRate);

% Downsample the signals
X_freq_rec = downsample(X_freq_rec, 4);
Y_freq_rec = downsample(Y_freq_rec, 4);

% Normalize the signal power
X_Power = mean(abs(X_freq_rec).^2);
X_freq_rec_norm = X_freq_rec / sqrt(X_Power / power_norm);

Y_Power = mean(abs(Y_freq_rec).^2);
Y_freq_rec_norm = Y_freq_rec / sqrt(Y_Power / power_norm);

% Equalization using LMS if selected
if EQ_mode == 'LMS'
    [X_eq, Y_eq, e_X, e_Y] = LMS(X_freq_rec_norm, Y_freq_rec_norm, EQ_mu, EQ_mu2, EQ_N_tap, TX_SYMB, M, EQ_N1);
    N = 3; % Number of repetitions to drop before calculating the error, to avoid high BER rate because LMS did not converge yet
    X_eq = X_eq(N*65536 + 1:end - 100);
    Y_eq = Y_eq(N*65536 + 1:end - 100);

    % Normalize the equalized signal power
    X_Power = mean(abs(X_eq).^2);
    X_eq = X_eq / sqrt(X_Power / power_norm);
    Y_Power = mean(abs(Y_eq).^2);
    Y_eq = Y_eq / sqrt(Y_Power / power_norm);

    % Demap the equalized symbols to bits
    [X_demappedBits, X_demappedSymb, Y_demappedBits, Y_demappedSymb] = Demapping(X_eq, Y_eq, SIG.Xpol, M);

    % Calculate BER for X and Y polarizations
    X_BER = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits), :)) / (length(X_demappedBits) * log2(M));
    Y_BER = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits), :)) / (length(Y_demappedBits) * log2(M));
else
    % If no equalization is selected, set BER to a high value (not implemented in this code)
    X_BER = 1;
    Y_BER = 1;
end

% Calculate total BER as the average of X and Y BERs
Ber_Tot = 0.5 * (min(X_BER) + min(Y_BER));

%--------------------------------------------------------------------------
% Scatter plot generation
if scatplot ~= 0
    num_bins = 50; % Number of bins for histogram
    num_points = 1e4; % Number of points for scatter plot
    jumps = round(0.5 * length(X_eq) / num_points); % Step size for plotting
    axlim = sqrt(M) + 1; % Axis limits for scatter plot
    pointsize = 10; % Size of points in scatter plot

    figure;

    % Scatter plot for CD compensated input
    subplot(2, 2, 1);
    [density, ~, ~, binX, binY] = histcounts2(real(X_distorted_AWGN(1:16*jumps:end)), imag(X_distorted_AWGN(1:16*jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_distorted_AWGN(1:16*jumps:end)), imag(X_distorted_AWGN(1:16*jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('CD comp. input');
    xlabel('X');
    ylabel('Y');
    grid on;
    axis square;
    hold off;

    % Scatter plot for frequency compensated output
    subplot(2, 2, 2);
    [density, ~, ~, binX, binY] = histcounts2(real(X_freq_rec_norm(1:4*jumps:end)), imag(X_freq_rec_norm(1:4*jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_freq_rec_norm(1:4*jumps:end)), imag(X_freq_rec_norm(1:4*jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('Freq. comp. output');
    xlabel('X');
    ylabel('Y');
    grid on;
    axis square;
    hold off;

    % Scatter plot for equalized output
    subplot(2, 2, 3);
    [density, ~, ~, binX, binY] = histcounts2(real(X_eq(1:2*jumps:end)), imag(X_eq(1:2*jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_eq(1:2*jumps:end)), imag(X_eq(1:2*jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('EQ output');
    xlabel('X');
    ylabel('Y');
    grid on;
    axis square;
    hold off;

    % Scatter plot for phase compensated output
    subplot(2, 2, 4);
    [density, ~, ~, binX, binY] = histcounts2(real(X_eq(round(end/2):jumps:end)), imag(X_eq(round(end/2):jumps:end)), [num_bins num_bins]);
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_eq(round(end/2):jumps:end)), imag(X_eq(round(end/2):jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0, 0], [-10, 10], 'k');
    plot([-10, 10], [0, 0], 'k');
    xlim([-axlim, axlim]);
    ylim([-axlim, axlim]);
    colormap(jet);
    title('Phase comp. output');
    xlabel('X');
    ylabel('Y');
    grid on;
    axis square;
    hold off;

    % Add a title to the entire figure
    sgtitle(sprintf('%s, BER=%0.1e', MODULATIONS(r), Ber_Tot));
end
end
