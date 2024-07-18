function [Ber_Tot] = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, EQ_N2, CarSync_DampFac, scatter_vec, savefigure,prjcname,GUIname)
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

% Construct the file name dynamically
fileName = sprintf('TXsequence_%s_%sGBaud.mat', MODULATIONS{r}, Baud_rate);
% Construct the full path to the .mat file
matFilePath = fullfile(fileparts(mfilename('fullpath')), '..', 'TXsequences', fileName);
% Load the .mat file
load(matFilePath);

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
[X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB);
[Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB);

% Downsample at 2SpS
X_distorted_AWGN = downsample(X_distorted_AWGN, 4);
Y_distorted_AWGN = downsample(Y_distorted_AWGN, 4);

% Compensate for chromatic dispersion
[X_CD_rec,Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, 2, 2, SIG.symbolRate);

% Frequency compensation
[X_freq_rec, Y_freq_rec] = freq_compensation(X_CD_rec, Y_CD_rec, 2, SIG.symbolRate);

% Equalization using LMS if selected
if isequal(EQ_mode, 'LMS')

    % Normalize the signal power
    X_Power = mean(abs(X_freq_rec).^2);
    X_freq_rec_norm = X_freq_rec / sqrt(X_Power / power_norm);

    Y_Power = mean(abs(Y_freq_rec).^2);
    Y_freq_rec_norm = Y_freq_rec / sqrt(Y_Power / power_norm);

    [X_eq, Y_eq, e_X, e_Y] = LMS(X_freq_rec_norm, Y_freq_rec_norm, EQ_mu, EQ_mu2, EQ_N_tap, TX_SYMB, M, EQ_N1,GUIname,r,Rs);
    N = 3; % Number of repetitions to drop before calculating the error, to avoid high BER rate because LMS did not converge yet
    X_eq = X_eq(N*65536 + 1:end - 100);
    Y_eq = Y_eq(N*65536 + 1:end - 100);

    X_sync = 0;
    Y_sync = 0;

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

elseif isequal(EQ_mode, 'CMA/RDE')

    X_freq_rec = X_freq_rec(65536*2+1:end);
    Y_freq_rec = Y_freq_rec(65536*2+1:end);

    if(rem(length(X_freq_rec),2) ~= 0)

        X_freq_rec = X_freq_rec(2:end);
        Y_freq_rec = Y_freq_rec(2:end);

    end

    % Normalize the signal power at unit power
    X_Power = mean(abs((X_freq_rec)).^2);
    X_CD_rec_norm = X_freq_rec/sqrt(X_Power);

    Y_Power = mean(abs((Y_freq_rec)).^2);
    Y_CD_rec_norm = Y_freq_rec/sqrt(Y_Power);

    TX_sig = [X_CD_rec_norm, Y_CD_rec_norm];

    [X_eq, Y_eq, e_X, e_Y] = CMA_RDE(TX_sig, r, EQ_mu, EQ_mu2, EQ_N_tap, EQ_N1, EQ_N2,GUIname,Rs);

    %------------------Delay&Phase recovery ---------------------
    if r==1
        carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r),"SamplesPerSymbol", 1, 'DampingFactor', CarSync_DampFac);
        [X_sync, phEstX] = carrSynch(X_eq);
        [Y_sync, phEstY] = carrSynch(Y_eq);
    else
        carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', CarSync_DampFac);
        [X_sync, phEstX] = carrSynch(X_eq);

        carrSynch2 = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', CarSync_DampFac);
        [Y_sync, phEstY] = carrSynch2(Y_eq);
    end
    X_BER_int = zeros(1,4);
    Y_BER_int = zeros(1,4);
    j=1;

    for i=0:pi/2:3/2*pi
        X_RX = X_sync*exp(1i*i);
        [~,transient_Xpol] = max(abs(xcorr(X_RX(1:65536*2), SIG.Xpol.txSymb)));
        X_RX = X_RX(transient_Xpol+1:end);

        Y_RX = Y_sync*exp(1i*i);
        [~,transient_Ypol] = max(abs(xcorr(Y_RX(1:65536*2), SIG.Ypol.txSymb)));
        Y_RX = Y_RX(transient_Ypol+1:end);

        % Normalize the equalized signal power
        X_Power = mean(abs(X_RX).^2);
        X_RX = X_RX / sqrt(X_Power / power_norm);
        Y_Power = mean(abs(Y_RX).^2);
        Y_RX = Y_RX / sqrt(Y_Power / power_norm);

        [X_demappedBits, X_demappedSymb, Y_demappedBits, Y_demappedSymb] = Demapping(X_RX, Y_RX, SIG.Xpol, M);

        X_BER_int(j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
        Y_BER_int(j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));

        j=j+1;
    end

    X_BER = min(X_BER_int);
    Y_BER = min(Y_BER_int);

else
    % If no equalization is selected, set BER to a high value (not implemented in this code)
    X_BER = 1;
    Y_BER = 1;
end

% Calculate total BER as the average of X and Y BERs
Ber_Tot = 0.5 * (min(X_BER) + min(Y_BER));

%--------------------------------------------------------------------------
if sum(scatter_vec)~=0
    scatterplot_func(X_distorted_AWGN,Y_distorted_AWGN,X_CD_rec,Y_CD_rec,X_freq_rec,Y_freq_rec,X_eq,Y_eq,X_sync,Y_sync,scatter_vec,prjcname,savefigure,M);
end
end
