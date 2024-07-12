function GUI_main_scatterplot(r, Rs, OSNR_dB, delta_nu, rad_sec, f_offset, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac)
% GUI_main_scatterplot Generates scatter plots of the optical communication system signals.
%
% Inputs:
%   r - Modulation index (1 for QPSK, 2 for 16QAM, 3 for 64QAM)
%   Rs - Symbol rate in GBaud
%   OSNR_dB - Optical Signal-to-Noise Ratio in dB
%   delta_nu - Laser linewidth
%   rad_sec - Phase noise
%   f_offset - Frequency offset
%   EQ_mode - Equalizer mode (e.g., 'LMS')
%   EQ_N_tap - Number of taps for the equalizer
%   EQ_mu - Learning rate for the equalizer training phase
%   EQ_mu2 - Learning rate for the equalizer tracking phase
%   EQ_N1 - Number of symbols used for training
%   CarSync_DampFac - Damping factor for carrier synchronization (unused in this code)
%
% Outputs:
%   None. Generates scatter plots of the signals.

% Define modulation schemes and parameters
MODULATIONS = ["QPSK", "16QAM", "64QAM"];
modulation = ["QPSK", "QAM", "QAM"];
Baud_rate = num2str(Rs);

% Load the transmitted sequence based on modulation and baud rate
load(strcat('TXsequences/TXsequence_', MODULATIONS(r), '_', Baud_rate, 'GBaud.mat'));

% Repeat the transmitted bits 10 times to simulate the original transmission
TX_BITS_Xpol = repmat(SIG.Xpol.bits, 10, 1);
TX_BITS_Ypol = repmat(SIG.Ypol.bits, 10, 1);

% Create delay and phase distorted signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate, f_offset);

% Add chromatic dispersion
[X_CD, Y_CD] = Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);

%% SIMULATION
% Run the core simulation and generate scatter plots
BER_Tot = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac, 1);

end
