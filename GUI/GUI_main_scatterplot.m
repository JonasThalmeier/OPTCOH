function GUI_main_scatterplot(r, Rs, OSNR_dB, delta_nu, rad_sec, f_offset, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, EQ_N2, CarSync_DampFac,scatter_vec,savefigure,prjcname)
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

% Construct the file name dynamically
fileName = sprintf('TXsequence_%s_%sGBaud.mat', MODULATIONS{r}, Baud_rate);
% Construct the full path to the .mat file
matFilePath = fullfile(fileparts(mfilename('fullpath')), '..', 'TXsequences', fileName);
% Load the .mat file
load(matFilePath);

% Repeat the transmitted bits 10 times to simulate the original transmission
TX_BITS_Xpol = repmat(SIG.Xpol.bits, 10, 1);
TX_BITS_Ypol = repmat(SIG.Ypol.bits, 10, 1);

% Create delay and phase distorted signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate, f_offset);

% Add chromatic dispersion
[X_CD, Y_CD] = Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1, SIG.symbolRate);

%% SIMULATION
% Run the core simulation and generate scatter plots
BER_Tot = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, EQ_N2, CarSync_DampFac, scatter_vec,savefigure,prjcname);
if savefigure == 1
    % Convert numeric values to strings
    param_values = {MODULATIONS(r), num2str(Rs), num2str(OSNR_dB), num2str(delta_nu), ...
        num2str(rad_sec), num2str(f_offset), EQ_mode, num2str(EQ_N_tap), ...
        num2str(EQ_mu), num2str(EQ_mu2), num2str(EQ_N1), num2str(EQ_N2), ...
        num2str(CarSync_DampFac)};

    % Define parameter names
    param_names = {'Modulation', 'Baud rate [GHz]', 'OSNR [dB]', 'delta nu [Hz]', ...
        'Pol. rotation [rad/sec]', 'freq. offset [Hz]', 'EQ mode', ...
        'EQ Num of Taps', 'EQ mu 1', 'EQ mu 2', 'EQ N1', 'EQ N2', ...
        'Damping factor Carrier Sync'};

    % Create a table
    T = table(param_names', param_values', 'VariableNames', {'Parameter', 'Value'});

    % Write the table to a text file
    writetable(T, 'OPTCOH_parameters.txt', 'Delimiter', '\t');
end
end
