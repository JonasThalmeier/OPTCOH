function GUI_main_BER(r, Rs, start_sweep, end_sweep, points_to_sweep, delta_nu, rad_sec, f_offset, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, EQ_N2, CarSync_DampFac,savefigure,prjcname,lbl)
% GUI_main_BER Simulates the BER vs. OSNR curve for different modulation schemes.
%
% Inputs:
%   r - Modulation index (1 for QPSK, 2 for 16QAM, 3 for 64QAM)
%   Rs - Symbol rate in GBaud
%   start_sweep - Start value for OSNR sweep
%   end_sweep - End value for OSNR sweep
%   points_to_sweep - Number of points to sweep between start and end values
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
%   None. Generates a plot of BER vs. OSNR.

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

% Generate fine-grained SNR values for theoretical BER calculation
SNR_fine = linspace(start_sweep, end_sweep, 1000);

% Set modulation-specific parameters and calculate theoretical BER
if r == 1
    % QPSK modulation
    M = 4;
    power_norm = 2;
    BER_TH = 0.5 * erfc(sqrt(10.^(SNR_fine / 10) / 2));
elseif r == 2
    % 16QAM modulation
    M = 16;
    power_norm = 10;
    BER_TH = 3/8 * erfc(sqrt(10.^(SNR_fine / 10) / 10));
else
    % 64QAM modulation
    M = 64;
    power_norm = 42;
    BER_TH = 7/24 * erfc(sqrt(10.^(SNR_fine / 10) / 42));
end

% Repeat the transmitted bits 10 times to simulate the original transmission
TX_BITS_Xpol = repmat(SIG.Xpol.bits, 10, 1);
TX_BITS_Ypol = repmat(SIG.Ypol.bits, 10, 1);

% Generate OSNR values for simulation
OSNR_dB = linspace(start_sweep, end_sweep, points_to_sweep);

% Initialize BER array for simulation results
Ber_Tot = zeros(1, points_to_sweep);

% Create delay and phase distorted signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate, f_offset);

% Add chromatic dispersion
[X_CD, Y_CD] = Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1, SIG.symbolRate);

%% SIMULATION
% Run the core simulation for each OSNR value
clc
for index = 1:points_to_sweep
    lbl.Text = sprintf('Simulation of OSNR value %d/%d',index,points_to_sweep);
    drawnow;
    Ber_Tot(index) = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB(index), EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, EQ_N2, CarSync_DampFac, [0,0,0,0,0],0,prjcname,'BER');
end

%% PLOT RESULTS
% Plot theoretical and simulated BER vs. OSNR curves
figure();
semilogy(SNR_fine, BER_TH, 'r', 'LineWidth', 1);
xlim([min(OSNR_dB), max(OSNR_dB)]);
grid on;
hold on;
semilogy(OSNR_dB, Ber_Tot, 'Marker', 'o', 'Color', "#77AC30", 'LineWidth', 1, 'LineStyle', '-.');
title(sprintf('%s BER curve', MODULATIONS(r)));
legend('Theoretical BER', 'Simulated BER', 'Interpreter', 'latex');
xlabel('OSNR [dB]', 'Interpreter', 'latex');
ylabel('BER', 'Interpreter', 'latex');
hold off;
if savefigure == 1
    if ~exist(prjcname, 'dir')
        mkdir(prjcname);
    end
    savefig(fullfile(prjcname,'OPTCOH_BER_plot'));
    % Convert numeric values to strings
    param_values = {MODULATIONS(r), num2str(Rs), num2str(delta_nu), ...
        num2str(rad_sec), num2str(f_offset), EQ_mode, num2str(EQ_N_tap), ...
        num2str(EQ_mu), num2str(EQ_mu2), num2str(EQ_N1), num2str(EQ_N2), ...
        num2str(CarSync_DampFac)};
    % Define parameter names
    param_names = {'Modulation', 'Baud rate [GHz]', 'delta nu [Hz]', ...
        'Pol. rotation [rad/sec]', 'freq. offset [Hz]', 'EQ mode', ...
        'EQ Num of Taps', 'EQ mu 1', 'EQ mu 2', 'EQ N1', 'EQ N2', ...
        'Damping factor Carrier Sync'};
    % Create a table
    T = table(param_names', param_values', 'VariableNames', {'Parameter', 'Value'});
    % Write the table to a text file
    writetable(T, fullfile(prjcname,'OPTCOH_parameters.txt'), 'Delimiter', '\t');
    % Create a table
    T = table(OSNR_dB', BER', 'VariableNames', {'OSNR [dB]', 'BER'});
    % Write the table to a text file
    writetable(T, fullfile(prjcname,'OPTCOH_BER.txt'), 'Delimiter', '\t');
end

end
