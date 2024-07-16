function GUI_main_deltaSNR(r, Rs, start_sweep, end_sweep, points_to_sweep, log_or_lin, value2sweep, limit_while, BER_goal, tol, delta_nu, rad_sec, f_offset, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, EQ_N2, CarSync_DampFac,savefigure,prjcname,lbl)
% GUI_main_deltaSNR Simulates the effect of various parameters on OSNR penalty at a given BER goal.
%
% Inputs:
%   r - Modulation index (1 for QPSK, 2 for 16QAM, 3 for 64QAM)
%   Rs - Symbol rate in GBaud
%   start_sweep - Start value for the parameter sweep
%   end_sweep - End value for the parameter sweep
%   points_to_sweep - Number of points to sweep between start and end values
%   log_or_lin - Sweep type: 'log' for logarithmic, 'lin' for linear
%   value2sweep - Parameter to sweep (e.g., 'delta_nu', 'rad_sec', etc.)
%   limit_while - Maximum number of iterations for while loop convergence
%   BER_goal - Bit Error Rate goal
%   tol - Tolerance for BER convergence
%   delta_nu - Laser linewidth for initial parameter
%   rad_sec - Phase noise for initial parameter
%   f_offset - Frequency offset for initial parameter
%   EQ_mode - Equalizer mode (e.g., 'LMS')
%   EQ_N_tap - Number of taps for the equalizer
%   EQ_mu - Learning rate for the equalizer training phase
%   EQ_mu2 - Learning rate for the equalizer tracking phase
%   EQ_N1 - Number of symbols used for training
%   CarSync_DampFac - Damping factor for carrier synchronization (unused in this code)
%
% Outputs:
%   None. Generates a plot of OSNR penalty vs. the swept parameter.

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

% Set modulation order and power normalization factor
if r == 1
    M = 4;
    power_norm = 2;
    SNR_opt = 10 * log10(2 * erfinv(1 - 2 * BER_goal)^2);
elseif r == 2
    M = 16;
    power_norm = 10;
    SNR_opt = 10 * log10(10 * erfinv(1 - 8/3 * BER_goal)^2);
else
    M = 64;
    power_norm = 42;
    SNR_opt = 10 * log10(42 * erfinv(1 - 24/7 * BER_goal)^2);
end

% Repeat the transmitted bits 10 times to simulate the original transmission
TX_BITS_Xpol = repmat(SIG.Xpol.bits, 10, 1);
TX_BITS_Ypol = repmat(SIG.Ypol.bits, 10, 1);

% Generate sweep values either logarithmically or linearly
if log_or_lin == 'log'
    sweep_values = logspace(log10(start_sweep), log10(end_sweep), points_to_sweep);
else
    sweep_values = linspace(start_sweep, end_sweep, points_to_sweep);
end

% Initialize Delta_SNR array
Delta_SNR = zeros(1, points_to_sweep);

%% SIMULATION
for index = 1:points_to_sweep
    % Update the parameter to be swept
    switch value2sweep
        case 'delta_nu'
            delta_nu = sweep_values(index);
        case 'rad_sec'
            rad_sec = sweep_values(index);
        case 'freq_offset'
            f_offset = sweep_values(index);
        case 'EQ_N_tap'
            EQ_N_tap = sweep_values(index);
        case 'EQ_mu'
            EQ_mu = sweep_values(index);
        case 'EQ_mu2'
            EQ_mu2 = sweep_values(index);
        case 'EQ_N1'
            EQ_N1 = sweep_values(index);
        case 'CarSync_DampFac'
            CarSync_DampFac = sweep_values(index);
    end

    %% IMPAIRMENTS PART
    % Create delay and phase distorted signals
    [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate, f_offset);

    % Add chromatic dispersion
    [X_CD, Y_CD] = Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1,SIG.symbolRate);

    % Initialize variables for the while loop
    cycle = 0;           % Number of iterations
    OSNR_calc = 0;       % OSNR adjustment value
    BER_Tot = 10;        % Initial BER (set to a high value to start the loop)
    OSNR_dB = SNR_opt;   % Initial OSNR value set to the optimal SNR for the modulation

    % While loop for OSNR adjustment to meet the BER goal within tolerance
    while (round(BER_Tot / BER_goal, 5) >= 1 + (tol / 100) || round(BER_Tot / BER_goal, 5) <= 1 - (tol / 100)) && (cycle < limit_while)
        cycle = cycle + 1;               % Increment the cycle count
        lbl.Text = sprintf('Sweep point %d/%d, Iteration %d/%d',index,points_to_sweep,cycle,limit_while);
        drawnow;
        OSNR_dB = OSNR_dB + OSNR_calc;   % Adjust the OSNR value

        % Run the core simulation to get the BER for the current OSNR
        BER_Tot = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, EQ_N2, CarSync_DampFac, [0,0,0,0,0],0);

        % Calculate the OSNR adjustment based on the BER and modulation scheme
        if r == 1
            % QPSK modulation
            OSNR_inv = 10 * log10(2 * erfinv(1 - 2 * BER_Tot)^2);    % Inverse OSNR calculation
            OSNR_calc = SNR_opt - OSNR_inv;                         % OSNR adjustment
        elseif r == 2
            % 16QAM modulation
            if round(BER_Tot - BER_goal, 5) >= 9e-4 && round(BER_Tot - BER_goal, 5) <= 9e-3
                OSNR_calc = 1.5;   % Small adjustment step
            elseif round(BER_Tot - BER_goal, 5) >= 9e-3
                OSNR_calc = 4.5;   % Larger adjustment step
            else
                OSNR_inv = 10 * log10(10 * erfinv(1 - 8/3 * BER_Tot)^2);   % Inverse OSNR calculation
                OSNR_calc = SNR_opt - OSNR_inv;                          % OSNR adjustment
            end
        else
            % 64QAM modulation
            if round(BER_Tot - BER_goal, 5) >= 9e-4 && round(BER_Tot - BER_goal, 5) <= 9e-3
                OSNR_calc = 1.5;   % Small adjustment step
            elseif round(BER_Tot - BER_goal, 5) >= 9e-3
                OSNR_calc = 4.5;   % Larger adjustment step
            else
                OSNR_inv = 10 * log10(42 * erfinv(1 - 24/7 * BER_Tot)^2);  % Inverse OSNR calculation
                OSNR_calc = SNR_opt - OSNR_inv;                          % OSNR adjustment
            end
        end
    end

    % Store the OSNR penalty
    Delta_SNR(index) = OSNR_dB - SNR_opt;

    % Check for convergence
    if cycle == limit_while
        Delta_SNR(index) = NaN;
        fprintf('Too much time to convergence, OSNR penalty too large\n');
    else
        fprintf('WHILE converged\n');
    end
end

%% PLOT RESULTS
figure;
if log_or_lin == 'log'
    semilogx(sweep_values, Delta_SNR, 'Color', 'r', 'LineWidth', 2);
else
    plot(sweep_values, Delta_SNR, 'Color', 'r', 'LineWidth', 2);
end
title(sprintf('%s OSNR penalty at BER=%.0d', MODULATIONS(r), BER_goal));
xlabel(value2sweep);
ylabel('OSNR penalty [dB]');
axis tight;
ylim([0, max(Delta_SNR)]);
grid on;
if savefigure == 1
    if ~exist(prjcname, 'dir')
        mkdir(prjcname);
    end
    savefig(fullfile(prjcname,'OPTCOH_SNR_penalty_plot'));
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
    writetable(T, fullfile(prjcname,'OPTCOH_parameters.txt'), 'Delimiter', '\t');

    T = table(sweep_values', Delta_SNR', 'VariableNames', {value2sweep, 'Delta SNR'});
    % Write the table to a text file
    writetable(T, fullfile(prjcname,'OPTCOH_Delta_SNR.txt'), 'Delimiter', '\t');
end
end
