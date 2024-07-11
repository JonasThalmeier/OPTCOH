function GUI_main_deltaSNR(r,Rs, start_sweep, end_sweep, points_to_sweep, log_or_lin, value2sweep, limit_while, BER_goal, tol, delta_nu, rad_sec, f_offset, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac)
% GUI_main_deltaSNR(1,64, 0, 100e3, 4, 'lin', 'delta_nu', 10, 10e-4, 50e-3, 0, 8, 8e-3, 8e-4, 1e3,3e5,150);
MODULATIONS = ["QPSK","16QAM","64QAM"];
modulation = ["QPSK","QAM","QAM"];
Baud_rate = num2str(Rs);
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));

if r == 1
    M = 4;
    power_norm = 2;
    SNR_opt = 10*log10(2*erfinv(1-2*BER_goal)^2);
elseif r==2
    M = 16;
    power_norm = 10;
    SNR_opt = 10*log10(10*erfinv(1-8/3*BER_goal)^2);
else
    M = 64;
    power_norm = 42;
    SNR_opt = 10*log10(42*erfinv(1-24/7*BER_goal)^2);
end

TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission

if log_or_lin == 'log'
    sweep_values = logspace(log10(start_sweep), log10(end_sweep), points_to_sweep);
else
    sweep_values = linspace(start_sweep, end_sweep, points_to_sweep);
end

OSNR_dB = SNR_opt;
Delta_SNR = zeros(1,points_to_sweep);

%% SIUMULATION
for index = 1:points_to_sweep

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
    % Create delay and phase convolved signals
    [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate, f_offset);
    %add chromatic dispersion
    [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);


    %% CMA mine

    cycle = 0;
    OSNR_calc = 0;
    BER_Tot = 10;
    while (round(BER_Tot/BER_goal,5)>=1+(tol/100) || round(BER_Tot/BER_goal,5)<=1-(tol/100)) && (cycle<limit_while)
        cycle = cycle+1;
        OSNR_dB = OSNR_dB + OSNR_calc;
        BER_Tot = core_simulation(X_CD,Y_CD,r,Rs, OSNR_dB, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac,0);
        if r==1
            OSNR_inv =  10*log10(2*erfinv(1-2*BER_Tot)^2);
            OSNR_calc = SNR_opt - OSNR_inv;
        elseif r==2
            if round(BER_Tot-BER_goal,5)>=9e-4 && round(BER_Tot-BER_goal,5)<=9e-3
                OSNR_calc = 1.5;
            elseif round(BER_Tot-BER_goal,5)>=9e-3
                OSNR_calc = 4.5;
            else
                OSNR_inv =  10*log10(10*erfinv(1-8/3*BER_Tot)^2);
                OSNR_calc = SNR_opt - OSNR_inv;
            end
        else
            if round(BER_Tot-BER_goal,5)>=9e-4 && round(BER_Tot-BER_goal,5)<=9e-3
                OSNR_calc = 1.5;
            elseif round(BER_Tot-BER_goal,5)>=9e-3
                OSNR_calc = 4.5;
            else
                OSNR_inv =  10*log10(42*erfinv(1-24/7*BER_Tot)^2);
                OSNR_calc = SNR_opt - OSNR_inv;
            end
        end
    end


    Delta_SNR(index) = OSNR_dB - SNR_opt;

    if cycle==limit_while
        Delta_SNR(index) = NaN;
        fprintf('Too much time to convergence, OSNR penalty too large\n');
        % break;
    else
        fprintf('WHILE converged\n');
    end


end
%%

%------------------FIGURES-------------

figure;
if log_or_lin == 'log'
    semilogx(sweep_values(1,1:length(Delta_SNR)),Delta_SNR, 'Color', 'r', 'LineWidth',2);
else
    plot(sweep_values(1,1:length(Delta_SNR)),Delta_SNR, 'Color', 'r', 'LineWidth',2);
end
title(sprintf('%s OSNR penalty at BER=%.0d', MODULATIONS(r), BER_goal));
xlabel(value2sweep);
ylabel('OSNR penalty [dB]');
axis tight;
ylim([0,max(Delta_SNR)]);
grid on;

end

