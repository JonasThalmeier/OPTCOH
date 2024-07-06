function GUI_main_deltaSNR(r,Rs, start_sweep, end_sweep, points_to_sweep, log_or_lin, value2sweep, limit_while, BER_goal, delta_nu, rad_sec, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, EQ_N2, CarSync_DampFac)
% GUI_main_deltaSNR(1,64, 0, 100e3, 4, 'lin', 'delta_nu', 10, 10e-4, 50e-3, 0, 8, 8e-3, 8e-4, 1e3,3e5,150);
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "QAM"];
Baud_rate = num2str(Rs);
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));

if r == 1
    M = 4;
    power_norm = 2;
    SNR_opt = 10*log10(2*erfinv(1-2*BER_goal)^2);
else
    M = 16;
    power_norm = 10;
    SNR_opt = 10*log10(10*erfinv(1-8/3*BER_goal)^2);
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
        case 'EQ_N_tap'
            EQ_N_tap = sweep_values(index);
        case 'EQ_mu'
            EQ_mu = sweep_values(index);
        case 'EQ_mu2'
            EQ_mu2 = sweep_values(index);
        case 'EQ_N1'
            EQ_N1 = sweep_values(index);
        case 'EQ_N2'
            EQ_N2 = sweep_values(index);
        case 'CarSync_DampFac'
            CarSync_DampFac = sweep_values(index);
    end


    %% IMPAIRMENTS PART
    % Create delay and phase convolved signals
    [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate);
    %add chromatic dispersion
    [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);


    %% CMA mine

    cycle = 0;
    OSNR_calc = 0;
    BER_Tot = 10;

    while (round(BER_Tot-BER_goal,5)>=0.05e-3 || round(BER_Tot-BER_goal,5)<=(-0.05e-3)) && (cycle<limit_while)

        cycle = cycle+1;


        OSNR_dB = OSNR_dB + OSNR_calc;

        [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB, SIG.symbolRate);
        [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB, SIG.symbolRate);


        % ----------------Compensation for CD-------------------

        [X_CD_rec,Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, SIG.Sps, 2);

        X_CD_rec = X_CD_rec(65536*8+1:end);
        Y_CD_rec = Y_CD_rec(65536*8+1:end);

        X_CD_rec = downsample(X_CD_rec, 4);
        Y_CD_rec = downsample(Y_CD_rec, 4);

        if(rem(length(X_CD_rec),2) ~= 0)

            X_CD_rec = X_CD_rec(2:end);
            Y_CD_rec = Y_CD_rec(2:end);

        end

        X_Power = mean(abs((X_CD_rec)).^2);
        X_CD_rec_norm = X_CD_rec/sqrt(X_Power);

        Y_Power = mean(abs((Y_CD_rec)).^2);
        Y_CD_rec_norm = Y_CD_rec/sqrt(Y_Power);

        TX_sig = [X_CD_rec_norm, Y_CD_rec_norm];

        % %CMA parameters
        % if r==1
        %     N_tap = 8; %13
        %     mu = 8e-3; %7e-3
        %     mu2 = 8e-4;
        %     N1 = 1e3; %5000
        %     N2 = 3e5; %300000
        % else
        %     N_tap = 11; %13
        %     mu = 8e-3; %7e-3
        %     mu2 = 8e-4;
        %     N1 = 1e3; %5000
        %     N2 = 3e5; %300000
        % end

        [X_out, Y_out, e_X, e_Y] = EQ_func_N(TX_sig, r, EQ_mu, EQ_mu2, EQ_N_tap, EQ_N1, EQ_N2);

        cut = floor(EQ_N_tap/2);

        X_eq_CMA = X_out;
        Y_eq_CMA = Y_out;


        %------------------Delay&Phase recovery ---------------------
        if r==1
            % CarSync_DampFac = 150;
            carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r),"SamplesPerSymbol", 1, 'DampingFactor', CarSync_DampFac);
            [X_eq, phEstX] = carrSynch(X_eq_CMA);
            [Y_eq, phEstY] = carrSynch(Y_eq_CMA);
        else
            % CarSync_DampFac = 31.6;
            carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', CarSync_DampFac);
            %             carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 230, 'NormalizedLoopBandwidth',1e-3);
            %             [X_test, phEstX_test] = carrSynch(X_eq_CMA);
            [X_eq, phEstX] = carrSynch(X_eq_CMA);

            carrSynch2 = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', CarSync_DampFac);
            %             [Y_test, phEstY_test] = carrSynch(Y_eq_CMA);
            [Y_eq, phEstY] = carrSynch2(Y_eq_CMA);

            %             [X_eq,~] = BPS_N(X_eq_CMA, 50, M, power_norm);
            %             [Y_eq,~] = BPS_N(Y_eq_CMA, 50, M, power_norm);
        end


        X_BER = zeros(1,4);
        Y_BER = zeros(1,4);
        j=1;

        for i=0:pi/2:3/2*pi
            X_RX = X_eq*exp(1i*i);
            [~,transient_Xpol] = max(abs(xcorr(X_RX(1:65536*2), SIG.Xpol.txSymb)));
            X_RX = X_RX(transient_Xpol+1:end);
            Y_RX = Y_eq*exp(1i*i);
            [~,transient_Ypol] = max(abs(xcorr(Y_RX(1:65536*2), SIG.Ypol.txSymb)));
            Y_RX = Y_RX(transient_Ypol+1:end);

            if r==2
                X_Power = mean(abs((X_RX)).^2);
                X_RX = X_RX/sqrt(X_Power/power_norm);

                Y_Power = mean(abs((Y_RX)).^2);
                Y_RX = Y_RX/sqrt(Y_Power/power_norm);
            end

            if r==1
                [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);

            else

                [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_RX,Y_RX);
            end

            X_BER(j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
            Y_BER(j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));

            j=j+1;

        end

        X_Ber_Tot = min(X_BER);
        Y_Ber_Tot = min(Y_BER);

        if r==1
            BER_Tot = (X_Ber_Tot+Y_Ber_Tot)./2;
            OSNR_inv =  10*log10(2*erfinv(1-2*BER_Tot)^2);
            OSNR_calc = SNR_opt - OSNR_inv;

        else
            BER_Tot = (X_Ber_Tot+Y_Ber_Tot)./2;
            if round(BER_Tot-BER_goal,5)>=9e-4 && round(BER_Tot-BER_goal,5)<=9e-3
                OSNR_calc = 1.5;
            elseif round(BER_Tot-BER_goal,5)>=9e-3
                OSNR_calc = 4.5;
            else
                OSNR_inv =  10*log10(10*erfinv(1-8/3*BER_Tot)^2);
                OSNR_calc = SNR_opt - OSNR_inv;
            end
        end

    end


    Delta_SNR(index) = OSNR_dB - SNR_opt;

    if cycle==limit_while
        Delta_SNR = Delta_SNR(1:index-1);
        fprintf('Too much time to convergence, OSNR penalty too large\n');
        break;
    else
        fprintf('WHILE converged\n');
    end

end
%%

%------------------FIGURES-------------

figure;
semilogx(sweep_values(1,1:length(Delta_SNR)),Delta_SNR, 'Color', 'r', 'LineWidth',2);
title(sprintf('%s OSNR penalty at BER=%.0d', MODULATIONS(r), BER_goal));
xlabel(value2sweep);
ylabel('OSNR penalty [dB]');
axis tight;
ylim([0,max(Delta_SNR)]);
grid on;

end

