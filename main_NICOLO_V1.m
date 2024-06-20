clear;
close all;
clc;

% Load the .mat file
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "QAM"];
Baud_rate = '128';
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
r = 2;
fprintf('The transmitted moduluation is: %s\n', modulation(r));
load(strcat('C:\Users\utente\Documents\GitHub\OPTCOH\TXsequences\TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));

%% PARAMETERS

BER_goal = 1e-3;
points_to_sweep = 10;
limit_while = 10;

if r == 1
    M = 4;
    power_norm = 2;
    SNR_opt =10*log10(2*erfinv(1-2*BER_goal)^2);
else
    M = 16;
    power_norm = 10;
    SNR_opt = 10*log10(10*erfinv(1-8/3*BER_goal)^2);
end

TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission


pol_sweep = logspace(1, 5, points_to_sweep);
lw_sweep = logspace(4,6,points_to_sweep);

disp(pol_sweep);

sweep_vector = [pol_sweep; lw_sweep];
standard_values = [50e3, 1e3]; %values to assign to the rotation not used

choise_rot = 1; % 1 or 2 to choose which sweep to apply (1-POL | 2-PHASE | >3-NO SWEEP)

if choise_rot<3 
    OSNR_dB = SNR_opt;
elseif choise_rot>=3 && r==1
    OSNR_dB = 3:11;
    points_to_sweep = 1;
else
    OSNR_dB = 17:20;
    points_to_sweep = 1;
end


Delta_Pol = zeros(1,points_to_sweep);

clear pol_sweep lw_sweep;

%% SIUMULATION
for index_rad_pol = 1:points_to_sweep

    %% IMPAIRMENTS PART
    % Create delay and phase convolved signals
    if choise_rot==1
        [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, standard_values(choise_rot), sweep_vector(choise_rot,index_rad_pol), SIG.symbolRate);
    elseif choise_rot==2
        [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, sweep_vector(choise_rot,index_rad_pol), standard_values(choise_rot), SIG.symbolRate);
    else
        [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, standard_values(1), standard_values(2), SIG.symbolRate);
    end

    %add chromatic dispersion
    [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);
    
    %%
    % Adding the noise
    
    X_Ber_Tot = zeros(1,length(OSNR_dB));
    Y_Ber_Tot = zeros(1,length(OSNR_dB));
    
    %% CMA mine

    index = 0;
    cycle = 0;
    OSNR_calc = 0;
    BER_Tot = 10;
    
    while (index<length(OSNR_dB) && (round(BER_Tot-BER_goal,5)>=0.05e-3 || round(BER_Tot-BER_goal,5)<=(-0.05e-3))) && (cycle<limit_while)

        cycle = cycle+1;

        if choise_rot<3
            index = 1;
            OSNR_dB(index) = OSNR_dB(index) + OSNR_calc;
        else
            index = index + 1;
        end
    
        [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
        [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
        
    
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
        
        %CMA parameters
        if r==1
            N_tap = 8; %13
            mu = 8e-3; %7e-3
            mu2 = 8e-4;
            N1 = 1e3; %5000
            N2 = 3e5; %300000
        else
            N_tap = 11; %13
            mu = 8e-3; %7e-3
            mu2 = 8e-4;
            N1 = 1e3; %5000
            N2 = 3e5; %300000
        end
    
        [X_out, Y_out, e_X, e_Y] = EQ_func_N(TX_sig, r, mu, mu2, N_tap, N1, N2);
        
        
        cut = floor(N_tap/2);
    
    
        if ((index==1 && index_rad_pol==1) || (index==length(OSNR_dB) && index_rad_pol==points_to_sweep))
        figure(), plot(e_X), title(sprintf('CMA error Xpol %d dB', OSNR_dB(index))), hold on, xline(cut, 'LineStyle','--', 'Color','r'), xlabel('N_samples');
        figure(), plot(e_Y), title(sprintf('CMA error Ypol %d dB', OSNR_dB(index))), hold on, xline(cut, 'LineStyle','--', 'Color','r'), xlabel('N_samples');
        end    
        
        
        X_eq_CMA = X_out;
        Y_eq_CMA = Y_out;
    
    
     %------------------Delay&Phase recovery ---------------------
        
        if r==1
            carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r),"SamplesPerSymbol", 1, 'DampingFactor', 150);
            [X_eq, phEstX] = carrSynch(X_eq_CMA);
            [Y_eq, phEstY] = carrSynch(Y_eq_CMA);
        else
            carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 1.6);
%             carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 230, 'NormalizedLoopBandwidth',1e-3);

%             [X_test, phEstX_test] = carrSynch(X_eq_CMA);
            [X_eq, phEstX] = carrSynch(X_eq_CMA);
    
%             carrSynch2 = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 240, 'NormalizedLoopBandwidth',.0002);%, 'ModulationPhaseOffset','Custom', 'CustomPhaseOffset', -pi/7);
%             [Y_test, phEstY_test] = carrSynch(Y_eq_CMA);
            [Y_eq, phEstY] = carrSynch(Y_eq_CMA);
    
%             [X_eq,~] = BPS_N(X_eq_CMA, 50, M, power_norm); 
%             [Y_eq,~] = BPS_N(Y_eq_CMA, 50, M, power_norm);
        end
        
        
        if ((index==1 && index_rad_pol==1) || (index==length(OSNR_dB) && index_rad_pol==points_to_sweep))
            scatterplot(X_eq);
            title(sprintf('%s Xpol after EQ and phase recovery, OSNR=%d dB',MODULATIONS(r), OSNR_dB(index)));
        end
    
        X_BER = zeros(1,4);
        Y_BER = zeros(1,4);
        j=1;
        
        for i=0:pi/2:3/2*pi
        
    %     fprintf('---------The phase tried is (degrees): %d-----------\n', (mean(i) *180 /pi));
        
        % fprintf('The total phase recovered is (degrees): %d\n', (mean(phEstX+i) *180 /pi));
           
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
    
         if(index==length(OSNR_dB) && j==1)
            fprintf('Transient Xpol: %d\n', transient_Xpol)
            fprintf('Transient Ypol: %d\n', transient_Ypol)
        end
    
        if (index==length(OSNR_dB) && j==1)
            scatterplot(X_RX);
            title(sprintf('%s Xpol after delay recovery, OSNR=%d dB',MODULATIONS(r), OSNR_dB(index)));
        end
        
        if r==1
            %fprintf('The tracked moduluation is: QPSK\n');
              [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);
      
        else
            %fprintf('The tracked moduluation is: 16-QAM\n');
    %         MyConst = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];
    %         X_demappedBits = qamdemod(X_RX, M, MyConst, OutputType='bit', PlotConstellation=true);
    %         N = length(X_demappedBits)/4;
    %         X_demappedBits = reshape(X_demappedBits, 4, N).';
    %        
            [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_RX,Y_RX);
        end
        
          X_BER(j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
          Y_BER(j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));
        
           j=j+1;
        
        end
        
        X_Ber_Tot(index) = min(X_BER);
        Y_Ber_Tot(index) = min(Y_BER);

        fprintf('The BER on Xpol is: %.6f\n', X_Ber_Tot(index));
        fprintf('The BER on Ypol is: %.6f\n', Y_Ber_Tot(index));

        if choise_rot<3
            BER_Tot = (X_Ber_Tot+Y_Ber_Tot)./2;
            fprintf('The BER on Ypol is: %.6f\n', Y_Ber_Tot(index));
            OSNR_inv =  10*log10(10*erfinv(1-8/3*BER_Tot)^2);
            OSNR_calc = SNR_opt - OSNR_inv;
            index = 0;
        end
    end
    
    if choise_rot>=3
        BER_Tot = (X_Ber_Tot+Y_Ber_Tot)./2;
    else
        Delta_Pol(index_rad_pol) = OSNR_dB - SNR_opt;

        if cycle==limit_while
            Delta_Pol = Delta_Pol(1:index_rad_pol-1);
            fprintf('Too much time to convergence, OSNR penalty too large\n');
            break;
        else
            fprintf('WHILE converged\n');
        end
    end
end
%%

%------------------FIGURES-------------

if choise_rot==1
    figure(), semilogx(sweep_vector(1,1:length(Delta_Pol)), Delta_Pol, 'Color', 'r', 'LineWidth',2), title('OSNR penalty vs polarization rotation'), xlabel('log(rad/sec)'), grid on;

elseif choise_rot==2 
    figure(), semilogx(sweep_vector(2,1:length(Delta_Pol)),Delta_Pol, 'Color', '#7E2F8E', 'LineWidth',2), title('OSNR penalty vs phase rotation'), xlabel('log($\Delta\nu$)','Intepreter','latex'), grid on;
    
else
    if r==1
        BER_TH = 0.5 * erfc(sqrt(10.^(OSNR_dB/10)/2));
     
        figure();
        semilogy(OSNR_dB, BER_TH, 'r', 'LineWidth', 1);
        xlim([min(OSNR_dB),max(OSNR_dB)]);
        grid on;
        hold on;
        semilogy(OSNR_dB, X_Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineWidth', 1, 'LineStyle','-.');  
        title(sprintf('%s BER curve of Xpol',MODULATIONS(r)));
        legend('Theoretical BER', 'Simulated BER - CMA', 'Interpreter', 'latex');
        xlabel('OSNR [dB]', 'Interpreter','latex');
        hold off;
        
        figure();
        semilogy(OSNR_dB,BER_TH, 'r','LineWidth', 1); 
        xlim([min(OSNR_dB),max(OSNR_dB)]);
        grid on;
        hold on;
        semilogy(OSNR_dB,Y_Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineWidth', 1, 'LineStyle','-.');
        title(sprintf('%s BER curve of Ypol',MODULATIONS(r)));
        legend('Theoretical BER', 'Simulated BER - CMA', 'Interpreter', 'latex');
        xlabel('OSNR [dB]', 'Interpreter','latex');
        hold off;
    
        figure();
        semilogy(OSNR_dB, BER_TH, 'r', 'LineWidth', 1);
        xlim([min(OSNR_dB),max(OSNR_dB)]);
        grid on;
        hold on;
        semilogy(OSNR_dB, BER_Tot, 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
        title(sprintf('%s BER curve',MODULATIONS(r)));
        legend('Theoretical BER', 'Simulated BER - CMA', 'Interpreter', 'latex');
        xlabel('OSNR [dB]', 'Interpreter','latex');
        hold off;

    else

        BER_TH = 3/8 * erfc(sqrt(10.^(OSNR_dB/10)/10));

        figure();
        semilogy(OSNR_dB, BER_TH, 'r');
        xlim([min(OSNR_dB),max(OSNR_dB)]);
        grid on;
        hold on;
        semilogy(OSNR_dB, X_Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
        title(sprintf('%s BER curve of Xpol',MODULATIONS(r)));
        legend('Theoretical BER', 'Simulated BER - RDE', 'Interpreter', 'latex');
        xlabel('OSNR [dB]', 'Interpreter','latex');
        hold off;
        
        figure();
        semilogy(OSNR_dB,BER_TH, 'r');
        xlim([min(OSNR_dB),max(OSNR_dB)]);
        grid on;
        hold on;
        semilogy(OSNR_dB,Y_Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
        title(sprintf('%s BER curve of Ypol',MODULATIONS(r)));
        legend('Theoretical BER', 'Simulated BER - RDE', 'Interpreter', 'latex');
        xlabel('OSNR [dB]', 'Interpreter','latex');
        hold off;   
    
        figure();
        semilogy(OSNR_dB, BER_TH, 'r', 'LineWidth', 1);
        xlim([min(OSNR_dB),max(OSNR_dB)]);
        grid on;
        hold on;
        semilogy(OSNR_dB, BER_Tot, 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
        title(sprintf('%s BER curve',MODULATIONS(r)));
        legend('Theoretical BER', 'Simulated BER - RDE', 'Interpreter', 'latex');
        xlabel('OSNR [dB]', 'Interpreter','latex');
        hold off;
        
    end
end

