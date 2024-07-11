function GUI_main_BER(r,Rs, start_sweep, end_sweep, points_to_sweep, delta_nu, rad_sec, f_offset, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac)
MODULATIONS = ["QPSK","16QAM","64QAM"];
modulation = ["QPSK","QAM","QAM"];
Baud_rate = num2str(Rs);
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));

SNR_fine = linspace(start_sweep,end_sweep,1000);
if r == 1
    M = 4;
    power_norm = 2;
    BER_TH = 0.5 * erfc(sqrt(10.^(SNR_fine/10)/2));
elseif r==2
    M = 16;
    power_norm = 10;
    BER_TH = 3/8 * erfc(sqrt(10.^(SNR_fine/10)/10));
else
    M = 64;
    power_norm = 42;
    BER_TH = 7/24 * erfc(sqrt(10.^(OSNR_dB/10)/42));
end

TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission


OSNR_dB = linspace(start_sweep, end_sweep, points_to_sweep);
Ber_Tot = zeros(1,points_to_sweep);

% Create delay and phase convolved signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate,f_offset);
%add chromatic dispersion
[X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);



%% SIUMULATION
for index = 1:points_to_sweep
    Ber_Tot(index) = core_simulation(X_CD,Y_CD,r,Rs, OSNR_dB(index), EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac,0);
end
%%

%------------------FIGURES-------------
figure();
semilogy(SNR_fine,BER_TH, 'r','LineWidth', 1);
xlim([min(OSNR_dB),max(OSNR_dB)]);
grid on;
hold on;
semilogy(OSNR_dB,Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineWidth', 1, 'LineStyle','-.');
title(sprintf('%s BER curve',MODULATIONS(r)));
legend('Theoretical BER', 'Simulated BER - CMA', 'Interpreter', 'latex');
xlabel('OSNR [dB]', 'Interpreter','latex');
hold off;

end

