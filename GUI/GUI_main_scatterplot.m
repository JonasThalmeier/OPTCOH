function GUI_main_scatterplot(r,Rs, OSNR_dB, delta_nu, rad_sec, f_offset, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac)
MODULATIONS = ["QPSK","16QAM","64QAM"];
modulation = ["QPSK","QAM","QAM"];
Baud_rate = num2str(Rs);
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));

TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission

% Create delay and phase convolved signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate, f_offset);
%add chromatic dispersion
[X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);

%% SIUMULATION
BER_Tot = core_simulation(X_CD,Y_CD,r,Rs, OSNR_dB, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac,1);
end



