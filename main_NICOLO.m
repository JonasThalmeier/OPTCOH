clear;
close all;
clc;

% Load the .mat file
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "QAM"];
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
r = 2;
fprintf('The transmitted moduluation is: %s\n', modulation(r));
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_64GBaud.mat'));
if r == 1
    M = 4;
else
    M = 16;
end

TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission


% Create delay and phase convolved signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig);

%add chromatic dispersion
[X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);

% Adding the noise
if r==1
    OSNR_dB = 4:10;
else
    OSNR_dB = 11:16;
end

X_Ber_Tot = zeros(1,length(OSNR_dB));
Y_Ber_Tot = zeros(1,length(OSNR_dB));
%%
for index = 1:length(OSNR_dB)

    [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
    [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
    
    % ----------------Compensation for CD-------------------
    
    [X_CD_rec,Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, SIG.Sps, 2);
    
    %fprintf('isequal = %d\n', isequal(round(X_CD,8), round(SIG.Xpol.txSig, 8)));
    
%     X_CD_rec = SIG.Xpol.txSig;
%     Y_CD_rec = SIG.Ypol.txSig;
    %------------------ Matched Flitering ---------------------
    X_CD_rec = downsample(X_CD_rec, 4);
    Y_CD_rec = downsample(Y_CD_rec, 4);
    
    [X_matched,Y_matched] = Matched_filtering(X_CD_rec, Y_CD_rec, PulseShaping.b_coeff);

      
    if (index==length(OSNR_dB))
        scatterplot(X_matched(1:2:end));
        title(sprintf('%s constellation of Xpol after matched filter %d dB of WGN',MODULATIONS(r), OSNR_dB(index)));
    end

    %------------------Delay&Phase recovery ---------------------
    
    if r==1
        carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r),"SamplesPerSymbol", 1, 'DampingFactor', 150);
        [X_eq, phEstX] = carrSynch(X_matched(1:2:end));
        [Y_eq, phEstY] = carrSynch(Y_matched(1:2:end));
    else
        carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 100);%, 'ModulationPhaseOffset','Custom', 'CustomPhaseOffset', -pi/7);
        [X_eq, phEstX] = carrSynch(X_matched(1:2:end));
        [Y_eq, phEstY] = carrSynch(Y_matched(1:2:end));
    end
    
    X_Power = mean(abs((X_eq)).^2);
    X_eq = X_eq/sqrt(X_Power/10);

    Y_Power = mean(abs((Y_eq)).^2);
    Y_eq = Y_eq/sqrt(Y_Power/10);

    if (index==length(OSNR_dB))
        scatterplot(X_eq);
        title(sprintf('%s constellation of Xpol after phase recovery',MODULATIONS(r)));
    end

%     X_eq = X_matched(1:2:end);
    %
    X_BER = zeros(1,4);
    Y_BER = zeros(1,4);
    j=1;
    
    for i=0:pi/2:3/2*pi
    
    fprintf('---------The phase tried is (degrees): %d-----------\n', (mean(i) *180 /pi));
    
    % fprintf('The total phase recovered is (degrees): %d\n', (mean(phEstX+i) *180 /pi));
    
    transient_Xpol = abs(finddelay(X_eq(1:65536), SIG.Xpol.txSymb));
    transient_Ypol = abs(finddelay(Y_eq(1:65536), SIG.Ypol.txSymb));
    
    fprintf('Transient Xpol: %d\n', transient_Xpol)
    fprintf('Transient Ypol: %d\n', transient_Ypol)
    
    X_RX = X_eq*exp(1i*i);
    X_RX = X_RX(transient_Xpol+1:end);
    Y_RX = Y_eq*exp(1i*i);
    Y_RX = Y_RX(transient_Ypol+1:end);

    if (index==length(OSNR_dB) && j==1)
        scatterplot(X_RX);
        title(sprintf('%s constellation of Xpol after delay recovery',MODULATIONS(r)));
    end
    
    if r==1
        fprintf('The tracked moduluation is: QPSK\n');
          [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);
    %       X_demappedBits = pskdemod(X_RX,M, pi/4*7); %it doesn't demodulate in the same way as our function
    %       X_demappedBits = From_MATLAB_pskdemod(X_demappedBits);
    else
        fprintf('The tracked moduluation is: 16-QAM\n');
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
end

%% 
%---------------CMA-------------------

X_Ber_Tot_CMA = zeros(1,length(OSNR_dB));
Y_Ber_Tot_CMA = zeros(1,length(OSNR_dB));
%%
for index = 1:length(OSNR_dB)

    [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
    [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
    
    % ----------------Compensation for CD-------------------
    
    [X_CD_rec,Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, SIG.Sps, 2);
    

    %---------------------------EQ---------------------------------------------
    %From 15 numtaps are the best for qpsk, also for qam, but 60 and 120 seem
    %better. The reference tap for qpsk are 1,3,65,121.
  
    if r == 1

        X_CD_rec = downsample(X_CD_rec, 4);
        Y_CD_rec = downsample(Y_CD_rec, 4);

        constellation = pskmod(0:3, 4, pi/4);
        stepsize = 1e-3; %best result
        numtaps = 9;
        referencetap = (numtaps-1)/2;

        EQ = comm.LinearEqualizer('Algorithm', 'CMA', 'StepSize', stepsize,'NumTaps', numtaps, 'InputSamplesPerSymbol', 2, 'Constellation', constellation, 'ReferenceTap', referencetap);
        [X_matched,errX] = EQ(X_CD_rec(1:end-mod(length(X_CD_rec),2))); 
        [Y_matched,errY] = EQ(Y_CD_rec(1:end-mod(length(Y_CD_rec),2))); 

    else

        X_CD_rec_1 = downsample(X_CD_rec, 4);
        Y_CD_rec_1 = downsample(Y_CD_rec, 4);

        [X_matched_1,Y_matched_1] = Matched_filtering(SIG.Xpol.txSig(1:4:end), SIG.Ypol.txSig(1:4:end), PulseShaping.b_coeff);

        X_Power = mean(abs((X_CD_rec)).^2);
        X_CD_rec_norm = X_CD_rec/sqrt(X_Power/10);
    
%         Y_Power = mean(abs((Y_eq)).^2);
%         Y_eq = Y_eq/sqrt(Y_Power/10);

%         X_Power = mean(abs(real(X_CD_rec)).^2);
%         X_CD_rec_norm = X_CD_rec/X_Power;
        %X_CD_rec_norm = downsample(X_CD_rec_norm, 4);
               
        transient_Xpol = abs(finddelay(X_CD_rec_norm(5:8:8*65536), SIG.Xpol.txSymb));

        X_CD_rec_norm = X_CD_rec_norm(transient_Xpol*8+1:end);

        MyConst = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];
        M = 16; % Size of the QAM constellation
        constellation = qammod(0:15, M, MyConst);
        stepsize = 2e-4;
        numtaps = 8; % bigger filter works better here
        referencetap = 2;%(numtaps-1)/2;

        
        EQ = comm.LinearEqualizer('Algorithm', 'LMS', 'StepSize', stepsize,'NumTaps', numtaps, 'InputSamplesPerSymbol', 8, 'ReferenceTap', referencetap, 'Constellation', constellation); % 'Constellation', constellation,
        [X_matched,errX] = EQ(X_CD_rec_norm, SIG.Xpol.txSymb);
        [Y_matched,errY] = EQ(Y_CD_rec, SIG.Ypol.txSymb);

    end
    
    transient_Xpol = abs(finddelay(X_matched(1:65536), SIG.Xpol.txSymb));
    X_matched = X_matched(transient_Xpol+1:end);

    if (index==length(OSNR_dB))
        scatterplot(X_matched);
        title(sprintf('%s constellation of Xpol after CMA %d dB of WGN',MODULATIONS(r), OSNR_dB(index)));
    end
   %%
    %------------------Delay&Phase recovery ---------------------
    
    if r==1
        carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r),"SamplesPerSymbol", 1, 'DampingFactor', 150);
        [X_eq, phEstX] = carrSynch(X_matched);
        [Y_eq, phEstY] = carrSynch(Y_matched);
    else
        carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 5, 'NormalizedLoopBandwidth', 5e-3);%, 'ModulationPhaseOffset','Custom', 'CustomPhaseOffset', -pi/7);
        [X_eq, phEstX] = carrSynch(X_matched);
        [Y_eq, phEstY] = carrSynch(Y_matched);
    end

    if (index==length(OSNR_dB))
        scatterplot(X_eq);
        title(sprintf('%s constellation of Xpol after phase recovery',MODULATIONS(r)));
    end    
    
    X_BER = zeros(1,4);
    Y_BER = zeros(1,4);
    j=1;
    
    for i=0:pi/2:3/2*pi
    
    fprintf('---------The phase tried is (degrees): %d-----------\n', (mean(i) *180 /pi));
    
    % fprintf('The total phase recovered is (degrees): %d\n', (mean(phEstX+i) *180 /pi));
    
    transient_Xpol = abs(finddelay(X_eq(1:65536), SIG.Xpol.txSymb));
    transient_Ypol = abs(finddelay(Y_eq(1:65536), SIG.Ypol.txSymb));
    
    fprintf('Transient Xpol: %d\n', transient_Xpol)
    fprintf('Transient Ypol: %d\n', transient_Ypol)
    
    X_RX = X_eq*exp(1i*i);
    X_RX = X_RX(transient_Xpol+1:end);
    Y_RX = Y_eq*exp(1i*i);
    Y_RX = Y_RX(transient_Ypol+1:end);

    if (index==length(OSNR_dB) && j==1)
        scatterplot(X_RX);
        title(sprintf('%s constellation of Xpol after delay recovery',MODULATIONS(r)));
    end
    
    if r==1
        fprintf('The tracked moduluation is: QPSK\n');
          [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);
    %       X_demappedBits = pskdemod(X_RX,M, pi/4*7); %it doesn't demodulate in the same way as our function
    %       X_demappedBits = From_MATLAB_pskdemod(X_demappedBits);
    else
        fprintf('The tracked moduluation is: 16-QAM\n');
       %     X_demappedBits = qamdemod(X_2Sps(1:2:end),M);
     [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_RX,Y_RX);
    end
    
      X_BER(j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
      Y_BER(j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));
    
       j=j+1;
    
    end
    
    X_Ber_Tot_CMA(index) = min(X_BER);
    Y_Ber_Tot_CMA(index) = min(Y_BER);
    fprintf('The BER on Xpol is: %.6f\n', X_Ber_Tot_CMA(index));
    fprintf('The BER on Ypol is: %.6f\n', Y_Ber_Tot_CMA(index));
end


%%
%------------------FIGURES-------------

if r==1
    BER_TH = 0.5 * erfc(sqrt(10.^(OSNR_dB/10)/2));
else
    BER_TH = 3/8 * erfc(sqrt(10.^(OSNR_dB/10)/10));
end

figure();
semilogy(OSNR_dB, X_Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineWidth', 1);
xlim([min(OSNR_dB),15]);
grid on;
hold on;
% semilogy(OSNR_dB, X_Ber_Tot_CMA, 'Marker','o', 'Color', 'b');
semilogy(OSNR_dB, BER_TH, 'r');
title(sprintf('%s BER curve of Xpol',MODULATIONS(r)));
% legend('Simulated BER - Matched filter','Simulated BER - CMA', 'Theoretical BER', 'Interpreter', 'latex');
xlabel('SNR', 'Interpreter','latex');
hold off;

figure();
semilogy(OSNR_dB,Y_Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineWidth', 1);
xlim([min(OSNR_dB),15]);
grid on;
hold on;
% semilogy(OSNR_dB, Y_Ber_Tot_CMA, 'Marker','o', 'Color', 'b');
semilogy(OSNR_dB,BER_TH, 'r');
title(sprintf('%s BER curve of Ypol',MODULATIONS(r)));
% legend('Simulated BER - Matched filter','Simulated BER - CMA', 'Theoretical BER', 'Interpreter', 'latex');
xlabel('SNR', 'Interpreter','latex');
hold off;