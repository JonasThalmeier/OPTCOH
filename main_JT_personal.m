clear;
close all;
clc;

% Load the .mat file
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "16-QAM"];
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
r = 2;
fprintf('The transmitted moduluation is: %s\n', modulation(r));
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_64GBaud.mat'));
if r == 1
    M = 4;
else
    M = 16;
end

TX_BITS_Xpol = [SIG.Xpol.bits; SIG.Xpol.bits; SIG.Xpol.bits; SIG.Xpol.bits; SIG.Xpol.bits; SIG.Xpol.bits; SIG.Xpol.bits; SIG.Xpol.bits; SIG.Xpol.bits; SIG.Xpol.bits];

%[rx_XPol, rx_YPol] = Matched_filtering(SIG.Xpol.txSig, SIG.Ypol.txSig);

% Create delay and phase convolved signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig);

%add chromatic dispersion
% [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);

% To see CD effect
% figure();
% scatter(real(X_CD), imag(X_CD));

% Adding the noise
OSNR_dB = 1:20;
X_Ber_Tot = zeros(1,length(OSNR_dB));

% X_distorted = SIG.Xpol.txSig;
% Y_distorted = SIG.Ypol.txSig;

for OSNR_dB_i = 1: length(OSNR_dB)

    OSNR_dB_i_b = OSNR_dB(OSNR_dB_i);
    [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_distorted, SIG.Sps, M, OSNR_dB_i_b, SIG.symbolRate);
    [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_distorted, SIG.Sps, M, OSNR_dB_i_b, SIG.symbolRate);

    % %----------------Compensation for CD-------------------
    %add cchromatic dispersion
    % [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN,SIG.Sps, 2);

    % To visually see the recovery
    % figure();
    % scatter(real(X_CD), imag(X_CD));
    %fprintf('isequal = %d\n', isequal(round(X_CD,8), round(SIG.Xpol.txSig, 8)));

    %------------------ Matched Flitering ---------------------
    X_distorted_AWGN = downsample(X_distorted_AWGN, 4);
    Y_distorted_AWGN = downsample(Y_distorted_AWGN, 4);

    [X_matched,Y_matched] = Matched_filtering(X_distorted_AWGN, Y_distorted_AWGN, PulseShaping.b_coeff);

    % scatterplot(X_matched(1:2:end));
    % if (OSNR_dB_i==20)
    %     scatterplot(X_matched(1:2:end));
    % pause;
    % end
    %------------------Delay&Phase recovery ---------------------
    %
    [X_matched,Y_matched] = samp_time_recovery(X_matched,Y_matched,8);
    carrSynch = comm.CarrierSynchronizer("Modulation", "QAM","SamplesPerSymbol", 1,'DampingFactor', 5, 'NormalizedLoopBandwidth', 5e-3);%, 'ModulationPhaseOffset','Custom', 'CustomPhaseOffset', -pi/7);
    [X_eq, phEstX] = carrSynch(5*X_matched(1:2:end));
    [Y_eq, phEstY] = carrSynch(5*Y_matched(1:2:end));
    % scatterplot(X_eq(1:end));
    % scatterplot(X_matched(1:2:end));

    
    % X_eq = X_matched(1:2:end); % era X_eq
    % X_eq = X_eq(1:2:end); % era X_eq

    X_BER = zeros(1,4);
    j=1;

    for i=0:pi/2:3/2*pi

        fprintf('---------The phase tried is (degrees): %d-----------\n', (mean(i) *180 /pi));

        % fprintf('The total phase recovered is (degrees): %d\n', (mean(phEstX+i) *180 /pi));

        transient_Xpol = abs(finddelay(X_eq(1:65536), SIG.Xpol.txSymb));
        fprintf('%d\n', transient_Xpol)
        X_RX = X_eq*exp(1i*i);
        X_RX = X_RX(transient_Xpol+1:end);
        % constell3 = comm.ConstellationDiagram('NumInputPorts', 1, 'SamplesPerSymbol', 1, 'ReferenceConstellation', constellation, 'Title', 'After transient correction');
        % constell3(X_eq);
        %scatterplot(X_RX);

        Y_2Sps =X_eq; %Just to give it a value, for the moment i test only the X_pol

        if r==1
            fprintf('The tracked moduluation is: QPSK\n');
            [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX,Y_2Sps);
            %       X_demappedBits = pskdemod(X_RX,M, pi/4*7); %it doesn't demodulate in the same way as our function
            %       X_demappedBits = From_MATLAB_pskdemod(X_demappedBits);
        else
            fprintf('The tracked moduluation is: 16-QAM\n');
            %     X_demappedBits = qamdemod(X_2Sps(1:2:end),M);
            [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_RX,Y_2Sps);
        end

        X_BER(j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
        j=j+1;

    end

    X_Ber_Tot(OSNR_dB_i) = min(X_BER);
    fprintf('The BER on Xpol is: %.26f\n', X_Ber_Tot(OSNR_dB_i));
end

% BER_TH = 0.5 * erfc(sqrt(10.^(OSNR_dB/10)/2));
BER_TH = 3/8 * erfc(sqrt(2*10.^(OSNR_dB/10)/5));
BER_TH2 = 3/8 * erfc(sqrt(10.^(OSNR_dB/10)/10));

figure();
semilogy(OSNR_dB,X_Ber_Tot, 'Marker','o');
% xlim([2,15]);
grid on;
hold on;
% semilogy(OSNR_dB,BER_TH);
semilogy(OSNR_dB,BER_TH2);
legend('sim', 'analytic');
















%%
%-------------------Demapping------------

if r == 1
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_eq,Y_eq);
else
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_eq,Y_eq);
end

%-------------------------BER calculation----------------------------------

if size(X_demappedBits) ~= size(SIG.Xpol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Xpol = sum(sum(X_demappedBits ~= SIG.Xpol.bits));

    % Calculate the total number of bits
    totalBits_Xpol = numel(X_demappedBits);

    % Calculate the Bit Error Rate (BER)
    BER_Xpol = bitErrors_Xpol / totalBits_Xpol;

    % Display the results
    fprintf('Total number of bits for X polarization: %d\n', totalBits_Xpol);
    fprintf('Number of bit errors for X polarization: %d\n', bitErrors_Xpol);
    fprintf('Bit Error Rate (BER) for X polariztion: %e\n', BER_Xpol);
end


%%
%---------------------------EQ---------------------------------------------
%From 15 numtaps are the best for qpsk, also for qam, but 60 and 120 seem
%better. The reference tap for qpsk are 1,3,65,121.
if r == 1
    constellation = pskmod(0:3, 4, pi/4);
    stepsize = 0.0001; %best result
    numtaps = 31;
    referencetap = 15;
    modulation = 'QPSK';
else
    constellation = qammod(0:15, 16);
    stepsize = 0.00001;
    numtaps = 31; % bigger filter works better here
    referencetap = 15;
    modulation = 'QAM';
end

aw = true;

EQ = comm.LinearEqualizer('Algorithm', 'CMA', 'StepSize', stepsize,'NumTaps', numtaps, 'InputSamplesPerSymbol', 2, 'Constellation', constellation, 'ReferenceTap', referencetap, 'InputDelay', abs(rx_Xpol_Delay));
[X_eq,err] = EQ(X_distorted_AGWN(1:end-mod(length(X_distorted_AGWN),2))); % here X_eq is already at 1SpS
% constell = comm.ConstellationDiagram('NumInputPorts', 1, 'SamplesPerSymbol', SpS_up, 'ReferenceConstellation', constellation, 'Title', 'Before phase correction');
% constell(X_eq);
scatterplot(X_eq);

carrSynch = comm.CarrierSynchronizer("Modulation", modulation,"SamplesPerSymbol", 1);

[X_eq, phEst] = carrSynch(X_eq);
fprintf('The random phase recovered is (degrees): %d\n', (mean(phEst) *180 /pi));
% constell2 = comm.ConstellationDiagram('NumInputPorts', 1, 'SamplesPerSymbol', 1, 'ReferenceConstellation', constellation, 'Title', 'After phase correction');
% constell2(X_eq);
scatterplot(X_eq);

%Majority voting

plot(abs(err))
xlabel('Symbols')
ylabel('Error Magnitude')
grid on
title('Time-Varying Channel Without Retraining')

for i=0:pi/2:3/2*pi

    fprintf('---------The phase tried is (degrees): %d-----------\n', (mean(i) *180 /pi));

    fprintf('The total phase recovered is (degrees): %d\n', (mean(phEst+i) *180 /pi));

    X_eq = X_eq*exp(1i*i);
    transient_Xpol = abs(finddelay(X_eq(1:65536), SIG.Xpol.txSymb));
    X_eq = X_eq(transient_Xpol+1:end-transient_Xpol);
    % constell3 = comm.ConstellationDiagram('NumInputPorts', 1, 'SamplesPerSymbol', 1, 'ReferenceConstellation', constellation, 'Title', 'After transient correction');
    % constell3(X_eq);
    scatterplot(X_eq);

    Y_2Sps =X_eq; %Just to give it a value, for the moment i test only the X_pol

    [counts, binEdges] = histcounts(angle(X_eq(1:2:end)), 12, 'Normalization', 'probability');
    if max(counts) > .17
        fprintf('The tracked moduluation is: QPSK\n');
        M = 4;
        [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_eq,Y_2Sps);
        %     X_demappedBits_2 = pskdemod(X_2Sps(1:2:end),M, pi/4); it doesn't demodulate in the same way as our function
    else
        fprintf('The tracked moduluation is: 16-QAM\n');
        M = 16;
        %     X_demappedBits = qamdemod(X_2Sps(1:2:end),M);
        [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_eq,Y_2Sps);
    end


    X_BER = sum(sum(X_demappedBits(1:length(SIG.Xpol.bits),:) ~= SIG.Xpol.bits))/(length(SIG.Xpol.bits)*4);
    fprintf('The BER on Xpol is: %.26f\n', X_BER);

end
% scatterplot(X_eq);
% eyediagram(X_eq,2*SpS_down);


%%

%-------------Remove Transient at the end of transmission-----------------
X_2Sps = X_2Sps(1:end-(length(PulseShaping.b_coeff)/(SIG.Sps)));
Y_2Sps = Y_2Sps(1:end-(length(PulseShaping.b_coeff)/(SIG.Sps)));


% Plot constellation
figure;
scatter(real(X_2Sps(1:2:end)), imag(X_2Sps(1:2:end)), ".", "k");
title('Xpol constellation');
grid on;

figure;
scatter(real(Y_2Sps(1:2:end)), imag(Y_2Sps(1:2:end)), ".", "k");
title('Ypol constellation');
grid on;

[counts, binEdges] = histcounts(angle(X_2Sps(1:2:end)), 12, 'Normalization', 'probability');
if max(counts) > .17
    fprintf('The tracked moduluation is: QPSK\n');
    M = 2;
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_2Sps,Y_2Sps);
else
    fprintf('The tracked moduluation is: 16-QAM\n');
    M = 4;
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_2Sps,Y_2Sps);
end

%----------------Majority voting over the repeated symbols-----------------
% Let's assume demappedBits_Xpol is your bit outcomes with size [N * Npp, 2]
% Where N is the number of unique symbols and Npp is the number of repetitions

N = length(SIG.Xpol.txSymb);  % Calculate the number of unique symbols

X_consolidated = zeros(N, size(SIG.Xpol.bits,2));
Y_consolidatedBits = zeros(N, size(SIG.Xpol.bits,2));

for i = 1:N
    X_consolidated(i, :) = mode(X_demappedBits(i:N:end,:),1);
    Y_consolidatedBits(i, :) = mode(Y_demappedBits(i:N:end,:),1);
end
% Now, consolidatedBits contains the 'averaged' bit decisions





% % ----- Recover from delay and phase -------------------------------------
% X_2Sps = Recover_Delay_Phase_Noise(X_2Sps,SIG.Xpol.txSymb);
% Y_2Sps = Recover_Delay_Phase_Noise(Y_2Sps,SIG.Ypol.txSymb);
%-------------------------BER calculation----------------------------------

if size(X_consolidated) ~= size(SIG.Xpol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Xpol = sum(sum(X_consolidated ~= SIG.Xpol.bits));

    % Calculate the total number of bits
    totalBits_Xpol = numel(X_consolidated);

    % Calculate the Bit Error Rate (BER)
    BER_Xpol = bitErrors_Xpol / totalBits_Xpol;

    % Display the results
    fprintf('Total number of bits for X polarization: %d\n', totalBits_Xpol);
    fprintf('Number of bit errors for X polarization: %d\n', bitErrors_Xpol);
    fprintf('Bit Error Rate (BER) for X polariztion: %f\n', BER_Xpol);
end

if size(Y_consolidatedBits) ~= size(SIG.Ypol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Ypol = sum(sum(Y_consolidatedBits ~= SIG.Ypol.bits));

    % Calculate the total number of bits
    totalBits_Ypol = numel(Y_consolidatedBits);

    % Calculate the Bit Error Rate (BER)
    BER_Ypol = bitErrors_Ypol / totalBits_Ypol;

    % Display the results
    fprintf('Total number of bits for Y polarization: %d\n', totalBits_Ypol);
    fprintf('Number of bit errors for Y polarization: %d\n', bitErrors_Ypol);
    fprintf('Bit Error Rate (BER) for Y polariztion: %f\n', BER_Ypol);
end

function y = phaseCorrection(y)
a = angle(y((real(y) > 0) & (imag(y) > 0)));
a(a < 0.1) = a(a < 0.1) + pi/2;
theta = mean(a) - pi/4;
y = y * exp(-1i*theta);
fprintf('The random phase recovered is (degrees): %d\n', (theta *180 /pi));
end
