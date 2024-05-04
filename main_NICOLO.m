clear;
close all;
clc;

% Load the .mat file
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "16-QAM"];
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
r = 1;
fprintf('The transmitted moduluation is: %s\n', modulation(r));
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_64GBaud.mat'));
if r == 1
    M = 4;
else
    M = 16;
end

%[rx_XPol, rx_YPol] = Matched_filtering(SIG.Xpol.txSig, SIG.Ypol.txSig);

% Create delay and phase convolved signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig);

%add chromatic dispersion
% [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);

% To see CD effect
% figure();
% scatter(real(X_CD), imag(X_CD));

% Adding the noise
[X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_distorted,SIG.Sps, M, 8);
[Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_distorted,SIG.Sps, M, 8);

% %----------------Compensation for CD-------------------
%add cchromatic dispersion
% [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN,SIG.Sps, 2);

% To visually see the recovery
% figure();
% scatter(real(X_CD), imag(X_CD));
%fprintf('isequal = %d\n', isequal(round(X_CD,8), round(SIG.Xpol.txSig, 8)));
%%
%------------------ Matched Flitering ---------------------
X_distorted_AWGN = downsample(X_distorted_AWGN, 4);
Y_distorted_AWGN = downsample(Y_distorted_AWGN, 4);

[X_matched,Y_matched] = Matched_filtering(X_distorted_AWGN, Y_distorted_AWGN, PulseShaping.b_coeff);
%X_delay = finddelay(X_matched(1:65536), SIG.Xpol.txSymb);
%fprintf("%d\n", abs(X_delay));
scatterplot(X_matched(:,1));

%%
%------------------Delay&Phase recovery ---------------------

carrSynch = comm.CarrierSynchronizer("Modulation", modulation,"SamplesPerSymbol", 1);
[X_eq, phEstX] = carrSynch(X_matched);
[Y_eq, phEstY] = carrSynch(Y_matched);

X_avg_amp = mean(abs(X_consolidated))/sqrt(2);
X_rotation = zeros(4,1);
X_delay = zeros(4,1);
Y_avg_amp = mean(abs(X_consolidated))/sqrt(2);
Y_rotation = zeros(4,1);
Y_delay = zeros(4,1);

for idx=1:4
    Xrotated = X_eq.*exp(1i*idx*pi/2);
    X_delay(idx,:) = finddelay(Xrotated, SIG.Xpol.txSymb);
    %Xrotated = circshift(Xrotated,Xdelay(idx));
    Xrotated = Xrotated(Xdelay(idx,:)+1:end);
    X_rotation(idx)=std(Xrotated-X_avg_amp.*SIG.Xpol.txSymb);
    Yrotated = Y_eq.*exp(1i*idx*pi/2);
    Y_delay(idx,:) = finddelay(Yrotated, SIG.Ypol.txSymb);
    %Yrotated = circshift(Yrotated,Ydelay(idx));
    Yrotated = Yrotated(Ydelay(idx,:)+1:end);
    Y_rotation(idx)=std(Yrotated-Y_avg_amp.*SIG.Ypol.txSymb);
end
[M I] = min(X_rotation);
X_consolidated_recovered = circshift(X_consolidated.*exp(1i*I*pi/2), Xdelay(I));
[M I] = min(X_rotation);
Y_consolidated_recovered = circshift(Y_consolidated.*exp(1i*I*pi/2), Ydelay(I));
scatterplot(X_consolidated_recovered)

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

if size(Y_demappedBits) ~= size(SIG.Ypol.bits)
    error('Arrays have different sizes.');
else
    % Calculate the number of bit errors
    bitErrors_Ypol = sum(sum(Y_demappedBits ~= SIG.Ypol.bits));

    % Calculate the total number of bits
    totalBits_Ypol = numel(Y_demappedBits);

    % Calculate the Bit Error Rate (BER)
    BER_Ypol = bitErrors_Ypol / totalBits_Ypol;

    % Display the results
    fprintf('Total number of bits for Y polarization: %d\n', totalBits_Ypol);
    fprintf('Number of bit errors for Y polarization: %d\n', bitErrors_Ypol);
    fprintf('Bit Error Rate (BER) for Y polariztion: %e\n', BER_Ypol);
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
