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

figure();
scatter(real(SIG.Xpol.txSig), imag(SIG.Xpol.txSig));

%add chromatic dispersion
[X_CD,Y_CD]=Chromatic_Dispersion(SIG.Xpol.txSig, SIG.Ypol.txSig, SIG.Sps, 1);

figure();
scatter(real(X_CD), imag(X_CD));



% Create delay and phase convolved signals
[X_distorted, Y_distorted] = DP_Distortion(X_CD,Y_CD);

% Adding the noise
[X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_distorted,SIG.Sps, M, 2);
[Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_distorted,SIG.Sps, M, 2);

% %----------------Compensation for CD-------------------
% %add cchromatic dispersion
% [X_CD,Y_CD]=Chromatic_Dispersion(X_CD, Y_CD,SIG.Sps, 2);
% 
% figure();
% scatter(real(X_CD), imag(X_CD));
% fprintf('isequal = %d\n', isequal(round(X_CD,8), round(SIG.Xpol.txSig, 8)));

%----------------Pulse shaping and Downsample the signal-------------------
X_distorted_AWGN = conv(PulseShaping.b_coeff, X_CD);
Y_distorted_AWGN = conv(PulseShaping.b_coeff, Y_CD);

%----------------Compensation for CD-------------------
%add cchromatic dispersion
[X_CD,Y_CD]=Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN,SIG.Sps, 2);

figure();
scatter(real(Y_CD), imag(Y_CD));
%%

%---------------------------EQ---------------------------------------------
%From 15 numtaps are the best for qpsk, also for qam, but 60 and 120 seem
%better. The reference tap for qpsk are 1,3,65,121.
if r == 1
    constellation = pskmod(0:3, 4, pi/4);
    stepsize = 1e-3; %best result
    numtaps = 9;
    referencetap = (numtaps-1)/2;
    modulation = 'QPSK';
else
    constellation = qammod(0:15, 16);
    stepsize = 5e-6;
    numtaps = 9; % bigger filter works better here
    referencetap = (numtaps-1)/2;
    modulation = 'QAM';
end

aw = true;
X_dist_2SpS = X_CD(1:4:end);
Y_dist_2SpS = Y_CD(1:4:end);

EQ = comm.LinearEqualizer('Algorithm', 'CMA', 'StepSize', stepsize,'NumTaps', numtaps, 'InputSamplesPerSymbol', 2, 'Constellation', constellation, 'ReferenceTap', referencetap, 'InputDelay', 0);
[X_eq,errX] = EQ(X_dist_2SpS(1:end-mod(length(X_dist_2SpS),2)));
[Y_eq,errY] = EQ(Y_dist_2SpS(1:end-mod(length(Y_dist_2SpS),2)));
% figure; plot(abs(err));
% constell = comm.ConstellationDiagram('NumInputPorts', 1, 'SamplesPerSymbol', SpS_up, 'ReferenceConstellation', constellation, 'Title', 'Before phase correction');
% constell(X_eq);
% scatterplot(X_eq);

carrSynch = comm.CarrierSynchronizer("Modulation", modulation,"SamplesPerSymbol", 1);
[X_eq, phEstX] = carrSynch(X_eq);
[Y_eq, phEstY] = carrSynch(Y_eq);

%-----------------Consolidation--------------------------------------------
N = length(SIG.Xpol.txSymb);  % Calculate the number of unique symbols
X_consolidated = zeros(N, 1);
Y_consolidated = zeros(N, 1);
for idx = 1:N
    X_consolidated(idx) = mean(X_eq(idx:N:end,:));
    Y_consolidated(idx) = mean(Y_eq(idx:N:end,:));
end

%-----------find right rotation and delay----------------------------------
X_avg_amp = mean(abs(X_consolidated))/sqrt(2);
X_rotation = zeros(4,1);
X_delay = zeros(4,1);
Y_avg_amp = mean(abs(X_consolidated))/sqrt(2);
Y_rotation = zeros(4,1);
Y_delay = zeros(4,1);
for idx=1:4
    Xrotated = X_consolidated.*exp(1i*idx*pi/2);
    Xdelay(idx) = finddelay(Xrotated, SIG.Xpol.txSymb);
    Xrotated = circshift(Xrotated,Xdelay(idx));
    X_rotation(idx)=std(Xrotated-X_avg_amp.*SIG.Xpol.txSymb);
    Yrotated = Y_consolidated.*exp(1i*idx*pi/2);
    Ydelay(idx) = finddelay(Yrotated, SIG.Ypol.txSymb);
    Yrotated = circshift(Yrotated,Ydelay(idx));
    Y_rotation(idx)=std(Yrotated-Y_avg_amp.*SIG.Ypol.txSymb);
end
[M I] = min(X_rotation);
X_consolidated_recovered = circshift(X_consolidated.*exp(1i*I*pi/2), Xdelay(I));
[M I] = min(X_rotation);
Y_consolidated_recovered = circshift(Y_consolidated.*exp(1i*I*pi/2), Ydelay(I));
scatterplot(X_consolidated_recovered)

if r == 1
    M = 2;
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_consolidated_recovered,Y_consolidated_recovered);
else
    M = 4;
    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_consolidated_recovered,Y_consolidated_recovered);
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