clear;
close all;
clc;

% Load the .mat file
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "QAM"];
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
r = 1;
fprintf('The transmitted moduluation is: %s\n', modulation(r));
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_64GBaud.mat'));
M=4
OSNR_dB=10;
%add noise
[X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(SIG.Xpol.txSig, SIG.Sps, M, OSNR_dB, SIG.symbolRate);
[Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(SIG.Ypol.txSig, SIG.Sps, M, OSNR_dB, SIG.symbolRate);
%go to 2 SpS --> needed for equalization
X_distorted_AWGN=downsample(X_distorted_AWGN, 4);
Y_distorted_AWGN=downsample(Y_distorted_AWGN, 4);
%needed to normalize the power to 1 
X_Power = mean(abs((X_distorted_AWGN)).^2);
X_eq = X_distorted_AWGN/sqrt(X_Power);
Y_Power = mean(abs((Y_distorted_AWGN)).^2);
Y_eq = Y_distorted_AWGN/sqrt(Y_Power);
%consider first repetition as training
X_eq_training=X_eq(1:65536*2); 
Y_eq_training=Y_eq(1:65536*2);

mu=0.01

for N_taps = 5:1:12
    test=EQ_func_1(X_eq,Y_eq,mu,N_taps,"CMA",X_eq_training,Y_eq_training);
    scatterplot(test(10000:end,1)); title(sprintf('Output EQ with mu=10^{%d}',N_taps));
end
%scatterplot(dist1(round(N/2):end)); title("EQ input");