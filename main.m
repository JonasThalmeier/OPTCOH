clear;
close all;

% Load the .mat file
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "QAM"];
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
r = 1;
fprintf('The transmitted moduluation is: %s\n', modulation(r));
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_64GBaud.mat'));

% Parameters
SpS_down = 4;
SpS_up = SIG.Sps/SpS_down;
OSNR_dB = 3:14;
seq_lenght = length(SIG.Xpol.txSymb);
mu = logspace(-3.4,-3,8);

if r == 1
    M = 4;
else
    M = 16;
end


TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission

% Create delay and phase convolved signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig);

% Adding chromatic dispersion
[X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);

X_Ber_Tot_CMA = zeros(1,length(OSNR_dB));
Y_Ber_Tot_CMA = zeros(1,length(OSNR_dB));
X_Ber_Tot_LMS = zeros(1,length(OSNR_dB));
Y_Ber_Tot_LMS = zeros(1,length(OSNR_dB));

for index = 1:length(OSNR_dB)
    % Adding Noise
    [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
    [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);

    % scatterplot(X_distorted_AWGN,8);
    % title('Noisy constellation, before filtering');

    % Recovering Chromatic Dispersion
    [X_CD_rec,Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, SIG.Sps, 2);

    % Downsampling and Matched Flitering
    X_CD_rec = downsample(X_CD_rec, 4);
    Y_CD_rec = downsample(Y_CD_rec, 4);

    [X_matched,Y_matched] = Matched_filtering(X_CD_rec, Y_CD_rec, PulseShaping.b_coeff);

    % scatterplot(X_matched(1:2:end),1);
    % title('Filtered constellation');
    %%
    % Carrier Synchroniazation and Normalization
    if r==1
        mu = 1e-3;
        NTaps = 9;
        N1 = 1e3;
        N2 = 1e4;
        XY_eq = EQ_func(X_CD_rec,Y_CD_rec,mu,NTaps,"CMA",N1,N2);

        Delta_nu = 50e3; % Laser line width
        Es = 1; % Symbol energy (=radius)
        Npol = 2;
        windowlen = 100;
        XY_vit = vit_n_vit(XY_eq, Delta_nu, SIG.symbolRate, OSNR_dB(index), Es, Npol, M, windowlen);
        X_eq = XY_vit(:,1);
        Y_eq = XY_vit(:,2);
    else
        mu = 1e-3;
        NTaps = 9;
        N1 = 1e3;
        N2 = 1e4;
        XY_eq = EQ_func(X_CD_rec,Y_CD_rec,mu,NTaps,"RDE",N1,N2);
        % scatterplot(XY_eq(1e5:end,1)); title('Post EQ, pre carrSync');
        carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 10, 'NormalizedLoopBandwidth',1e-2);
        [X_eq, phEstX] = carrSynch(XY_eq(:,1));
        [Y_eq, phEstY] = carrSynch(XY_eq(:,2));
    end

    X_Power = mean(abs((X_eq)).^2);
    X_eq = X_eq/sqrt(X_Power/10);
    Y_Power = mean(abs((Y_eq)).^2);
    Y_eq = Y_eq/sqrt(Y_Power/10);
    X_eq = X_eq(length(SIG.Xpol.bits):end);
    Y_eq = Y_eq(length(SIG.Ypol.bits):end);
    % scatterplot(X_eq(1e5:end));
    % title('Constellation after Carrier Synchronization');
    %%
    % Finding the constellations right orientation
    X_BER = zeros(4,4);
    Y_BER = zeros(4,4);
    j=1;
    k=1;
    for rot=0:pi/2:3/2*pi
    for i=0:pi/2:3/2*pi
        fprintf('---------The phase tried is (degrees): %d-----------\n', (mean(i) *180 /pi));

        transient_Xpol = abs(finddelay(X_eq(1:seq_lenght), SIG.Xpol.txSymb));
        transient_Ypol = abs(finddelay(Y_eq(1:seq_lenght), SIG.Ypol.txSymb));

        fprintf('Transient Xpol: %d\n', transient_Xpol)
        fprintf('Transient Ypol: %d\n', transient_Ypol)

        N = 5;
        X_RX = X_eq*exp(1i*i);
        X_RX = X_RX((N-1)*seq_lenght+transient_Xpol+1:N*seq_lenght+transient_Xpol+transient_Xpol);
        Y_RX = Y_eq*exp(1i*i);
        Y_RX = Y_RX((N-1)*seq_lenght+transient_Ypol+1:N*seq_lenght+transient_Xpol+transient_Ypol);
        
        X_RX = X_RX.*cos(rot)+Y_RX.*sin(rot);
        Y_RX = -X_RX.*sin(rot)+Y_RX.*cos(rot);

        if r==1
            fprintf('The tracked moduluation is: QPSK\n');
            [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);
        else
            fprintf('The tracked moduluation is: 16-QAM\n');
            %         MyConst = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];
            %         X_demappedBits = qamdemod(X_RX, M, MyConst, OutputType='bit', PlotConstellation=true);
            %         N = length(X_demappedBits)/4;
            %         X_demappedBits = reshape(X_demappedBits, 4, N).';
            %
            [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_RX,Y_RX);
        end
        X_BER(k,j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
        Y_BER(k,j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));
        j=j+1;
        k=k+1;
    end
    end
    X_Ber_Tot(index) = min(X_BER,[],'all');
    Y_Ber_Tot(index) = min(Y_BER,[],'all');
end
%%
if r == 1
    BER_TH = 0.5 * erfc(sqrt(10.^(OSNR_dB/10)/2));
    BER_MED_MF = 0.5 * (X_Ber_Tot + Y_Ber_Tot);
    figure();
    semilogy(OSNR_dB, BER_TH, 'r', 'LineWidth', 1);
    xlim([1,14]);
    grid on;
    hold on;
    semilogy(OSNR_dB, BER_MED_MF, 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
    title(sprintf('%s BER curve',MODULATIONS(r)));
    legend('Theoretical BER','Simulated BER - CMA', 'Interpreter', 'latex');
    xlabel('OSNR [dB]', 'Interpreter','latex');
    hold off;
    fprintf('The BER on Xpol is: %.6f\n', X_Ber_Tot);
    fprintf('The BER on Ypol is: %.6f\n', Y_Ber_Tot);
else
    BER_TH = 3/8 * erfc(sqrt(10.^(OSNR_dB/10)/10));
    BER_MED_MF = 0.5 * (X_Ber_Tot + Y_Ber_Tot);
    figure();
    semilogy(OSNR_dB, BER_TH, 'r', 'LineWidth', 1);
    % xlim([1,14]);
    grid on;
    hold on;
    semilogy(OSNR_dB, BER_MED_MF, 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
    title(sprintf('%s BER curve',MODULATIONS(r)));
    legend('Theoretical BER','Simulated BER - CMA', 'Interpreter', 'latex');
    xlabel('OSNR [dB]', 'Interpreter','latex');
    hold off;
    fprintf('The BER on Xpol is: %.6f\n', X_Ber_Tot);
    fprintf('The BER on Ypol is: %.6f\n', Y_Ber_Tot);
end
