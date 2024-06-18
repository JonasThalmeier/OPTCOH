clear;
close all;
runtime = tic;
% Load the .mat file
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "QAM"];
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
r = 2;
fprintf('The transmitted moduluation is: %s\n', modulation(r));
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_64GBaud.mat'));

% Parameters
SpS_down = 4;
SpS_up = SIG.Sps/SpS_down;
seq_lenght = length(SIG.Xpol.txSymb);
B = 30;
N = 50;
sweep_par = linspace(0,1e5,4);
BER_goal = 1e-3;


if r == 1
    M = 4;
    power_norm = 2;
    SNR_opt = 2*log10(10*erfinv(1-2*BER_goal)^2);
else
    M = 16;
    power_norm = 10;
    SNR_opt = 10*log10(10*erfinv(1-8/3*BER_goal)^2);
end


SNR = SNR_opt;
BER = ones(length(sweep_par),2);
delta_SNR = zeros(length(sweep_par),1);
for idx_sweep = 1:length(sweep_par)
    SNR = SNR_opt;
    while BER(idx_sweep,2)>BER_goal
        TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
        TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission

        % Create delay and phase convolved signals
        % [X_distorted, Y_distorted] = DP_Distortion_N(SIG.Xpol.txSig, SIG.Ypol.txSig);
        % halfleng = round(1*length(SIG.Xpol.txSig));
        [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, 50e3, sweep_par(idx_sweep));

        % seq_lenght = length(SIG.Xpol.txSig(1:halfleng));

        % Adding chromatic dispersion
        [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);


        % Adding Noise
        [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, SNR, SIG.symbolRate);
        [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, SNR, SIG.symbolRate);

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
            Rs = 64e9;
            Es = 1; % Symbol energy (=radius)
            Npol = 2;
            windowlen = 100;
            XY_vit = vit_n_vit(XY_eq, Delta_nu, SIG.symbolRate, SNR, Es, Npol, M, windowlen);
            X_eq = XY_vit(:,1);
            Y_eq = XY_vit(:,2);
        else
            mu = 1e-3;
            mu2 = 1e-5;
            NTaps = 9;
            N1 = 5e4;
            N2 = 5e3;
            [X_eq, Y_eq, e_X, e_Y] = EQ_func_J([X_CD_rec,Y_CD_rec],r,mu,mu2,NTaps,N1,N2);
            XY_eq = [X_eq, Y_eq];
            carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 100, 'NormalizedLoopBandwidth',1e-3);
            [X_eq, phEstX] = carrSynch(XY_eq(:,1));
            [Y_eq, phEstY] = carrSynch(XY_eq(:,2));
            X_eq = XY_eq(:,1).*exp(-1i*movmean(phEstX,1e3));
            Y_eq = XY_eq(:,2).*exp(-1i*movmean(phEstY,1e3));
            % [X_eq,~] = BPS_N(XY_eq(:,1),100,10,M,power_norm);
            % [Y_eq,~] = BPS_N(XY_eq(:,2),100,10,M,power_norm);
        end


        % [X_eq,~] = BPS_N(XY_eq(:,1),Bvec(idxB),Nvec(idxN),M,power_norm);
        % [Y_eq,~] = BPS_N(XY_eq(:,2),Bvec(idxB),Nvec(idxN),M,power_norm);

        X_Power = mean(abs((X_eq)).^2);
        X_eq = X_eq/sqrt(X_Power/10);
        Y_Power = mean(abs((Y_eq)).^2);
        Y_eq = Y_eq/sqrt(Y_Power/10);
        % X_eq = X_eq(2*length(SIG.Xpol.bits):end);
        % Y_eq = Y_eq(2*length(SIG.Ypol.bits):end);
        % scatterplot(X_eq(1e5:end));
        % title('Constellation after Carrier Synchronization');
        %%
        % Finding the constellations right orientation
        X_BER = zeros(1,4);
        Y_BER = zeros(1,4);
        k=1;
        for rot=0:pi/2:3/2*pi
            j=1;
            for i=0:pi/2:3/2*pi
                % fprintf('---------The phase tried is (degrees): %d-----------\n', (mean(i) *180 /pi));



                % fprintf('Transient Xpol: %d\n', transient_Xpol)
                % fprintf('Transient Ypol: %d\n', transient_Ypol)

                N = 5;
                X_RX = X_eq*exp(1i*i);
                Y_RX = Y_eq*exp(1i*i);
                X_RX = X_RX.*cos(rot)-Y_RX.*sin(rot);
                Y_RX = X_RX.*sin(rot)+Y_RX.*cos(rot);
                transient_Xpol = abs(finddelay(X_eq(1:seq_lenght), SIG.Xpol.txSymb));
                transient_Ypol = abs(finddelay(Y_eq(1:seq_lenght), SIG.Ypol.txSymb));
                X_RX = X_RX((N-1)*seq_lenght+transient_Xpol+1:N*seq_lenght+transient_Xpol+transient_Xpol);
                Y_RX = Y_RX((N-1)*seq_lenght+transient_Ypol+1:N*seq_lenght+transient_Xpol+transient_Ypol);


                if r==1
                    % fprintf('The tracked moduluation is: QPSK\n');
                    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);
                else
                    % fprintf('The tracked moduluation is: 16-QAM\n');
                    %         MyConst = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];
                    %         X_demappedBits = qamdemod(X_RX, M, MyConst, OutputType='bit', PlotConstellation=true);
                    %         N = length(X_demappedBits)/4;
                    %         X_demappedBits = reshape(X_demappedBits, 4, N).';
                    %
                    [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_RX,Y_RX);
                end
                X_BER(k,j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
                Y_BER(k,j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));
                YX_BER(k,j) = biterr(X_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));
                XY_BER(k,j) = biterr(Y_demappedBits, TX_BITS_Xpol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));
                j=j+1;
            end
            k=k+1;
        end
        X_Ber = min(X_BER,[],'all');
        Y_Ber = min(Y_BER,[],'all');
        BER(idx_sweep,:) = [BER(idx_sweep,2:end), .5*(X_Ber+Y_Ber)];
        SNR = SNR+.5;
    end
    SNR_coarse = [SNR-1, SNR-.5];
    SNR_fine = linspace(SNR-1,SNR-.5,1e4);
    delta_BER = interp1(SNR_coarse,BER(idx_sweep,:),SNR_fine)-BER_goal;
    [~,min_idx] = min(abs(delta_BER));
    delta_SNR(idx_sweep) = SNR_fine(min_idx)-SNR_opt;
end
%%



figure();
plot(sweep_par, delta_SNR, 'k', 'LineWidth', 2);
grid on;
hold on;
title(sprintf('SNR penalty (%s)',MODULATIONS(r)));
% legend('Theoretical BER','Simulated BER - CMA', 'X Pol', 'Y Pol', 'Interpreter', 'latex');
xlabel('Sweep Parameter', 'Interpreter','latex');
ylabel('$\Delta SNR$', 'Interpreter','latex');
hold off;


runtime_out = toc(runtime)