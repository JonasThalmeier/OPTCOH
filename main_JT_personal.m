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
OSNR_dB = 6:4:18;
seq_lenght = length(SIG.Xpol.txSymb);
B = 30;
N = 50;
sweep_par = 1:3;


if r == 1
    M = 4;
    power_norm = 2;
else
    M = 16;
    power_norm = 10;
end

BERmat = zeros(length(sweep_par), length(OSNR_dB));
for idx_sweep = 1:length(sweep_par)
    for idx_OSNR = 1:length(OSNR_dB)
        TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
        TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission

        % Create delay and phase convolved signals
        % [X_distorted, Y_distorted] = DP_Distortion_N(SIG.Xpol.txSig, SIG.Ypol.txSig);
        % halfleng = round(1*length(SIG.Xpol.txSig));
        [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, sweep_par(idx_sweep), 1e3);

        % seq_lenght = length(SIG.Xpol.txSig(1:halfleng));

        % Adding chromatic dispersion
        [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);

        X_Ber_Tot_CMA = zeros(1,length(OSNR_dB));
        Y_Ber_Tot_CMA = zeros(1,length(OSNR_dB));
        X_Ber_Tot_LMS = zeros(1,length(OSNR_dB));
        Y_Ber_Tot_LMS = zeros(1,length(OSNR_dB));


        % Adding Noise
        [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB(idx_OSNR), SIG.symbolRate);
        [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB(idx_OSNR), SIG.symbolRate);

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
            XY_vit = vit_n_vit(XY_eq, Delta_nu, SIG.symbolRate, OSNR_dB(idx_OSNR), Es, Npol, M, windowlen);
            X_eq = XY_vit(:,1);
            Y_eq = XY_vit(:,2);
        else
            mu = 1e-3;
            mu2 = 5e-5;
            % mu = 5e-3;
            NTaps = 9;
            N1 = 5e4;
            N2 = 5e3;
            if sweep_par(idx_sweep) == 1
                XY_eq = EQ_func(X_CD_rec,Y_CD_rec,mu,NTaps,"RDE",N1,N2);
            elseif sweep_par(idx_sweep) == 2
                mu2 = 5e-4;
                [X_eq, Y_eq, e_X, e_Y] = EQ_func_N([X_CD_rec,Y_CD_rec],r,mu,mu2,NTaps,N1,N2);
                XY_eq = [X_eq, Y_eq];
            else
                [X_eq, Y_eq, e_X, e_Y] = EQ_func_J([X_CD_rec,Y_CD_rec],r,mu,mu2,NTaps,N1,N2);
                XY_eq = [X_eq, Y_eq];
            end
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
        BERmat(idx_sweep,idx_OSNR) = .5*(min(X_BER,[],'all')+min(Y_BER,[],'all'));
        X_Ber_Tot(idx_OSNR) = min(X_BER,[],'all');
        Y_Ber_Tot(idx_OSNR) = min(Y_BER,[],'all');
    end
end
save('matrixBER.mat', 'BERmat');
%%
% figure;
% h1 = heatmap(BERmat);
% h1.Title = 'BER at 8dB';
% h1.XLabel = 'N';
% h1.YLabel = 'B';
% h1.Colormap = parula;  % Ändert das Farbschema
% h1.XDisplayLabels = Nvec;
% h1.YDisplayLabels = Bvec;
% figure;
% h2 = heatmap(log10(squeeze(BERmat(1,:,:))));
% h2.Title = 'log10(BER)';
% h2.YLabel = 'N';
% h2.XLabel = 'SNR';
% h2.Colormap = parula;  % Ändert das Farbschema
% h2.YDisplayLabels = Nvec;
% h2.XDisplayLabels = OSNR_dB;

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
    legend('Theoretical BER','Simulated BER - CMA', 'X Pol', 'Y Pol', 'Interpreter', 'latex');
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
    semilogy(OSNR_dB, BERmat(1,:), 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
    semilogy(OSNR_dB, BERmat(2,:), 'Marker','o', 'LineStyle','-.', 'LineWidth', 1);
    semilogy(OSNR_dB, BERmat(3,:), 'Marker','o', 'LineStyle','-.', 'LineWidth', 1);
    % title(sprintf('%s BER curve (B=%i, N=%i)',MODULATIONS(r),Bvec,Nvec));
    legend('Theoretical BER','Simulated BER - CMA Book', 'Simulated BER - CMA Nicolo', 'Simulated BER - CMA Jonas', 'Interpreter', 'latex');
    xlabel('OSNR [dB]', 'Interpreter','latex');
    hold off;
    fprintf('The BER on Xpol is: %.6f\n', X_Ber_Tot);
    fprintf('The BER on Ypol is: %.6f\n', Y_Ber_Tot);


    % figure();
    % hold on
    % for idx_OSNR = 1:length(OSNR_dB)
    %         SNR = 10*log10(10*erfinv(-8/3.*squeeze(BERmat(:,idx_OSNR))+1).^2);
    %         plot(log10(sweep_par), OSNR_dB(idx_OSNR)-SNR, 'LineWidth', 1, 'DisplayName', sprintf('OSNR = %d dB', OSNR_dB(idx_OSNR)));
    % end
    % % xlim([1,14]);
    % grid on;
    % title(sprintf('%s SNR penalty',MODULATIONS(r)));
    % legend show;
    % xlabel('$log_{10}\Delta\nu$', 'Interpreter','latex');
    % ylabel('$\Delta$OSNR [dB]', 'Interpreter','latex');
    % hold off;
    % fprintf('The BER on Xpol is: %.6f\n', X_Ber_Tot);
    % fprintf('The BER on Ypol is: %.6f\n', Y_Ber_Tot);
end
runtime_out = toc(runtime)