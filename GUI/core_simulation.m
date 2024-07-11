function [Ber_Tot] = core_simulation(X_CD,Y_CD,r,Rs, OSNR_dB, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, EQ_N2, CarSync_DampFac,scatplot)
MODULATIONS = ["QPSK","16QAM","64QAM"];
modulation = ["QPSK","QAM","QAM"];
Baud_rate = num2str(Rs);
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));
TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
TX_SYMB= [repmat(SIG.Xpol.txSymb,10,1), repmat(SIG.Ypol.txSymb,10,1)]; % repeat the bits 10 times to simulate the original transmission in symb for LMS
if r == 1
    M = 4;
    power_norm = 2;
elseif r == 2
    M = 16;
    power_norm = 10;
else
    M = 64;
    power_norm = 42;
end
[X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB, SIG.symbolRate);
[Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB, SIG.symbolRate);


% ----------------Compensation for CD-------------------

[X_CD_rec,Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, SIG.Sps, 2);

[X_CD_rec,Y_CD_rec] = freq_compensation(X_CD_rec, Y_CD_rec, SIG.Sps, SIG.symbolRate);

X_CD_rec = X_CD_rec(65536*8+1:end);
Y_CD_rec = Y_CD_rec(65536*8+1:end);

X_CD_rec = downsample(X_CD_rec, 4);
Y_CD_rec = downsample(Y_CD_rec, 4);

if(rem(length(X_CD_rec),2) ~= 0)

    X_CD_rec = X_CD_rec(2:end);
    Y_CD_rec = Y_CD_rec(2:end);

end

X_Power = mean(abs((X_CD_rec)).^2);
X_CD_rec_norm = X_CD_rec/sqrt(X_Power/power_norm);

Y_Power = mean(abs((Y_CD_rec)).^2);
Y_CD_rec_norm = Y_CD_rec/sqrt(Y_Power/power_norm);

TX_sig = [X_CD_rec_norm, Y_CD_rec_norm];

[X_eq, Y_eq, e_X, e_Y] = LMS_N(TX_sig, modulation(r), EQ_mu, EQ_mu2, EQ_N_tap, TX_SYMB, M, 1*65536);

N=3; % Number of repetitions to drop before calculating the error, to avoid high BER rate because LMS did not converge yet
X_eq = X_eq(N*65536+1:end-100);
Y_eq = Y_eq(N*65536+1:end-100);

X_Power = mean(abs((X_eq)).^2);
X_eq = X_eq/sqrt(X_Power/power_norm);
Y_Power = mean(abs((Y_eq)).^2);
Y_eq = Y_eq/sqrt(Y_Power/power_norm);

[X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = Demapping(X_eq,Y_eq,SIG.Xpol,M);

% Plot for the Distribution of errors
% X_dif = X_demappedBits - TX_BITS_Xpol(1:length(X_demappedBits),:);
% X_dif = sum(abs(X_dif),2);
% figure;plot(cumsum(X_dif));

X_BER = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
Y_BER = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));

Ber_Tot = .5*(min(X_BER)+min(Y_BER));














%--------------------------------------------------------------------------
% Only applies for the "scatterplot" mode of the simulator

if scatplot ~= 0
    num_bins = 50;
    num_points = 1e4;
    jumps = round(.5*length(X_eq)/num_points);
    axlim = 2;
    pointsize = 10;
    figure;
    subplot(2,2,1)
    [density,~,~,binX,binY] = histcounts2(real(X_distorted_AWGN(1:16*jumps:end)), imag(X_distorted_AWGN(1:16*jumps:end)), [num_bins num_bins]);  % 30x30 Gitter für die Dichte
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_distorted_AWGN(1:16*jumps:end)), imag(X_distorted_AWGN(1:16*jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0,0],[-10,10], 'k');
    plot([-10,10],[0,0], 'k');
    xlim([-axlim,axlim]); ylim([-axlim,axlim]);
    colormap(jet);
    title('CD comp. input');
    xlabel('X');
    ylabel('Y');
    grid on;
    axis square;
    hold off;

    subplot(2,2,2)
    [density,~,~,binX,binY] = histcounts2(real(X_CD_rec_norm(1:4*jumps:end)), imag(X_CD_rec_norm(1:4*jumps:end)), [num_bins num_bins]);  % 30x30 Gitter für die Dichte
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_CD_rec_norm(1:4*jumps:end)), imag(X_CD_rec_norm(1:4*jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0,0],[-10,10], 'k');
    plot([-10,10],[0,0], 'k');
    xlim([-axlim,axlim]); ylim([-axlim,axlim]);
    colormap(jet);
    title('CD comp. output');
    xlabel('X');
    ylabel('Y');
    grid on;
    axis square;
    hold off;

    if r==2
        axlim = axlim/1.3;
    end

    subplot(2,2,3)
    [density,~,~,binX,binY] = histcounts2(real(X_out(1:2*jumps:end)), imag(X_out(1:2*jumps:end)), [num_bins num_bins]);  % 30x30 Gitter für die Dichte
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_out(1:2*jumps:end)), imag(X_out(1:2*jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0,0],[-10,10], 'k');
    plot([-10,10],[0,0], 'k');
    xlim([-axlim,axlim]); ylim([-axlim,axlim]);
    colormap(jet);
    title('EQ output');
    xlabel('X');
    ylabel('Y');
    grid on;
    axis square;
    hold off;

    subplot(2,2,4)
    [density,~,~,binX,binY] = histcounts2(real(X_eq(round(end/2):jumps:end)), imag(X_eq(round(end/2):jumps:end)), [num_bins num_bins]);  % 30x30 Gitter für die Dichte
    idx = sub2ind(size(density), binX, binY);
    pointDensity = density(idx);
    scatter(real(X_eq(round(end/2):jumps:end)), imag(X_eq(round(end/2):jumps:end)), pointsize, pointDensity, 'filled');
    hold on;
    plot([0,0],[-10,10], 'k');
    plot([-10,10],[0,0], 'k');
    xlim([-axlim,axlim]); ylim([-axlim,axlim]);
    colormap(jet);
    title('Phase comp. output');
    xlabel('X');
    ylabel('Y');
    grid on;
    axis square;
    hold off;

    sgtitle(sprintf('%s, BER=%0.1e',MODULATIONS(r),Ber_Tot));
end
end

