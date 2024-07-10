function GUI_main_scatterplot(r,Rs, OSNR_dB, delta_nu, rad_sec, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, EQ_N2, CarSync_DampFac)
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "QAM"];
Baud_rate = num2str(Rs);
load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));

if r == 1
    M = 4;
    power_norm = 2;
else
    M = 16;
    power_norm = 10;
end

TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission

% Create delay and phase convolved signals
[X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate);
%add chromatic dispersion
[X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);



%% SIUMULATION

[X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB, SIG.symbolRate);
[Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB, SIG.symbolRate);


% ----------------Compensation for CD-------------------

[X_CD_rec,Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, SIG.Sps, 2);

X_CD_rec = X_CD_rec(65536*8+1:end);
Y_CD_rec = Y_CD_rec(65536*8+1:end);

X_CD_rec = downsample(X_CD_rec, 4);
Y_CD_rec = downsample(Y_CD_rec, 4);

if(rem(length(X_CD_rec),2) ~= 0)

    X_CD_rec = X_CD_rec(2:end);
    Y_CD_rec = Y_CD_rec(2:end);

end

X_Power = mean(abs((X_CD_rec)).^2);
X_CD_rec_norm = X_CD_rec/sqrt(X_Power);

Y_Power = mean(abs((Y_CD_rec)).^2);
Y_CD_rec_norm = Y_CD_rec/sqrt(Y_Power);z

TX_sig = [X_CD_rec_norm, Y_CD_rec_norm];

% %CMA parameters
% if r==1
%     N_tap = 8; %13
%     mu = 8e-3; %7e-3
%     mu2 = 8e-4;
%     N1 = 1e3; %5000
%     N2 = 3e5; %300000
% else
%     N_tap = 11; %13
%     mu = 8e-3; %7e-3
%     mu2 = 8e-4;
%     N1 = 1e3; %5000
%     N2 = 3e5; %300000
% end

[X_out, Y_out, e_X, e_Y] = EQ_func_N(TX_sig, r, EQ_mu, EQ_mu2, EQ_N_tap, EQ_N1, EQ_N2);

cut = floor(EQ_N_tap/2);

X_eq_CMA = X_out;
Y_eq_CMA = Y_out;


%------------------Delay&Phase recovery ---------------------
if r==1
    % CarSync_DampFac = 150;
    carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r),"SamplesPerSymbol", 1, 'DampingFactor', CarSync_DampFac);
    [X_eq, phEstX] = carrSynch(X_eq_CMA);
    [Y_eq, phEstY] = carrSynch(Y_eq_CMA);
else
    % CarSync_DampFac = 31.6;
    carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', CarSync_DampFac);
    %             carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 230, 'NormalizedLoopBandwidth',1e-3);
    %             [X_test, phEstX_test] = carrSynch(X_eq_CMA);
    [X_eq, phEstX] = carrSynch(X_eq_CMA);

    carrSynch2 = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', CarSync_DampFac);
    %             [Y_test, phEstY_test] = carrSynch(Y_eq_CMA);
    [Y_eq, phEstY] = carrSynch2(Y_eq_CMA);

    %             [X_eq,~] = BPS_N(X_eq_CMA, 50, M, power_norm);
    %             [Y_eq,~] = BPS_N(Y_eq_CMA, 50, M, power_norm);
end

X_BER = zeros(1,4);
Y_BER = zeros(1,4);
j=1;

for i=0:pi/2:3/2*pi
    X_RX = X_eq*exp(1i*i);
    [~,transient_Xpol] = max(abs(xcorr(X_RX(1:65536*2), SIG.Xpol.txSymb)));
    X_RX = X_RX(transient_Xpol+1:end);
    Y_RX = Y_eq*exp(1i*i);
    [~,transient_Ypol] = max(abs(xcorr(Y_RX(1:65536*2), SIG.Ypol.txSymb)));
    Y_RX = Y_RX(transient_Ypol+1:end);

    if r==2
        X_Power = mean(abs((X_RX)).^2);
        X_RX = X_RX/sqrt(X_Power/power_norm);

        Y_Power = mean(abs((Y_RX)).^2);
        Y_RX = Y_RX/sqrt(Y_Power/power_norm);
    end

    if r==1
        [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);

    else

        [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_RX,Y_RX);
    end

    X_BER(j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
    Y_BER(j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));

    j=j+1;

end

Ber_Tot = .5*(min(X_BER)+min(Y_BER));


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



