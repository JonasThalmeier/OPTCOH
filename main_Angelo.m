clear;
close all;

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
OSNR_dB = 15;
seq_lenght = length(SIG.Xpol.txSymb);

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
%%

for index = 1:length(OSNR_dB)
    % Adding Noise
    [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
    [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);

    % Recovering Chromatic Dispersion
    [X_CD_rec,Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, SIG.Sps, 2);

    % Downsampling and Matched Flitering
    X_CD_rec = downsample(X_CD_rec, 4);
    Y_CD_rec = downsample(Y_CD_rec, 4);
    
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
        XY_vit = vit_n_vit(XY_eq, Delta_nu, SIG.symbolRate, OSNR_dB(index), Es, Npol, M, windowlen);
        X_eq = XY_vit(:,1);
        Y_eq = XY_vit(:,2);
    else
        mu = 1e-3;
        NTaps = 9;
        N1 = 1e3;
        N2 = 1e4;
        XY_eq = EQ_func(X_CD_rec,Y_CD_rec,mu,NTaps,"RDE",N1,N2);
        X_eq=XY_eq(:,1);
        Y_eq=XY_eq(:,2);
    end
    
    scatterplot(X_eq(40000:end));
    X_eq=X_eq(40000:end);
    
    %considering pi/2 rotation
    p=pi/2;
    B=100; %B must be even
    b=-B/2:1:B/2-1;
    Theta_b=b/B*p;
    distances=zeros(length(Theta_b),1);
    costFunction = @(x, y) abs(x - y).^2; % Cost function (Mean Squared Error)
    
    X_Power = mean(abs((X_eq)).^2);
    X_eq = X_eq/sqrt(X_Power/10);
    Y_Power = mean(abs((Y_eq)).^2);
    Y_eq = Y_eq/sqrt(Y_Power/10);
    transient_Xpol = abs(finddelay(X_eq(1:2*65536), SIG.Xpol.txSymb));
    transient_Ypol = abs(finddelay(Y_eq(1:2*65536), SIG.Ypol.txSymb));
    X_RX = X_eq(transient_Xpol+1:end);
    signal_compensated=zeros(length(X_RX),1);
 
    for k=1:length(X_RX)
        
      Rotations=zeros(length(Theta_b),1);
      
      for i=1:length(Theta_b) 
        Rotations(i,:) = X_RX(k,:)*exp(1i*Theta_b(i));
      end
        
      if r==1
          fprintf('The tracked moduluation is: QPSK\n');
          [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);
      else
          [X_demappedBits,X_demappedSymb,demappedConfig_Xpol] = QAM_16_demapping(Rotations);
      end
      
      distances=costFunction(Rotations, demappedConfig_Xpol);
      
%       zeros_vector = zeros(length(distances) / 2, 1);
%       
%       distances_1=[distances; zeros_vector];
%       %
%       N_window=30;
%       
%       for w=1:length(distances)
%         distances(w,:)=sum(distances_1(w:w+2*N_window,:));
%       end
      
      %distances=distances(1:length(distances)-2*N_window,:);
      
      [M,I]=min(distances);
      signal_compensated(k,:)=Rotations(I,:);
      
    end
    
    scatterplot(signal_compensated)
    
   % Finding the constellations right orientation
    X_BER = zeros(1,4);
    Y_BER = zeros(1,4);
    j=1;
    M=16
    
    for i=0:pi/2:3/2*pi
        
        signal_compensated_1=signal_compensated*exp(1i*i);
        
        if r==1
            fprintf('The tracked moduluation is: QPSK\n');
            [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);
        else
            [X_demappedBits,X_demappedSymb,demappedConfig_Xpol] = QAM_16_demapping(signal_compensated_1);
        end
        
        X_BER(:,j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
        %Y_BER(k,j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));
        j=j+1;
    end
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