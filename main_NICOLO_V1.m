clear;
close all;
clc;

% Load the .mat file
MODULATIONS = ["QPSK","16QAM"];
modulation = ["QPSK" "QAM"];
% r = randi([1, 2], 1); % Get a 1 or 2 randomly.
r = 2;
fprintf('The transmitted moduluation is: %s\n', modulation(r));
load(strcat('C:\Users\utente\Documents\GitHub\OPTCOH\TXsequences\TXsequence_', MODULATIONS(r) , '_64GBaud.mat'));
if r == 1
    M = 4;
    power_norm = 2;
else
    M = 16;
    power_norm = 10;
end

TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission

rad_pol = [1e1 1e2 1e3 1e4 1e5 1e6];
Delta_Pol = zeros(1,length(rad_pol));

for index_rad_pol = 1:length(rad_pol)

    %% IMPAIRMENTS PART
    % Create delay and phase convolved signals
    [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, 5e2, rad_pol(index_rad_pol));
    
    %----------------Jones Matrix------------------
    
    %dynamic
    
    % amplitude = 1;
    % period = 8.5; %with 4.5 is a rotation from 0 to 90°(almost), and it is not recovered
    % t = linspace(0, 1, length(X_distorted));
    % phi = pi/3;
    % 
    % pol_xx = amplitude * cos(2 * pi * t / period);
    % pol_xy = exp(-1i*phi)*amplitude * sin(2 * pi * t / period);
    % pol_yx = -exp(1i*phi)*amplitude * sin(2 * pi * t / period);
    % pol_yy = amplitude * cos(2 * pi * t / period);
    % 
    % 
    % X_distorted_Jones = pol_xx'.*X_distorted + pol_xy'.*Y_distorted;
    % Y_distorted_Jones = pol_yx'.*X_distorted + pol_yy'.*Y_distorted;
    % 
    % 
    % plot(t, pol_xx, 'b-', 'LineWidth', 2);
    % xlabel('samples');
    % ylabel('cos(n)');
    % grid on;
    
    %     theta = linspace(pi/4,pi/5,length(X_distorted));
    %     phi = linspace(pi/3,pi/4,length(X_distorted));
    % 
    %     Sig_distorted_Jones = zeros(2,length(X_distorted));
    %     XY_Dist = [X_distorted'; Y_distorted'];
    %     R = zeros(2,2,length(theta));
    % 
    %     for l=1:length(theta)
    %         R(:,:,l) = [cos(theta(l)) exp(-1i*phi(l))*sin(theta(l));-exp(1i*phi(l))*sin(theta(l)) cos(theta(l))];  
    %         Sig_distorted_Jones(:,l) = R(:,:,l)*XY_Dist(:,l);
    %     end
    
    %static
    % theta = pi/4;
    % phi = pi/3;
    % 
    % XY_Dist = [X_distorted'; Y_distorted'];
    % 
    % R = [cos(theta) exp(-1i*phi)*sin(theta);-exp(1i*phi)*sin(theta) cos(theta)];  
    % Sig_distorted_Jones = R*XY_Dist;
    % 
    % 
    % X_distorted_Jones = Sig_distorted_Jones(1,:)';
    % Y_distorted_Jones = Sig_distorted_Jones(2,:)';
    
    %     X_distorted_AWGN_Jones = X_distorted_AWGN * cos(theta) + Y_distorted_AWGN * exp(-1i*phi)*sin(theta);
    %     Y_distorted_AWGN_Jones = Y_distorted_AWGN * cos(theta) - X_distorted_AWGN * exp(1i*phi)*sin(theta);
    
    %add chromatic dispersion
    [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1);
    
    %%
    % Adding the noise
    if r==1
        OSNR_dB = 8:11;
    else
        OSNR_dB = 16;
    end
    
    X_Ber_Tot = zeros(1,length(OSNR_dB));
    Y_Ber_Tot = zeros(1,length(OSNR_dB));
    
    %% CMA mine
    
    for index = 1:length(OSNR_dB)
    
        [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
        [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
        
    
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
        Y_CD_rec_norm = Y_CD_rec/sqrt(Y_Power);
    
        TX_sig = [X_CD_rec_norm'; Y_CD_rec_norm'];
        
        %CMA parameters
        if r==1
            N_tap = 15;
            mu = 1e-4;
            mu2 = 8e-5;
            N1 = 5000;
            N2 = 300000;
        else
            N_tap = 11; %13
            mu = 8e-3; %7e-3
            mu2 = 8e-4;
            N1 = 5000;
            N2 = 300000;
        end
    
        [X_out, Y_out, e_X, e_Y] = EQ_func_N(TX_sig,r,mu,mu2,N_tap,N1,N2);
        
        X_out = X_out';
        Y_out = Y_out';
        cut = 15;
    
    
        figure(), plot(e_X), title(sprintf('CMA error Xpol %d dB', OSNR_dB(index))), hold on, xline(cut, 'LineStyle','--', 'Color','r'), xlabel('N_samples');
        figure(), plot(e_Y), title(sprintf('CMA error Ypol %d dB', OSNR_dB(index))), hold on, xline(cut, 'LineStyle','--', 'Color','r'), xlabel('N_samples');
    
        
        
        X_eq_CMA = X_out(1:end-cut,:);
        Y_eq_CMA = Y_out(1:end-cut,:);
    
    
     %------------------Delay&Phase recovery ---------------------
        
        if r==1
            carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r),"SamplesPerSymbol", 1, 'DampingFactor', 150);
            [X_eq, phEstX] = carrSynch(X_eq_CMA);
            [Y_eq, phEstY] = carrSynch(Y_eq_CMA);
        else
    %         carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 450, 'NormalizedLoopBandwidth',.0005);%, 'ModulationPhaseOffset','Custom', 'CustomPhaseOffset', -pi/7);
    %         %[X_test, phEstX_test] = carrSynch(X_eq_CMA);
    %         [X_eq, phEstX] = carrSynch(X_eq_CMA);
    
    %         carrSynch2 = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 240, 'NormalizedLoopBandwidth',.0002);%, 'ModulationPhaseOffset','Custom', 'CustomPhaseOffset', -pi/7);
    %         [Y_test, phEstY_test] = carrSynch(Y_eq_CMA);
    %         [Y_eq, phEstY] = carrSynch2(Y_eq_CMA);
    
            [X_eq,~] = BPS_N(X_eq_CMA, 50, M, power_norm); 
            [Y_eq,~] = BPS_N(Y_eq_CMA, 50, M, power_norm);
        end
        
        
        if (index==length(OSNR_dB))
            scatterplot(X_eq);
            title(sprintf('%s Xpol after MF and phase recovery, OSNR=%d dB',MODULATIONS(r), OSNR_dB(index)));
        end
    
        X_BER = zeros(1,4);
        Y_BER = zeros(1,4);
        j=1;
        
        for i=0:pi/2:3/2*pi
        
    %     fprintf('---------The phase tried is (degrees): %d-----------\n', (mean(i) *180 /pi));
        
        % fprintf('The total phase recovered is (degrees): %d\n', (mean(phEstX+i) *180 /pi));
           
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
    
         if(index==length(OSNR_dB))
            fprintf('Transient Xpol: %d\n', transient_Xpol)
            fprintf('Transient Ypol: %d\n', transient_Ypol)
        end
    
        if (index==length(OSNR_dB) && j==1)
            scatterplot(X_RX);
            title(sprintf('%s Xpol after delay recovery, OSNR=%d dB',MODULATIONS(r), OSNR_dB(index)));
        end
        
        if r==1
            %fprintf('The tracked moduluation is: QPSK\n');
              [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);
      
        else
            %fprintf('The tracked moduluation is: 16-QAM\n');
    %         MyConst = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];
    %         X_demappedBits = qamdemod(X_RX, M, MyConst, OutputType='bit', PlotConstellation=true);
    %         N = length(X_demappedBits)/4;
    %         X_demappedBits = reshape(X_demappedBits, 4, N).';
    %        
            [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_RX,Y_RX);
        end
        
          X_BER(j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
          Y_BER(j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));
        
           j=j+1;
        
        end
        
        X_Ber_Tot(index) = min(X_BER);
        Y_Ber_Tot(index) = min(Y_BER);
        fprintf('The BER on Xpol is: %.6f\n', X_Ber_Tot(index));
        fprintf('The BER on Ypol is: %.6f\n', Y_Ber_Tot(index));
    end
    BER_Tot = (X_Ber_Tot+Y_Ber_Tot)/2;
    OSNR_th = 10*log10(10*erfinv(-8/3*BER_Tot+1)^2);
    Delta_Pol(index_rad_pol) = OSNR_dB - OSNR_th;
end

figure(), semilogx(Delta_Pol);

fprintf('\n');
% %%
% %--------------------MATCHED FILTER ----------------------------------------
% 
% for index = 1:length(OSNR_dB)
% 
%     [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
%     [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
%     
%     if (index==length(OSNR_dB))
%         scatterplot(X_distorted_AWGN(5:8:end));
%         title(sprintf('%s constellation Xpol after CD and WGN, OSNR=%d dB',MODULATIONS(r), OSNR_dB(index)));
%     end
% 
%     % ----------------Compensation for CD-------------------
%     
%     [X_CD_rec,Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, SIG.Sps, 2);
% 
%     if (index==length(OSNR_dB))
%         scatterplot(X_CD_rec(5:8:end));
%         title(sprintf('%s constellation Xpol with WGN, OSNR=%d dB', MODULATIONS(r), OSNR_dB(index)));
%     end
% 
%     X_CDrec_powers_per_sample = zeros(1, SIG.Sps);
%     Y_CDrec_powers_per_sample = zeros(1, SIG.Sps);
% 
%     for rep = (1:SIG.Sps)
%         X_CDrec_powers_per_sample(rep) = mean(abs((X_CD_rec(rep:8:end))).^2);
%         Y_CDrec_powers_per_sample(rep) = mean(abs((Y_CD_rec(rep:8:end))).^2);
%     end
% 
%     [X_Max_Power, X_index_Max_Power] = min(X_CDrec_powers_per_sample);
%     [Y_Max_Power, Y_index_Max_Power] = min(Y_CDrec_powers_per_sample);
% 
%     X_CD_rec = X_CD_rec(X_index_Max_Power:end);
%     Y_CD_rec = Y_CD_rec(Y_index_Max_Power:end);
%  
%     
%     %------------------ Matched Flitering ---------------------
%     X_CD_rec = downsample(X_CD_rec, 4);
%     Y_CD_rec = downsample(Y_CD_rec, 4);
%    
%     [X_matched,Y_matched] = Matched_filtering(X_CD_rec, Y_CD_rec, PulseShaping.b_coeff);
% 
%       
%     if (index==length(OSNR_dB))
%         scatterplot(X_matched(1:2:end));
%         title(sprintf('%s Xpol after matched filter, OSNR=%d dB',MODULATIONS(r), OSNR_dB(index)));
%     end
% 
%     %------------------Delay&Phase recovery ---------------------
%     
%     if r==1
%         carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r),"SamplesPerSymbol", 1, 'DampingFactor', 150);
%         [X_eq, phEstX] = carrSynch(X_matched(1:2:end));
%         [Y_eq, phEstY] = carrSynch(Y_matched(1:2:end));
%     else
%         carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r), "SamplesPerSymbol", 1,'DampingFactor', 150, 'NormalizedLoopBandwidth',.0005);%, 'ModulationPhaseOffset','Custom', 'CustomPhaseOffset', -pi/7);
%         [X_test, phEstX_test] = carrSynch(X_matched(1:2:end));
%         [X_eq, phEstX] = carrSynch(X_matched(1:2:end));
%         [Y_eq, phEstY] = carrSynch(Y_matched(1:2:end));
%     end
%     
%     X_Power = mean(abs((X_eq)).^2);
%     X_eq = X_eq/sqrt(X_Power/10);
% 
%     Y_Power = mean(abs((Y_eq)).^2);
%     Y_eq = Y_eq/sqrt(Y_Power/10);
% 
%     if (index==length(OSNR_dB))
%         scatterplot(X_eq);
%         title(sprintf('%s Xpol after MF and phase recovery, OSNR=%d dB',MODULATIONS(r), OSNR_dB(index)));
%     end
% 
% %     X_eq = X_matched(1:2:end);
%     %
%     X_BER = zeros(1,4);
%     Y_BER = zeros(1,4);
%     j=1;
%     
%     for i=0:pi/2:3/2*pi
%     
% %     fprintf('---------The phase tried is (degrees): %d-----------\n', (mean(i) *180 /pi));
%     
%     % fprintf('The total phase recovered is (degrees): %d\n', (mean(phEstX+i) *180 /pi));
%     
%     transient_Xpol = abs(finddelay(X_eq(1:65536), SIG.Xpol.txSymb));
%     transient_Ypol = abs(finddelay(Y_eq(1:65536), SIG.Ypol.txSymb));
%     
%     if(index==length(OSNR_dB))
%         fprintf('Transient Xpol: %d\n', transient_Xpol)
%         fprintf('Transient Ypol: %d\n', transient_Ypol)
%     end
%     
%     X_RX = X_eq*exp(1i*i);
%     X_RX = X_RX(transient_Xpol+1:end);
%     Y_RX = Y_eq*exp(1i*i);
%     Y_RX = Y_RX(transient_Ypol+1:end);
% 
%     if (index==length(OSNR_dB) && j==1)
%         scatterplot(X_RX);
%         title(sprintf('%s Xpol after delay recovery, OSNR=%d dB',MODULATIONS(r), OSNR_dB(index)));
%     end
%     
%     if r==1
%         %fprintf('The tracked moduluation is: QPSK\n');
%           [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);
%   
%     else
%         %fprintf('The tracked moduluation is: 16-QAM\n');
% %         MyConst = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];
% %         X_demappedBits = qamdemod(X_RX, M, MyConst, OutputType='bit', PlotConstellation=true);
% %         N = length(X_demappedBits)/4;
% %         X_demappedBits = reshape(X_demappedBits, 4, N).';
% %        
%         [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_RX,Y_RX);
%     end
%     
%       X_BER(j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
%       Y_BER(j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));
%     
%        j=j+1;
%     
%     end
%     
%     X_Ber_Tot(index) = min(X_BER);
%     Y_Ber_Tot(index) = min(Y_BER);
%     fprintf('The BER on Xpol is: %.6f\n', X_Ber_Tot(index));
%     fprintf('The BER on Ypol is: %.6f\n', Y_Ber_Tot(index));
% end
% 
% 
% fprintf('\n');
% 
% %% 
% %------------------------- CMA/LMS ----------------------------
% 
% fprintf('-----------------CMA/LMS------------------\n');
% 
% X_Ber_Tot_CMA = zeros(1,length(OSNR_dB));
% Y_Ber_Tot_CMA = zeros(1,length(OSNR_dB));
% X_Ber_Tot_LMS = zeros(1,length(OSNR_dB));
% Y_Ber_Tot_LMS = zeros(1,length(OSNR_dB));
% 
% for index = 1:length(OSNR_dB)
% 
%     [X_distorted_AWGN, NoiseX] = WGN_Noise_Generation(X_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
%     [Y_distorted_AWGN, NoiseY] = WGN_Noise_Generation(Y_CD, SIG.Sps, M, OSNR_dB(index), SIG.symbolRate);
%     
%     % ----------------Compensation for CD-------------------
%     
%     [X_CD_rec,Y_CD_rec] = Chromatic_Dispersion(X_distorted_AWGN, Y_distorted_AWGN, SIG.Sps, 2);
%     
%     %---------------------------EQ---------------------------------------------
%     %From 15 numtaps are the best for qpsk, also for qam, but 60 and 120 seem
%     %better. The reference tap for qpsk are 1,3,65,121.
% 
%     X_CD_rec = X_CD_rec(65536*8+1:end);
%     Y_CD_rec = Y_CD_rec(65536*8+1:end);
%    
%     X_CD_rec = downsample(X_CD_rec, 4);
%     Y_CD_rec = downsample(Y_CD_rec, 4);
% 
%     if(rem(length(X_CD_rec),2) ~= 0)
% 
%         X_CD_rec = X_CD_rec(2:end);
%         Y_CD_rec = Y_CD_rec(2:end);
% 
%     end
% 
%     X_Power = mean(abs((X_CD_rec)).^2);
%     X_CD_rec_norm = X_CD_rec/sqrt(X_Power/power_norm);
% 
%     Y_Power = mean(abs((Y_CD_rec)).^2);
%     Y_CD_rec_norm = Y_CD_rec/sqrt(Y_Power/power_norm);
% 
%     if(rem(length(X_CD_rec),2) ~= 0)
% 
%         X_CD_rec_norm = X_CD_rec_norm(2:end);
%         Y_CD_rec_norm = Y_CD_rec_norm(2:end);
% 
%     end
%   
%     if r == 1
%   
%         %-------------CMA------------------
%         constellation = pskmod(0:3, 4, pi/4);
%         stepsize = 6e-4; %best result
%         numtaps = 9;
%         referencetap = (numtaps-1)/2;
%         algorithm = 'CMA';
% 
%         EQ = comm.LinearEqualizer('Algorithm', algorithm, 'StepSize', stepsize,'NumTaps', numtaps, 'InputSamplesPerSymbol', 2, 'Constellation', constellation, 'ReferenceTap', referencetap);
%         [X_matched_CMA,errX_CMA] = EQ(X_CD_rec(1:end-mod(length(X_CD_rec_norm),2))); 
%         [Y_matched_CMA,errY_CMA] = EQ(Y_CD_rec(1:end-mod(length(Y_CD_rec_norm),2)));  
% 
%         if (index==length(OSNR_dB))
%             scatterplot(X_matched_CMA);
%             title(sprintf('%s constellation of Xpol after %s, OSNR= %d dB', MODULATIONS(r), algorithm, OSNR_dB(index)));
%         end       
% 
%         %------------------Delay&Phase recovery ---------------------
% 
%         carrSynch = comm.CarrierSynchronizer("Modulation", modulation(r),"SamplesPerSymbol", 1, 'DampingFactor', 150);
%         [X_eq, phEstX] = carrSynch(X_matched_CMA);
%         [Y_eq, phEstY] = carrSynch(Y_matched_CMA);
%         
%         if (index==length(OSNR_dB))
%             scatterplot(X_eq);
%             title(sprintf('%s constellation of Xpol after phase recovery %s, OSNR= %d dB', MODULATIONS(r), algorithm, OSNR_dB(index)));
%         end
% 
%         [X_Ber_Tot_CMA(index), Y_Ber_Tot_CMA(index)] = Demapping_function(X_eq, Y_eq, OSNR_dB, r, index, TX_BITS_Xpol, TX_BITS_Ypol, M, SIG.Xpol.txSymb, SIG.Ypol.txSymb, MODULATIONS(r), algorithm);
% 
% 
%         %---------LMS------------------
%         constellation = pskmod(0:3, 4, pi/4);
%         stepsize = 2e-4; %best result
%         numtaps = 8;
%         referencetap = floor((numtaps-1)/2);
%         algorithm = 'LMS';
% 
%         EQ = comm.LinearEqualizer('Algorithm', algorithm, 'StepSize', stepsize,'NumTaps', numtaps, 'InputSamplesPerSymbol', 2, 'ReferenceTap', referencetap, 'Constellation', constellation); % 'Constellation', constellation,
%         [X_test_LMS,errX_test_LMS] = EQ(X_CD_rec_norm, SIG.Xpol.txSymb);
%         [X_matched_LMS,errX_LMS] = EQ(X_CD_rec_norm, SIG.Xpol.txSymb);
%         [Y_test_LMS,errY_test_LMS] = EQ(Y_CD_rec_norm, SIG.Ypol.txSymb);
%         [Y_matched_LMS,errY_LMS] = EQ(Y_CD_rec_norm, SIG.Ypol.txSymb);
% 
%         X_matched_LMS=X_matched_LMS(80000:end);
%         Y_matched_LMS=Y_matched_LMS(80000:end);
% 
%         if (index==length(OSNR_dB))
%             scatterplot(X_matched_LMS);
%             title(sprintf('%s constellation of Xpol after %s, OSNR= %d dB',MODULATIONS(r), algorithm, OSNR_dB(index)));
%         end 
%         
% 
%         X_Power = mean(abs((X_matched_LMS)).^2);
%         X_matched_LMS_norm = X_matched_LMS/sqrt(X_Power/2);
% 
%         Y_Power = mean(abs((Y_matched_LMS)).^2);
%         Y_matched_LMS_norm = Y_matched_LMS/sqrt(Y_Power/2);
% 
%         [X_Ber_Tot_LMS(index), Y_Ber_Tot_LMS(index)] = Demapping_function(X_matched_LMS_norm, Y_matched_LMS_norm, OSNR_dB, r, index, TX_BITS_Xpol, TX_BITS_Ypol, M, SIG.Xpol.txSymb, SIG.Ypol.txSymb, MODULATIONS(r), algorithm);
% 
% 
%     else
% 
%         MyConst = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];
%         M = 16; % Size of the QAM constellation
%         constellation = qammod(0:15, M, MyConst);
%         stepsize = 2.15e-5;%7.5025e-4;
%         numtaps = 11; 
%         referencetap = floor((numtaps-1)/2);
%         algorithm = 'LMS';
%         
%         EQ = comm.LinearEqualizer('Algorithm', algorithm, 'StepSize', stepsize,'NumTaps', numtaps, 'InputSamplesPerSymbol', 2, 'ReferenceTap', referencetap, 'Constellation', constellation); % 'Constellation', constellation,
%         [X_matched,errX] = EQ(X_CD_rec_norm, SIG.Xpol.txSymb);
%         [Y_matched,errY] = EQ(Y_CD_rec_norm, SIG.Ypol.txSymb);
% 
%         X_matched=X_matched(0.85e5:end);
%         Y_matched=Y_matched(1e5:end);
% 
%         if (index==length(OSNR_dB))
%             scatterplot(X_matched);
%             title(sprintf('%s constellation of Xpol after %s %d dB of WGN',MODULATIONS(r), algorithm, OSNR_dB(index)));
%         end
% 
%         
%         [X_Ber_Tot_LMS(index), Y_Ber_Tot_LMS(index)] = Demapping_function(X_matched, Y_matched, OSNR_dB, r, index, TX_BITS_Xpol, TX_BITS_Ypol, M, SIG.Xpol.txSymb, SIG.Ypol.txSymb, MODULATIONS(r), algorithm);
% 
%     end
% end    
% 

%%

%------------------FIGURES-------------

% if r==1
%     BER_TH = 0.5 * erfc(sqrt(10.^(OSNR_dB/10)/2));
% 
%     BER_MED_MF = 0.5 * (X_Ber_Tot + Y_Ber_Tot);
% %     BER_MED_CMA = 0.5 * (X_Ber_Tot_CMA + Y_Ber_Tot_CMA);
% %     BER_MED_LMS = 0.5 * (X_Ber_Tot_LMS + Y_Ber_Tot_LMS);
% 
%     figure();
%     semilogy(OSNR_dB, BER_TH, 'r', 'LineWidth', 1);
%     xlim([min(OSNR_dB),max(OSNR_dB)]);
%     grid on;
%     hold on;
%     semilogy(OSNR_dB, X_Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineWidth', 1, 'LineStyle','-.');
% %     semilogy(OSNR_dB, X_Ber_Tot_CMA, 'Marker','o', 'Color', 'b', 'LineStyle','-.');
% %     semilogy(OSNR_dB, X_Ber_Tot_LMS, 'Marker','o', 'Color', 'm', 'LineStyle','-.');    
%     title(sprintf('%s BER curve of Xpol',MODULATIONS(r)));
%     legend('Theoretical BER', 'Simulated BER - Matched filter','Simulated BER - CMA', 'Simulated BER - LMS', 'Interpreter', 'latex');
%     xlabel('OSNR [dB]', 'Interpreter','latex');
%     hold off;
%     
%     figure();
%     semilogy(OSNR_dB,BER_TH, 'r','LineWidth', 1); 
%     xlim([min(OSNR_dB),max(OSNR_dB)]);
%     grid on;
%     hold on;
%     semilogy(OSNR_dB,Y_Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineWidth', 1, 'LineStyle','-.');
%     semilogy(OSNR_dB, Y_Ber_Tot_CMA, 'Marker','o', 'Color', 'b', 'LineStyle','-.');
%     semilogy(OSNR_dB, Y_Ber_Tot_LMS, 'Marker','o', 'Color', 'm', 'LineStyle','-.');
%     title(sprintf('%s BER curve of Ypol',MODULATIONS(r)));
%     legend('Theoretical BER', 'Simulated BER - Matched filter', 'Simulated BER - CMA', 'Simulated BER - LMS', 'Interpreter', 'latex');
%     xlabel('OSNR [dB]', 'Interpreter','latex');
%     hold off;
% 
%     figure();
%     semilogy(OSNR_dB, BER_TH, 'r', 'LineWidth', 1);
%     xlim([min(OSNR_dB),max(OSNR_dB)]);
%     grid on;
%     hold on;
%     semilogy(OSNR_dB, BER_MED_MF, 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
%     semilogy(OSNR_dB, BER_MED_CMA, 'Marker','o', 'Color', 'b', 'LineStyle','-.');
%     semilogy(OSNR_dB, BER_MED_LMS, 'Marker','o', 'Color', 'm', 'LineStyle','-.');
%     title(sprintf('%s BER curve',MODULATIONS(r)));
%     legend('Theoretical BER', 'Simulated BER - Matched filter','Simulated BER - CMA', 'Simulated BER - LMS', 'Interpreter', 'latex');
%     xlabel('OSNR [dB]', 'Interpreter','latex');
%     hold off;

% else
% 
%     BER_TH = 3/8 * erfc(sqrt(10.^(OSNR_dB/10)/10));
% 
%     BER_MED_MF = 0.5 * (X_Ber_Tot + Y_Ber_Tot);
% %     BER_MED_LMS = 0.5 * (X_Ber_Tot_LMS + Y_Ber_Tot_LMS);
% 
%     figure();
%     semilogy(OSNR_dB, BER_TH, 'r');
%     xlim([min(OSNR_dB),max(OSNR_dB)]);
%     grid on;
%     hold on;
%     semilogy(OSNR_dB, X_Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
% %     semilogy(OSNR_dB, X_Ber_Tot_LMS, 'Marker','o', 'Color', 'b', 'LineStyle','-.');
%     title(sprintf('%s BER curve of Xpol',MODULATIONS(r)));
% %     legend('Theoretical BER', 'Simulated BER - Matched filter',sprintf('Simulated BER - %s', algorithm), 'Interpreter', 'latex');
%     xlabel('OSNR [dB]', 'Interpreter','latex');
%     hold off;
%     
%     figure();
%     semilogy(OSNR_dB,BER_TH, 'r');
%     xlim([min(OSNR_dB),max(OSNR_dB)]);
%     grid on;
%     hold on;
%     semilogy(OSNR_dB,Y_Ber_Tot, 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
%     semilogy(OSNR_dB, Y_Ber_Tot_LMS, 'Marker','o', 'Color', 'b', 'LineStyle','-.');
%     semilogy(OSNR_dB,BER_TH, 'r');
%     title(sprintf('%s BER curve of Ypol',MODULATIONS(r)));
%     legend('Theoretical BER', 'Simulated BER - Matched filter', sprintf('Simulated BER - %s', algorithm), 'Interpreter', 'latex');
%     xlabel('OSNR [dB]', 'Interpreter','latex');
%     hold off;   

%     figure();
%     semilogy(OSNR_dB, BER_TH, 'r', 'LineWidth', 1);
%     xlim([min(OSNR_dB),max(OSNR_dB)]);
%     grid on;
%     hold on;
%     semilogy(OSNR_dB, BER_MED_MF, 'Marker','o', 'Color', "#77AC30", 'LineStyle','-.', 'LineWidth', 1);
%     semilogy(OSNR_dB, BER_MED_LMS, 'Marker','o', 'Color', 'b', 'LineStyle','-.');
%     title(sprintf('%s BER curve',MODULATIONS(r)));
%     legend('Theoretical BER', 'Simulated BER - Matched filter',sprintf('Simulated BER - %s', algorithm), 'Interpreter', 'latex');
%     xlabel('OSNR [dB]', 'Interpreter','latex');
%     hold off;
%     
% end

