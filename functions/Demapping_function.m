function [X_Ber_Tot_CMA, Y_Ber_Tot_CMA] = Demapping_function(X_eq, Y_eq, OSNR_dB, r, index, TX_BITS_Xpol, TX_BITS_Ypol, M, SIG_Xpol_Symb, SIG_Ypol_Symb, MODULATIONS, algorithm)

    X_Power = mean(abs((X_eq)).^2);
    X_eq = X_eq/sqrt(X_Power/10);

    Y_Power = mean(abs((Y_eq)).^2);
    Y_eq = Y_eq/sqrt(Y_Power/10);
      
    X_BER = zeros(1,4);
    Y_BER = zeros(1,4);
    j=1;
    
    for i=0:pi/2:3/2*pi
    
    %fprintf('---------The phase tried is (degrees): %d-----------\n', (mean(i) *180 /pi));
    
    % fprintf('The total phase recovered is (degrees): %d\n', (mean(phEstX+i) *180 /pi));
    
    
    transient_Xpol = abs(finddelay(X_eq(1:65536*2), SIG_Xpol_Symb));
    transient_Ypol = abs(finddelay(Y_eq(1:65536*2), SIG_Ypol_Symb));
    
    if(index==length(OSNR_dB))
        fprintf('Transient Xpol: %d\n', transient_Xpol)
        fprintf('Transient Ypol: %d\n', transient_Ypol)
    end
    
    X_RX = X_eq*exp(1i*i);
    X_RX = X_RX(transient_Xpol+1:end);
    Y_RX = Y_eq*exp(1i*i);
    Y_RX = Y_RX(transient_Ypol+1:end);

    if (index==length(OSNR_dB) && j==1)
        scatterplot(X_RX);
        title(sprintf('%s Xpol after delay recovery %s, OSNR=%d dB',MODULATIONS, algorithm, OSNR_dB(index)));
    end
    
    if r==1
        %fprintf('The tracked moduluation is: QPSK\n');
          [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QPSK_demapping(X_RX, Y_RX);
    %       X_demappedBits = pskdemod(X_RX,M, pi/4*7); %it doesn't demodulate in the same way as our function
    %       X_demappedBits = From_MATLAB_pskdemod(X_demappedBits);
    else
        %fprintf('The tracked moduluation is: 16-QAM\n');
       %     X_demappedBits = qamdemod(X_2Sps(1:2:end),M);
     [X_demappedBits,X_demappedSymb,Y_demappedBits, Y_demappedSymb] = QAM_16_demapping(X_RX,Y_RX);
    end
    
      X_BER(j) = biterr(X_demappedBits, TX_BITS_Xpol(1:length(X_demappedBits),:))/(length(X_demappedBits)*(log2(M)));
      Y_BER(j) = biterr(Y_demappedBits, TX_BITS_Ypol(1:length(Y_demappedBits),:))/(length(Y_demappedBits)*(log2(M)));
    
       j=j+1;
    
    end
    
    X_Ber_Tot_CMA = min(X_BER);
    Y_Ber_Tot_CMA = min(Y_BER);
    fprintf('The BER at %d dB on Xpol is: %.6f\n', OSNR_dB(index), X_Ber_Tot_CMA);
    fprintf('The BER at %d dB on Ypol is: %.6f\n', OSNR_dB(index), Y_Ber_Tot_CMA);
