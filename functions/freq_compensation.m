function [X_CD_rec,Y_CD_rec] = freq_compensation(TX_Xpol,TX_Ypol, SIG_SpS, SIG_symbolRate)
    
 Bs = SIG_SpS * SIG_symbolRate;
 f = (-Bs/2:Bs/length(TX_Xpol):Bs/2 - Bs/length(TX_Xpol)); %centered in 0
 T = 1/SIG_symbolRate; 
 N = length(TX_Xpol); 
 t = (1:N)* T;

 y = TX_Xpol.^4;
 Y = fftshift(abs(fft(y)));
 f_est = abs((1/4)*f(Y == max(Y)))/SIG_SpS;
 X_CD_rec = TX_Xpol .* exp(-1j*2*pi*f_est*t)';
 
 y = TX_Ypol.^4;
 Y = fftshift(abs(fft(y)));
 f_est = abs((1/4)*f(Y == max(Y)))/SIG_SpS;
 Y_CD_rec = TX_Ypol .* exp(-1j*2*pi*f_est*t)';

end