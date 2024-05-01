function [X_CD,Y_CD] = Chromatic_Dispersion(SIG_Xpol, SIG_Ypol, SpS, flag)

f0 = 193.4e12; % Considering 1550 nm as the center frequency for c-band
c = 299792458; % Speed of light in vacuum (m/s)
lambda = c / f0; % Wavelength of light in meters
D = 17; % Dispersion parameter in ps/(nm*km)
L = 50; % Length of the fiber in km

beta_2 = -(D * lambda^2) / (2 * pi * c); 

Bs = SpS * 64e9; 
sample_period=1/Bs;

f = (-Bs/2:Bs/length(SIG_Xpol):Bs/2 - Bs/length(SIG_Xpol)).';

if (flag==1)
    H_cd = exp(-1i * 2 * pi^2 * beta_2 * (f - f0).^2 * L);
    
    SIG_Xpol_txSig_disp_f = fft(SIG_Xpol).* H_cd;
    SIG_Ypol_txSig_disp_f = fft(SIG_Ypol).* H_cd;
    
    X_CD = ifft(SIG_Xpol_txSig_disp_f);
    Y_CD = ifft(SIG_Ypol_txSig_disp_f);
    
elseif (flag==2)
    fft_length=2^15; %try different parameters for it
    [X_CD,Y_CD]=CD_fil_fde(SIG_Xpol, SIG_Ypol, fft_length, D, lambda, sample_period);
    
end