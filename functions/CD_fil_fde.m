function [data_IQx,data_IQy]=CD_fil_fde(datax, datay, fft_length, dispersion, wavelength, sample_period)
m=size(datax, 1);

if m == 1
datax = datax(:);
datay = datay(:);
end

m = size(datax,1);
if m == 1
datay = datay(:);
end

ndata = size(datax,1);

G = fft_length/4;
L = fft_length/2;

C=3e5 ; 
q=[-fft_length/2: fft_length/2-1]';
omega=2*pi*(q)/(sample_period*fft_length) ;
filt=exp(-j*(dispersion*(wavelength)^2/(4*pi*C))*omega.^2) ;

fil= struct('filter', filt, 'nfft', fft_length, 'window', L, 'ndata', ndata,'overlap', G);
data_IQx = overlap_add(fil,datax) ;
data_IQy = overlap_add(fil,datay) ;


end