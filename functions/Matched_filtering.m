function [X_matched,Y_matched] = Matched_filtering(X_pol, Y_pol, Pulse_Shaping)
%----------------Pulse shaping and Downsample the signal-------------------
Pulse_Shaping = downsample(Pulse_Shaping, 4);

% X_pol_f = ifft(fft(X_pol));
% Y_pol_f = ifft(fft(Y_pol));
% Pulse_f = ifft(fft(Pulse_Shaping));
% 
% N_Xpol = length(X_pol_f);
% N_Ypol = length(Y_pol_f);
% N_Pulse = length(Pulse_Shaping);

% 
% X_matched_f = X_pol_f .* conj(Pulse_f);
% Y_matched_f = Y_pol_f .* conj(Pulse_f);
 
% X_matched = ifft(X_matched_f);
% Y_matched = ifft(Y_matched_f);

X_matched = fftfilt(conj(fliplr(Pulse_Shaping)), X_pol);
Y_matched = fftfilt(conj(fliplr(Pulse_Shaping)), Y_pol);


end