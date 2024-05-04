function [X_matched,Y_matched] = Matched_filtering(X_pol, Y_pol, Pulse_Shaping)
%----------------Pulse shaping and Downsample the signal-------------------
Pulse_Shaping = downsample(Pulse_Shaping, 4);

X_pol_f = fft(X_pol);
Y_pol_f = fft(Y_pol);
Pulse_f = fft(Pulse_Shaping);

N_Xpol = length(X_pol_f);
N_Ypol = length(Y_pol_f);
N_Pulse = length(Pulse_Shaping);


X_matched_f = X_pol_f(1:floor(N_Xpol/2)+1,:) .* conj(Pulse_f(1:floor(N_Pulse/2)+1));
Y_matched_f = Y_pol_f(1:floor(N_Ypol/2)+1,:) .* conj(Pulse_f(1:floor(N_Pulse/2)+1));

X_matched = ifft(X_matched_f);
Y_matched = ifft(Y_matched_f);

end