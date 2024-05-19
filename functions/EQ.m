function [XPol_out] = EQ(Xpol_in,mu)
num_taps = 9;  % Number of filter taps
fir = zeros(1, num_taps);  % Initialize the FIR filter coefficients
N = length(Xpol_in);  % Length of the input signal
XPol_out = zeros(N - num_taps + 1, 1);  % Output signal
fir(ceil(num_taps / 2)) = 1;  % Initialize the center tap to 1
XPol_out(1) = fir*Xpol_in(1:num_taps);
for idx=1:N-num_taps
    % Calculate the error term for the CMA update
    er = conj(XPol_out(idx))*(abs(XPol_out(idx))^2-1);
    % Update the FIR filter coefficients using CMA
    fir = fir - mu*er.*Xpol_in(idx:idx+num_taps-1)';
    XPol_out(idx+1) = fir*Xpol_in(idx+1:idx+num_taps);
end