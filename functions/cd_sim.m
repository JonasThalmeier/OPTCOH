function [Xpol_dc, Ypol_dc] = cd_sim(Xpol, Ypol, Sps, Rb, D, L, lambda0)
    % applyChromaticDispersionFD Applies chromatic dispersion to a signal
    % using the frequency domain approach, considering samples per symbol and
    % symbol rate.
    %
    % Inputs:
    %   Xpol    - Input time-domain signal
    %   Ypol    - Input time-domain signal
    %   Sps     - Samples per symbol
    %   Rb      - Symbol rate in symbols per second
    %   D       - Dispersion parameter in ps/(nm*km)
    %   L       - Length of the fiber in km
    %   lambda0 - Central wavelength of the signal in meters
    %
    % Output:
    %   Xpol_dc - Dispersed signal after traveling through the fiber
    %   Ypol_dc - Dispersed signal after traveling through the fiber

    % Calculate sampling frequency
    Fs = Sps * Rb;               % Sampling frequency in Hz

    % Constants
    c = 3e8;                     % Speed of light in vacuum (m/s)

    % Frequency vector setup
    N = length(Xpol);            % Number of points in the signal
    df = Fs / N;                 % Frequency resolution
    f = (-N/2:N/2-1)' * df;      % Frequency vector centered at zero

    % Convert D to s/m^2 for use with the frequency vector in Hz
    D_s = D * 1e-6;   % Convert D from ps/(nm*km) to s/m^2
    beta2 = -(lambda0^2 * D_s) / (2 * pi * c);  % Dispersion parameter in s^2/m

    % Dispersion phase shift
    phi = -pi * beta2 * L * 1e3 * f.^2;  % Dispersion phase factor in radians

    % Apply dispersion in frequency domain
    Xpol_f = fftshift(fft(Xpol));         % Transform to frequency domain
    Xpol_f_cd = Xpol_f .* exp(1i * phi);  % Apply dispersion
    Xpol_dc = ifft(ifftshift(Xpol_f_cd));  % Transform back to time domain
    Ypol_f = fftshift(fft(Ypol));         % Transform to frequency domain
    Ypol_f_cd = Ypol_f .* exp(1i * phi);  % Apply dispersion
    Ypol_dc = ifft(ifftshift(Ypol_f_cd));  % Transform back to time domain

    % Output the final dispersed signal
end
