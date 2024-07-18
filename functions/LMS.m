function [X_out, Y_out, e_X, e_Y] = LMS(Xpol, Ypol, mu, mu2, N_tap, training_sequence, M, train_length,GUIname,r,Rs)
% LMS Performs LMS equalization for M-QAM modulated signals.
%
% Inputs:
%   Xpol - Input signal for X polarization
%   Ypol - Input signal for Y polarization
%   mu - Learning rate for the training phase
%   mu2 - Learning rate for the tracking phase
%   N_tap - Number of taps for the adaptive filter
%   training_sequence - Known training sequence for the equalizer training
%   M - Modulation order (e.g., 16 for 16-QAM)
%   train_length - Number of symbols used for training
%
% Outputs:
%   X_out - Equalized output signal for X polarization
%   Y_out - Equalized output signal for Y polarization
%   e_X - Error signal for X polarization
%   e_Y - Error signal for Y polarization

% Concatenate the input signals
TX_sig = [Xpol, Ypol];

% Reference vector for decoding (QAM constellation points)
distances_vector = qammod(0:(M-1), M);

% Taps initialization
h_xx = zeros(N_tap, 1);
h_xy = zeros(N_tap, 1);
h_yx = zeros(N_tap, 1);
h_yy = zeros(N_tap, 1);

% Set the main (middle) tap to 1 for the initial filter coefficients
h_xx(ceil(N_tap/2), 1) = 1;
h_yy(ceil(N_tap/2), 1) = 1;
h_xy(ceil(N_tap/2), 1) = 0;
h_yx(ceil(N_tap/2), 1) = 0;

% Align the received sequence with the reference one using finddelay
[transient_Xpol] = 2 * abs(finddelay(TX_sig(1:2:65536*2, 1), training_sequence(1:65536, 1)));
[transient_Ypol] = 2 * abs(finddelay(TX_sig(1:2:65536*2, 2), training_sequence(1:65536, 2)));

% Check consistency of the transients and align the signals accordingly
if transient_Xpol ~= transient_Ypol
    fprintf('Transients in LMS are different\n');
    diff = abs(transient_Ypol - transient_Xpol);
    if transient_Ypol > transient_Xpol
        TX_sig_1 = TX_sig(transient_Xpol + 1:end - diff, 1);
        TX_sig_2 = TX_sig(transient_Ypol + 1:end, 2);
    else
        TX_sig_1 = TX_sig(transient_Xpol + 1:end, 1);
        TX_sig_2 = TX_sig(transient_Ypol + 1:end - diff, 2);
    end
else
    TX_sig_1 = TX_sig(transient_Xpol + 1:end, 1);
    TX_sig_2 = TX_sig(transient_Ypol + 1:end, 2);
end

% Zero padding to align the first symbol with the main tap of the filter
TX_sig = [zeros(floor(N_tap/2), 2); TX_sig_1, TX_sig_2; zeros(floor(N_tap/2), 2)];

% Variables initialization
e_X = zeros(1, floor(size(TX_sig, 1)/2) - floor(N_tap/2) + 1);
e_Y = zeros(1, floor(size(TX_sig, 2)/2) - floor(N_tap/2) + 1);

X_out = zeros(floor(size(TX_sig, 1)/2) - floor(N_tap/2) + 1, 1);
Y_out = zeros(floor(size(TX_sig, 2)/2) - floor(N_tap/2) + 1, 1);

%% TRAINING PHASE
for i = N_tap:2:2 * train_length
    k = floor(i/2) - floor(N_tap/2) + 1;
    x_in = TX_sig(i:-1:i-N_tap+1, 1);
    y_in = TX_sig(i:-1:i-N_tap+1, 2);
    X_out(k, 1) = h_xx' * x_in + h_xy' * y_in;
    Y_out(k, 1) = h_yx' * x_in + h_yy' * y_in;
    e_X(k) = training_sequence(k, 1) - X_out(k, 1);
    if isnan(e_X(k))
        errordlg('LMS error diverged in training, try different parameters');
        feval(GUIname,r,Rs,'LMS');
        error('LMS error diverged in training, try different parameters');
    end
    e_Y(k) = training_sequence(k, 2) - Y_out(k, 1);
    if isnan(e_Y(k))
        errordlg('LMS error diverged in training, try different parameters');
        feval(GUIname,r,Rs,'LMS');
        error('LMS error diverged in training, try different parameters');
    end
    % Update the filter coefficients using LMS algorithm
    h_xx = h_xx + mu * conj(e_X(k)) .* x_in;
    h_xy = h_xy + mu * conj(e_X(k)) .* y_in;
    h_yx = h_yx + mu * conj(e_Y(k)) .* x_in;
    h_yy = h_yy + mu * conj(e_Y(k)) .* y_in;
end

%% TRACKING PHASE
for j = i + 2:2:size(TX_sig, 1)
    k = floor(j/2) - floor(N_tap/2) + 1;
    x_in = TX_sig(j:-1:j-N_tap+1, 1);
    y_in = TX_sig(j:-1:j-N_tap+1, 2);
    X_out(k, 1) = h_xx' * x_in + h_xy' * y_in;
    Y_out(k, 1) = h_yx' * x_in + h_yy' * y_in;
    % Estimate the transmitted symbols using minimum distance criterion
    [~, ind] = min(abs(distances_vector - X_out(k, 1)));
    X_est = distances_vector(ind);
    [~, ind] = min(abs(distances_vector - Y_out(k, 1)));
    Y_est = distances_vector(ind);
    e_X(k) = X_est - X_out(k, 1);
    if isnan(e_X(k))
        errordlg('LMS error diverged in tracking, try different parameters');
        feval(GUIname,r,Rs,'LMS');
        error('LMS error diverged in tracking, try different parameters');
    end
    e_Y(k) = Y_est - Y_out(k, 1);
    if isnan(e_Y(k))
        errordlg('LMS error diverged in tracking, try different parameters');
        feval(GUIname,r,Rs,'LMS');
        error('LMS error diverged in tracking, try different parameters');
    end
    % Update the filter coefficients using LMS algorithm
    h_xx = h_xx + mu2 * conj(e_X(k)) .* x_in;
    h_xy = h_xy + mu2 * conj(e_X(k)) .* y_in;
    h_yx = h_yx + mu2 * conj(e_Y(k)) .* x_in;
    h_yy = h_yy + mu2 * conj(e_Y(k)) .* y_in;
end
% Uncomment the following line to plot the FFT of the h_xx filter coefficients
% figure(), plot(abs(fft(h_xx))); title('FFT(h_xx)');
end