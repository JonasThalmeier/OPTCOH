function [X_out, Y_out, e_X, e_Y] = LMS_N(TX_sig, modulation_r, mu, mu2, N_tap, training_sequence, M)

% for the moment consider as training sequence the whole transmission
% (65536*10)

% Reference vector for decoding
distances_vector = qammod(0:(M-1), M, 'UnitAveragePower', true);

%Taps initialization
h_xx = zeros(N_tap,1);
h_xy = zeros(N_tap,1);
h_yx = zeros(N_tap,1);
h_yy = zeros(N_tap,1);

h_xx(ceil(N_tap/2),1) = 1;
h_yy(ceil(N_tap/2),1) = 1;
h_xy(ceil(N_tap/2),1) = 0;
h_yx(ceil(N_tap/2),1) = 0;

% Align the received sequence with the reference one
% [~,transient_Xpol] = max(abs(xcorr(TX_sig(1:2:65536*2,1), training_sequence(1:65536,1))));
[transient_Xpol] = abs(finddelay(TX_sig(1:2:65536*2,1), training_sequence(1:65536,1)));
% TX_sig_1 = TX_sig(transient_Xpol+1:end,1);

% [~,transient_Ypol] = max(abs(xcorr(TX_sig(1:2:65536*2,2), training_sequence(1:65536,2))));
[transient_Ypol] = abs(finddelay(TX_sig(1:2:65536*2,2), training_sequence(1:65536,2)));
% TX_sig_2 = TX_sig(transient_Ypol+1:end,2);

% Check of consistency of the transients
if transient_Xpol ~= transient_Ypol
    fprintf('Transients in LMS are different\n');

    diff = abs(transient_Ypol-transient_Xpol);

    if transient_Ypol > transient_Xpol
        TX_sig_1 = TX_sig(transient_Xpol+1:end-diff,1);
        TX_sig_2 = TX_sig(transient_Ypol+1:end,2);
    else
        TX_sig_1 = TX_sig(transient_Xpol+1:end,1);
        TX_sig_2 = TX_sig(transient_Ypol+1:end-diff,2);
    end

else
        TX_sig_1 = TX_sig(transient_Xpol+1:end,1);
        TX_sig_2 = TX_sig(transient_Ypol+1:end,2);
end

TX_sig = [TX_sig_1, TX_sig_2];

% Variables initialization
e_X = zeros(1, floor(size(TX_sig,1)/2) - floor(N_tap/2) + 1);
e_Y = zeros(1, floor(size(TX_sig,2)/2) - floor(N_tap/2) + 1);

X_out = zeros(floor(size(TX_sig,1)/2) - floor(N_tap/2) + 1, 1);
Y_out = zeros(floor(size(TX_sig,2)/2) - floor(N_tap/2) + 1, 1);

% Cut training sequence to same length of RX sequence and starting at N_tap
training_sequence = training_sequence(floor(N_tap/2):floor(size(TX_sig,1)/2), :);

for i = N_tap:2:size(TX_sig,1)   %size(training_sequence,2) when we don't know the whole sequence

    k = floor(i/2) - floor(N_tap/2) + 1;

    x_in = TX_sig(i:-1:i-N_tap+1, 1);
    y_in = TX_sig(i:-1:i-N_tap+1, 2);
   
    X_out(k,1) = h_xx' * x_in + h_xy' * y_in;
    Y_out(k,1) = h_yx' * x_in + h_yy' * y_in;
        

    e_X(k)  = training_sequence(k, 1) - X_out(k,1);
    if isnan(e_X(k))
        fprintf('X IS NAN\n')
        pause;
    end
    e_Y(k)  = training_sequence(k, 2) - Y_out(k,1);
    if isnan(e_Y(k))
        fprintf('Y IS NAN\n')
        pause;
    end

    h_xx  = h_xx + mu * conj(e_X(k)) .* x_in;
    h_xy  = h_xy + mu * conj(e_X(k)) .* y_in;
    h_yx  = h_yx + mu * conj(e_Y(k)) .* x_in;
    h_yy  = h_yy + mu * conj(e_Y(k)) .* y_in;

end

figure(), plot(abs(fft(h_xx)));
    

% for i = len(transient_Xpol)+1:2:size(TX_sig,1)
%     
%     k = floor(i/2) - floor(N_tap/2) + 1;
% 
%     x_in = TX_sig(i:-1:i-N_tap+1, 1);
%     y_in = TX_sig(i:-1:i-N_tap+1, 2);
%    
%     X_out(k,1) = h_xx' * x_in + h_xy' * y_in;
%     Y_out(k,1) = h_yx' * x_in + h_yy' * y_in;
% 
%     for B_points = 1:B
%         [~, min_i]=min(abs(rotated_X_in(B_points,:) - QAM_points.').^2);
%         D(B_points,:) = QAM_points(min_i);
%     end
% 
%     e_X(k)  = RX_2 - rX;
%     if isnan(e_X(k))
%         fprintf('X IS NAN\n')
%         pause;
%     end
%     e_Y(k)  = RY_2 - rY;
%     if isnan(e_Y(k))
%         fprintf('Y IS NAN\n')
%         pause;
%     end
% 
%     h_xx  = h_xx + mu2 * conj(e_X(k)) .* x_in;
%     h_xy  = h_xy + mu2 * conj(e_X(k)) .* y_in;
%     h_yx  = h_yx + mu2 * conj(e_Y(k)) .* x_in;
%     h_yy  = h_yy + mu2 * conj(e_Y(k)) .* y_in;
% 
% end

