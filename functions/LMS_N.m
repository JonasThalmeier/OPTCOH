function [X_out, Y_out, e_X, e_Y] = LMS_N(TX_sig, modulation_r, mu, mu2, N_tap, training_sequence, M)

% for the moment consider as training sequence the whole transmission
% (65536*10)

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


e_X = zeros(1, (size(TX_sig,2)-N_tap+1)/2);
e_Y = zeros(1, (size(TX_sig,2)-N_tap+1)/2);

X_out = zeros((size(TX_sig,2)-N_tap+1)/2, 1);
Y_out = zeros((size(TX_sig,2)-N_tap+1)/2, 1);

[~,transient_Xpol] = max(abs(xcorr(TX_sig(1:65536*2,1), training_sequence(1:65536,1))));
TX_sig = TX_sig(transient_Xpol+1:end,1);

[~,transient_Ypol] = max(abs(xcorr(TX_sig(1:65536*2,2), training_sequence(1:65536,2))));
TX_sig = TX_sig(transient_Ypol+1:end,2);

if transient_Xpol ~= transient_Ypol
    fprintf('Transients in LMS are different\n');
end

training_sequence = training_sequence(1:size(TX_sig,2), :);

for i = N_tap:2:size(TX_sig,2)   %size(training_sequence,2) when we don't know the whole sequence

    k = floor(i/2) - floor(N_tap/2) + 1;

    x_in = TX_sig(i:-1:i-N_tap+1, 1);
    y_in = TX_sig(i:-1:i-N_tap+1, 2);
   
    X_out(k,1) = h_xx' * x_in + h_xy' * y_in;
    Y_out(k,1) = h_yx' * x_in + h_yy' * y_in;

    carrSynch = comm.CarrierSynchronizer("Modulation", modulation_r, "SamplesPerSymbol", 1,'DampingFactor', 31.6);
    [X_eq, ~] = carrSynch(X_out);
    
    carrSynch2 = comm.CarrierSynchronizer("Modulation", modulation_r, "SamplesPerSymbol", 1,'DampingFactor', 31.6);
    [Y_eq, ~] = carrSynch2(Y_out);

    error_tempX = zeros(1,4);
    error_tempY = zeros(1,4);

    for ind=0:3
                  
        X_RX = X_eq*exp(1i*ind*pi/2);
        Y_RX = Y_eq*exp(1i*ind*pi/2);

        x_demod = zeros(size(X_RX,1), size(X_RX,2));
        y_demod = zeros(size(Y_RX,1), size(Y_RX,2));

        for points = 1:length(X_RX)
            [~, min_x]=min(abs(X_RX(points,:) - distances_vector.').^2);
            x_demod(B_points,:) = distances_vector(min_x);

            [~, min_y]=min(abs(Y_RX(points,:) - distances_vector.').^2);
            y_demod(points,:) = distances_vector(min_y);
        end


        
    end

    e_X(k)  = RX_2 - rX;
    if isnan(e_X(k))
        fprintf('X IS NAN\n')
        pause;
    end
    e_Y(k)  = RY_2 - rY;
    if isnan(e_Y(k))
        fprintf('Y IS NAN\n')
        pause;
    end

    h_xx  = h_xx + mu * conj(e_X(k)) .* x_in;
    h_xy  = h_xy + mu * conj(e_X(k)) .* y_in;
    h_yx  = h_yx + mu * conj(e_Y(k)) .* x_in;
    h_yy  = h_yy + mu * conj(e_Y(k)) .* y_in;

end
    

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

