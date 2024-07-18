function [X_out, Y_out, e_X, e_Y] = CMA_RDE(TX_sig,r,mu,mu2,N_tap,N1,N2,GUIname,Rs)
% CMA_RDE performs equalization on the transmitted signal using CMA/RDE algorithms.
%
% Inputs:
%   TX_sig - Transmitted signal
%   r - Modulation index (1 for QPSK, 2 for 16QAM, 3 for 64QAM)
%   mu - Initial learning rate
%   mu2 - Learning rate after transition
%   N_tap - Number of taps for the FIR filter
%   N1 - Number of symbols when to change from CMA to RDE
%   N2 - Number of symbols when to change from mu to mu2
%   GUIname - Name of the GUI function to call if error occurs
%   Rs - Symbol rate in GBaud
%
% Outputs:
%   X_out - Equalized signal for X polarization
%   Y_out - Equalized signal for Y polarization
%   e_X - Error signal for X polarization
%   e_Y - Error signal for Y polarization

%Taps initialization
h_xx = zeros(N_tap,1);
h_xy = zeros(N_tap,1);
h_yx = zeros(N_tap,1);
h_yy = zeros(N_tap,1);

% Set the main (middle) tap to 1 for diagonal components
h_xx(ceil(N_tap/2),1) = 1;
h_yy(ceil(N_tap/2),1) = 1;
h_xy(ceil(N_tap/2),1) = 0;
h_yx(ceil(N_tap/2),1) = 0;

% Repeat the process 3 times
for rep = 1:3
    e_X = zeros(1, (size(TX_sig,2)-N_tap+1)/2);
    e_Y = zeros(1, (size(TX_sig,2)-N_tap+1)/2);

    X_out = zeros((size(TX_sig,2)-N_tap+1)/2, 1);
    Y_out = zeros((size(TX_sig,2)-N_tap+1)/2, 1);

    RDE_flag = 0;

    % Iterate through the transmitted signal with a step of 2
    for i = N_tap:2:size(TX_sig,1)

        k = floor(i/2) - floor(N_tap/2) + 1;

        % Extract input signal segments
        x_in = TX_sig(i:-1:i-N_tap+1, 1);
        y_in = TX_sig(i:-1:i-N_tap+1, 2);

        % Calculate the output signals
        X_out(k,1) = h_xx' * x_in + h_xy' * y_in;
        Y_out(k,1) = h_yx' * x_in + h_yy' * y_in;

        % Calculate the squared magnitudes of the output signals
        rX = abs(X_out(k,1))^2;
        rY = abs(Y_out(k,1))^2;

        % Transition to RDE if conditions are met
        if i>=N2 && RDE_flag==0 && rep==1
            RDE_flag = 1;
            mu = mu2;
        end

        % Adjust learning rate for subsequent repetitions
        if rep==2 && i==N_tap
            mu = 8e-5;
        elseif rep==3 && i==N_tap
            mu = 5e-5;
        end

        % Determine the desired modulus for CMA or RDE
        if r==1 || (RDE_flag==0 && i<N1 && rep==1)
            if r==1
                RX_2  = 1;
                RY_2  = 1;
            else
                RX_2  = 1.32;
                RY_2  = 1.32;
            end
        else
            if rX < 0.6
                RX_2 = 0.2;
            elseif rX > 1.4
                RX_2 = 1.8;
            else
                RX_2 = 1.0;
            end
            if rY < 0.6
                RY_2 = 0.2;
            elseif rY > 1.4
                RY_2 = 1.8;
            else
                RY_2 = 1.0;
            end
        end

        % Calculate the error signals
        e_X(k)  = RX_2 - rX;
        if isnan(e_X(k))
            % Display error message and restart GUI if error diverges
            errordlg('CMA/RDE error diverged, try different parameters');
            feval(GUIname,r,Rs,'CMA/RDE');
            error('CMA/RDE error diverged, try different parameters');
        end
        e_Y(k)  = RY_2 - rY;
        if isnan(e_Y(k))
            % Display error message and restart GUI if error diverges
            errordlg('CMA/RDE error diverged, try different parameters');
            feval(GUIname,r,Rs,'CMA/RDE');
            error('CMA/RDE error diverged, try different parameters');
        end
        % Update the filter taps
        h_xx  = h_xx + mu * e_X(k) .* x_in.*conj(X_out(k,1));
        h_xy  = h_xy + mu * e_X(k) .* y_in.*conj(X_out(k,1));
        h_yx  = h_yx + mu * e_Y(k) .* x_in.*conj(Y_out(k,1));
        h_yy  = h_yy + mu * e_Y(k) .* y_in.*conj(Y_out(k,1));


    end
end

