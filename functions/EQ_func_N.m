function [X_out, Y_out, e_X, e_Y] = EQ_func_N(TX_sig,r,mu,mu2,N_tap,N1,N2)
% N1 when to change from CMA to RDE
% N2 when to go from mu to mu2

%Taps initialization
h_xx = zeros(N_tap,1);
h_xy = zeros(N_tap,1);
h_yx = zeros(N_tap,1);
h_yy = zeros(N_tap,1);

h_xx(ceil(N_tap/2),1) = 1;
h_yy(ceil(N_tap/2),1) = 1;
h_xy(ceil(N_tap/2),1) = 0;
h_yx(ceil(N_tap/2),1) = 0;

for rep = 1:3
    e_X = zeros(1, (size(TX_sig,2)-N_tap+1)/2);
    e_Y = zeros(1, (size(TX_sig,2)-N_tap+1)/2);

    X_out = zeros((size(TX_sig,2)-N_tap+1)/2, 1);
    Y_out = zeros((size(TX_sig,2)-N_tap+1)/2, 1);

    RDE_flag = 0;

    for i = N_tap:2:size(TX_sig,1)
        
        k = floor(i/2) - floor(N_tap/2) + 1;

        x_in = TX_sig(i:-1:i-N_tap+1, 1);
        y_in = TX_sig(i:-1:i-N_tap+1, 2);
       
        X_out(k,1) = h_xx' * x_in + h_xy' * y_in;
        Y_out(k,1) = h_yx' * x_in + h_yy' * y_in;

        rX = abs(X_out(k,1))^2;
        rY = abs(Y_out(k,1))^2;
        
        if i>=N2 && RDE_flag==0 && rep==1 
            RDE_flag = 1;
            mu = mu2;
        end

        if rep==2 && i==N_tap
            mu = 8e-5;
        elseif rep==3 && i==N_tap
            mu = 5e-5;
        end

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

       % fprintf('%d \n', RX_2)
        
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

        h_xx  = h_xx + mu * e_X(k) .* x_in.*conj(X_out(k,1));
        h_xy  = h_xy + mu * e_X(k) .* y_in.*conj(X_out(k,1));
        h_yx  = h_yx + mu * e_Y(k) .* x_in.*conj(Y_out(k,1));
        h_yy  = h_yy + mu * e_Y(k) .* y_in.*conj(Y_out(k,1));
    

    end
end

