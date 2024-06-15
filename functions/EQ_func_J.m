function [X_out, Y_out, e_X, e_Y] = EQ_func_J(TX_sig,r,mu,mu2,N_tap,N1,N2)
% N1 when to change from CMA to RDE
% N2 when to go from mu to mu2

%Taps initialization
h_xx = zeros(N_tap,1);
h_xy = zeros(N_tap,1);
h_yx = zeros(N_tap,1);
h_yy = zeros(N_tap,1);

if r==1
    h_xx(1,1) = cos(pi/6);
    h_yy(1,1) = cos(pi/6);
    h_xy(1,1) = exp(-1i*pi/5)*sin(pi/6);
    h_yx(1,1) = -exp(1i*pi/5)*sin(pi/6);
else
    h_xx(ceil(N_tap/2),1) = 1;
    h_yy(ceil(N_tap/2),1) = 1;
    h_xy(ceil(N_tap/2),1) = 0;
    h_yx(ceil(N_tap/2),1) = 0;
end

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

    if i==N1
        % h_yy = conj(h_xx(end:-1:1,1));
        % h_yx = -conj(h_xy(end:-1:1,1));
        if r==2
            RDE_flag = 1;
        end
    end
    if i==N2
        mu = mu2;
    end

    if RDE_flag==0
    % CMA
        if r==1
            RX_2  = sqrt(2);
            RY_2  = sqrt(2);
        else
            RX_2  = sqrt(2);
            RY_2  = sqrt(2);
        end
    else
    % RDE
        if rX < 1.0233
            RX_2 = sqrt(2/5);
        elseif rX > 1.6558
            RX_2 = 3*sqrt(2/5);
        else
            RX_2 = sqrt(2);
        end
        if rY < 1.0233
            RY_2 = sqrt(2/5);
        elseif rY > 1.6558
            RY_2 = 3*sqrt(2/5);
        else
            RY_2 = sqrt(2);
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

    h_xx  = h_xx + mu * e_X(k) .* x_in.*conj(X_out(k,1));
    h_xy  = h_xy + mu * e_X(k) .* y_in.*conj(X_out(k,1));
    h_yx  = h_yx + mu * e_Y(k) .* x_in.*conj(Y_out(k,1));
    h_yy  = h_yy + mu * e_Y(k) .* y_in.*conj(Y_out(k,1));

end

