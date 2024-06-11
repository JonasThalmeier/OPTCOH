function [X_out, Y_out, e_X, e_Y] = EQ_func_N(TX_sig,r,mu,mu2,N_tap,N1,N2)
% N1 when to change from CMA to RDE
% N2 when to go from mu to mu2

%Taps initialization
h_xx = zeros(1,N_tap);
h_xy = zeros(1,N_tap);
h_yx = zeros(1,N_tap);
h_yy = zeros(1,N_tap);

if r==1
    h_xx(ceil(N_tap/2)) = cos(pi/6);
    h_yy(ceil(N_tap/2)) = cos(pi/6);
    h_xy(ceil(N_tap/2)) = exp(-1i*pi/5)*sin(pi/6);
    h_yx(ceil(N_tap/2)) = -exp(1i*pi/5)*sin(pi/6);
else
    h_xx(ceil(N_tap/2)) = exp(1i*pi/6)*cos(pi/6);
    h_yy(ceil(N_tap/2)) = exp(1i*pi/6)*cos(pi/6);
    h_xy(ceil(N_tap/2)) = exp(-1i*pi/6)*sin(pi/6);
    h_yx(ceil(N_tap/2)) = -exp(1i*pi/6)*sin(pi/6);
end

for rep = 1:3
    e_X = zeros(1, size(TX_sig,2)/2);
    e_Y = zeros(1, size(TX_sig,2)/2);

    X_out = zeros(1,size(TX_sig,2)/2);
    Y_out = zeros(1,size(TX_sig,2)/2);

    out_index = 0;
    RDE_flag = 0;

    for i = 1:2:size(TX_sig,2)

        if (i+N_tap-1>=size(TX_sig,2))
            break;
        else
            out_index = out_index + 1;
            X_out(out_index) = sum(conj(h_xx).*fliplr(TX_sig(1,i:i+N_tap-1))) + sum(conj(h_xy).*fliplr(TX_sig(2,i:i+N_tap-1)));
            Y_out(out_index) = sum(conj(h_yx).*fliplr(TX_sig(1,i:i+N_tap-1))) + sum(conj(h_yy).*fliplr(TX_sig(2,i:i+N_tap-1)));

            rX = abs(X_out(out_index))^2;
            rY = abs(Y_out(out_index))^2;
            
            if out_index==N2 && RDE_flag==0 && rep==1
                RDE_flag = 1;
                mu = mu2;
            end

            if rep==2 && i==1
                mu = 8e-5;
            elseif rep==3 && i==1
                mu = 5e-5;
            end

            if r==1 || (RDE_flag==0 && out_index<N1 && rep==1)
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

            e_X(out_index)  = RX_2 - rX;
            if isnan(e_X(out_index))
                fprintf('X IS NAN\n')
                pause;
            end
            e_Y(out_index)  = RY_2 - rY;
            if isnan(e_Y(out_index))
                fprintf('Y IS NAN\n')
                pause;
            end

            h_xx  = h_xx + mu * e_X(out_index) .* fliplr(TX_sig(1,i:i+N_tap-1)).*conj(X_out(out_index));
            h_xy  = h_xy + mu * e_X(out_index) .* fliplr(TX_sig(2,i:i+N_tap-1)).*conj(X_out(out_index));
            h_yx  = h_yx + mu * e_Y(out_index) .* fliplr(TX_sig(1,i:i+N_tap-1)).*conj(Y_out(out_index));
            h_yy  = h_yy + mu * e_Y(out_index) .* fliplr(TX_sig(2,i:i+N_tap-1)).*conj(Y_out(out_index));
        end

    end
end

