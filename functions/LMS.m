function [X_out2,e_xx] = LMS(TX_sig,Pilot,mu,N_tap)
% N1 when to change from CMA to RDE
% N2 when to go from mu to mu2
% TX_sig = TX_sig/mean(abs(TX_sig));
%Taps initialization
leng = 65536;
%% CMA RDE
h_xx = zeros(N_tap,1);
h_xy = zeros(N_tap,1);
h_yx = zeros(N_tap,1);
h_yy = zeros(N_tap,1);

h_xx(ceil(N_tap/2),1) = 1;
h_yy(ceil(N_tap/2),1) = 1;
h_xy(ceil(N_tap/2),1) = 0;
h_yx(ceil(N_tap/2),1) = 0;
TX_sig = [zeros(floor(N_tap/2),2);TX_sig;zeros(floor(N_tap/2),2)];
X_out = zeros((size(TX_sig,2)-N_tap+1)/2, 2);
for i = N_tap:2:size(TX_sig,1)
    k = floor(i/2) - floor(N_tap/2) + 1;
    x_in = TX_sig(i-N_tap+1:1:i, 1);
    y_in = TX_sig(i-N_tap+1:1:i, 2);

    X_out(k,1) = h_xx' * x_in + h_xy' * y_in;
    X_out(k,2) = h_yx' * x_in + h_yy' * y_in;
    rX = abs(X_out(k,1));
    rY = abs(X_out(k,2));
    if rX > 1.236
        R_X = 1.416;
    elseif rX > 0.7641
        R_X = 1.056;
    else
        R_X = 0.4721;
    end

    if rY > 1.236
        R_Y = 1.416;
    elseif rX > 0.7641
        R_Y = 1.056;
    else
        R_Y = 0.4721;
    end

    h_xx  = h_xx + mu * (R_X^2-rX^2) .* conj(X_out(k,1)) .* x_in;
    h_yx  = h_yx + mu * (R_Y^2-rY^2) .* conj(X_out(k,2)) .* x_in;
    h_yy  = h_yy + mu * (R_Y^2-rY^2) .* conj(X_out(k,2)) .* y_in;
    h_xy  = h_xy + mu * (R_X^2-rX^2) .* conj(X_out(k,1)) .* y_in;
end

scatterplot(X_out(round(end/2):end,1));title("After CMA RDE");
% figure;stem(abs(X_out));title("After RDE");
%% LMS
delay = finddelay(Pilot(1:65000,1),X_out(end-2*(leng-floor(Ntap/2)):end,1));
TX_sig = TX_sig(2*delay+floor(Ntap/2)+1:end,:);
X_out2 = zeros((size(TX_sig,2)-N_tap+1)/2, 2);
e_xx = zeros((size(TX_sig,2)-N_tap+1)/2, 2);
N_tap = N_tap;
mu = mu;
h_xx = zeros(N_tap,1);
h_xy = zeros(N_tap,1);
h_yx = zeros(N_tap,1);
h_yy = zeros(N_tap,1);

h_xx(ceil(N_tap/2),1) = 1;
h_yy(ceil(N_tap/2),1) = 1;
h_xy(ceil(N_tap/2),1) = 0;
h_yx(ceil(N_tap/2),1) = 0;

for i = N_tap:2:size(TX_sig,1)

    k = floor(i/2) - floor(N_tap/2) + 1;
    x_in = TX_sig(i-N_tap+1:1:i, 1);
    y_in = TX_sig(i-N_tap+1:1:i, 2);
    if k==3e5
        mu = mu/2;
        figure;stem(abs(h_xx));
    elseif k==6e5
        mu=mu/2;
        figure;stem(abs(h_xx));
    end
    X_out2(k,1) = h_xx' * x_in + h_xy' * y_in;
    X_out2(k,2) = h_yx' * x_in + h_yy' * y_in;
    e(k,:) = (Pilot(k,:)-X_out2(k,:));
    h_xx  = h_xx + mu * conj(e(k,1)) .* (x_in);
    h_yx  = h_yx + mu * conj(e(k,2)) .* (x_in);
    h_yy  = h_yy + mu * conj(e(k,2)) .* (y_in);
    h_xy  = h_xy + mu * conj(e(k,1)) .* (y_in);
end
figure;stem(abs(h_xx));
end

