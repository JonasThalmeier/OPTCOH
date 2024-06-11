function [out] = EQ_func_1(Xpol_in,Ypol_in,mu,NTaps,alg,N1,N2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EQ_func performs equalization on a dual polarization signal using
% either the CMA or LMS algorithm.
%
% Inputs:
% - Xpol_in: Input signal for the X polarization (column vector)
% - Ypol_in: Input signal for the Y polarization (column vector)
% - mu: Step-size for coefficient update
% - NTaps: Number of taps for the filters
% - alg: Algorithm to be used ('CMA' or 'LMS')
% - Xorg: Original X polarization signal for training
% - Yorg: Original Y polarization signal for training
%
% Output:
% - out: Equalized signal (matrix with two columns, each corresponding
%        to one polarization orientation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants

R_CMA = sqrt(2);
% R_RDE = .25*(sqrt(2)+2*sqrt(10)+sqrt(18)).*[sqrt(2),sqrt(10),sqrt(18)];
R_RDE = sqrt(2).*[1/sqrt(5) 1 3/sqrt(5)];
SpS = 2;
NOut = 0;
% train_len = length(Xorg);

% Input preparation
x = [Xpol_in, Ypol_in];
x = [x(end-floor(NTaps/2)+1:end,:) ; x ; x(1:floor(NTaps/2),:)];
xV = convmtx(x(:,1).', NTaps);
xH = convmtx(x(:,2).', NTaps);
xV = xV(:,NTaps:SpS:end-NTaps+1);
xH = xH(:,NTaps:SpS:end-NTaps+1);

% Calculate output length
OutLength = floor((size(x,1) - NTaps + 1)/2);
clearvars x

% Initialize outputs
y1 = zeros(OutLength, 1);
y2 = zeros(OutLength, 1);

% Initialize filter coefficients

w1V = zeros(NTaps, 1);
w1H = zeros(NTaps, 1);
w2V = zeros(NTaps, 1);
w2H = zeros(NTaps, 1);
w1V(floor(NTaps/2) + 1) = 1;
% w2V(floor(NTaps/2) + 1) = 1;
% w2H(floor(NTaps/2) + 1) = 1;
% w1H(floor(NTaps/2) + 1) = 1;

error = zeros(1,OutLength);

for i = 1:OutLength
    % % Calculate the outputs
    y1(i) = w1V'*xV(:,i) + w1H'*xH(:,i);
    y2(i) = w2V'*xV(:,i) + w2H'*xH(:,i);
    % y1(i) = w1V' * xV(:, i);
    % y2(i) = w2H' * xH(:, i);


    if i == N1
        w2H = conj(w1V(end:-1:1,1)) ; w2V = -conj(w1H(end:-1:1,1));
        % w2H = conj(w1V) ; w2V = -conj(w1H);
        % w2H = conj(gram_schmidt(w1V,NTaps));
        % w2V = -conj(gram_schmidt(w1H,NTaps));
    end
    % Update filter coefficients
    if alg == "RDE" && i>N2
        [w1V, w1H, w2V, w2H, error(i)] = RDE(xV(:, i), xH(:, i), y1(i), y2(i), w1V, w1H, w2V, w2H, R_RDE, mu);
    elseif alg == "CMA" || alg == "RDE"
        [w1V, w1H, w2V, w2H, error(i)] = CMA(xV(:, i), xH(:, i), y1(i), y2(i), w1V, w1H, w2V, w2H, R_CMA, mu);
    elseif alg == "LMS"
        [w1V, w1H, w2V, w2H] = LMS(xV(:, i), xH(:, i), y1(i), y2(i), w1V, w1H, w2V, w2H, mu, Xorg(i), Yorg(i));
        if i >= train_len
            alg = "";
        end
    end
end
%
% temp=w1Vmat./max(w1Vmat,[],1);
% mesh(abs(w1Vmat(:,1:100:train_len)));

% Output samples
out = [y1 y2];
out = out(1 + NOut:end, :);
% figure,plot(error);
end

function orthogonal_vector = gram_schmidt(reference_vector, NTaps)
orthogonal_vector = zeros(size(reference_vector));
orthogonal_vector(floor(NTaps/2) + 1) = 1;  % Initialize with a unit vector

% Subtract the projection of the reference vector
projection = (dot(reference_vector, orthogonal_vector) / norm(reference_vector)^2) * reference_vector;
orthogonal_vector = orthogonal_vector - projection;

% Normalize the orthogonal vector
orthogonal_vector = orthogonal_vector / norm(orthogonal_vector);
end

function [w1V, w1H, w2V, w2H, error] = CMA(xV, xH, y1, y2, w1V, w1H, w2V, w2H, R, Mu)
% CMA algorithm for filter coefficient update
error = R - abs(y1).^2;
w1V = w1V + Mu * xV * (R - abs(y1).^2) * conj(y1);
w1H = w1H + Mu * xH * (R - abs(y1).^2) * conj(y1);
w2V = w2V + Mu * xV * (R - abs(y2).^2) * conj(y2);
w2H = w2H + Mu * xH * (R - abs(y2).^2) * conj(y2);
end

function [w1V, w1H, w2V, w2H,error] = RDE(xV, xH, y1, y2, w1V, w1H, w2V, w2H, R, Mu)
% Radius for output y1 and output y2:
[~,r1] = min(abs(R-abs(y1))) ; [~,r2] = min(abs(R-abs(y2)));
%Updating the filters:
error = R(r1)^2-abs(y1).^2;
w1V = w1V + Mu*xV*(R(r1)^2-abs(y1).^2)*conj(y1);
w1H = w1H + Mu*xH*(R(r1)^2-abs(y1).^2)*conj(y1);
w2V = w2V + Mu*xV*(R(r2)^2-abs(y2).^2)*conj(y2);
w2H = w2H + Mu*xH*(R(r2)^2-abs(y2).^2)*conj(y2);
end

function [w1V, w1H, w2V, w2H] = LMS(xV, xH, y1, y2, w1V, w1H, w2V, w2H, Mu, Xorg, Yorg)
% LMS algorithm for filter coefficient update
w1V = w1V + Mu * xV * y1 * conj(Xorg - y1);
w1H = w1H + Mu * xH * y1 * conj(Yorg - y1);
w2V = w2V + Mu * xV * y2 * conj(Xorg - y2);
w2H = w2H + Mu * xH * y2 * conj(Yorg - y2);
end