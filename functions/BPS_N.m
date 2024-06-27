function [X_rec, theta] = BPS_N(X_in, B, M, power_norm)

    L_x = length(X_in);
    theta = zeros(L_x, 1);
    X_rec = zeros(size(X_in));
    D = zeros(B,L_x);
    N = 50;
    a = zeros(1,L_x);
    phi_k = zeros(1,2);

    % Normalize input signal power
    X_Power = mean(abs(X_in).^2);
    X_in = X_in / sqrt(X_Power / power_norm);

    % Phase candidates
    phi = (0:B-1) / B * pi/2;

    % Precompute rotated symbols
    rotated_X_in = exp(-1i * phi.') * X_in.';

    QAM_points=[1+1i 1-1i -1+1i -1-1i 1+3i 3+1i 3-1i 1-3i -1-3i -3-1i -3+1i -1+3i 3+3i 3-3i -3-3i -3+3i]; 


    for B_points = 1:B
        [~, min_i]=min(abs(rotated_X_in(B_points,:) - QAM_points.').^2);
        D(B_points,:) = QAM_points(min_i);
    end
        

    for k = 1:L_x

        % Determine window boundaries
        kmN = max(1, k-N);
        kpN = min(L_x, k+N);

        % Reshape rotated symbols within the window
        X_window = reshape(rotated_X_in(:, kmN:kpN), B, []);

        % Calculate distances to QAM-16 constellation points
%         D = QAM_16_demapping_BPS(reshape(X_window, [], 1));
%         D = reshape(D, B, []);

        % Compute MSE for each phase candidate
        mse = sum(abs(X_window - D(:, kmN:kpN)).^2, 2);

        % Find minimum MSE and corresponding phase
        [~, min_index] = min(mse);
        phi_k(2) = phi(min_index);
        

%         if k-1>0
%             a(k) = a(k-1) + floor(1/2 + 1/(2*pi)*(phi_k(2) - phi_k(1))); %everything is moved by one to have 0 at a(1)
%         else
%             a(k) = floor(1/2 + 1/(2*pi)*(phi_k(2)));
%         end
% 
%         phi_k(1) = phi_k(2);
% 
%         theta(k) = phi_k(2) + a(k);
        % Correct the phase
        theta(k) = phi_k(2); % Optionally add unwrapping logic here
        
        % Recover the symbol with corrected phase
        X_rec(k) = X_in(k) * exp(-1i * theta(k));

        % Progress logging
        if k == floor(L_x/4)
            fprintf('BPS 25/100\n');
        elseif k == floor(L_x/2)
            fprintf('BPS 50/100\n');
        elseif k == floor(L_x*3/4)
            fprintf('BPS 75/100\n');
        end
    end
end