function [X_rec,theta] = BPS_N(X_in, B, M, power_norm)

    L_x = length(X_in);
    a = zeros(1,L_x+1);
    phi_k = zeros(L_x,1);
    theta = zeros(L_x,1);
    X_rec = zeros(size(X_in));
    N=10;

    X_Power = mean(abs((X_in)).^2);
    X_in = X_in/sqrt(X_Power/power_norm);

    b = 0:B-1;
    phi = b/B.*pi/2;
    
    for k = 1:L_x

        if(k-N <= 0)
            kmN = 1;
        else
            kmN = k-N;
        end

         if(k+N > L_x)
            kpN = L_x;
        else
            kpN = k+N;
        end

        
       D = reshape(QAM_16_demapping_BPS(reshape(exp(-1i*phi).* X_in(kmN:kpN), length(X_in(kmN:kpN))*B,1)), length(X_in(kmN:kpN)),B);
      
        
        [~, phi_k(k)] = min(sum(abs( exp(-1i*phi).* X_in(kmN:kpN) - D).^2, 1));
        
%         if k-1>0
%             a(k+1) = a(k) + floor(1/2 + M/(2*pi)*(phi(phi_k(k))-phi(phi_k(k-1)))); %everything is moved by one to have 0 at a(1)
%         else
%             a(k+1) = a(k) + floor(1/2 + M/(2*pi)*(phi(phi_k(k))));
%         end

        
        theta(k) = phi(phi_k(k));% + a(k+1);
        
        X_rec(k) = X_in(k)*exp(-1i*theta(k));

        if (k==floor(L_x/4))
            fprintf('BPS 25/100\n');
        elseif (k==floor(L_x/2))
            fprintf('BPS 50/100\n');
        elseif((k==floor(L_x/4*3)))
            fprintf('BPS 75/100\n');
        end
        
    end

end