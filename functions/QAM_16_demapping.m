function [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QAM_16_demapping(downsampledSig_Xpol,downsampledSig_Ypol)
downsampledSig_Xpol = downsampledSig_Xpol(2:2:end); %downsample to 1SpS
downsampledSig_Ypol = downsampledSig_Ypol(2:2:end);
division = (abs(real(max(downsampledSig_Xpol)))+abs(real(min(downsampledSig_Xpol))))/2;

% Symbol Demapping
demappedBits_Xpol = zeros(length(downsampledSig_Xpol),4); % Adjust size accordingly for bit pairs, etc.
demappedSymb_Xpol = zeros(length(downsampledSig_Xpol),1);
demappedBits_Ypol = zeros(length(downsampledSig_Ypol),4); % Adjust size accordingly for bit pairs, etc.
demappedSymb_Ypol = zeros(length(downsampledSig_Ypol),1);

% Assuming QPSK and not accounting for noise, just a direct mapping
for i = 1:length(downsampledSig_Xpol)
    % This is a simplistic approach; real demapping would consider noise, etc.
    if real(downsampledSig_Xpol(i)) > 0
        if imag(downsampledSig_Xpol(i)) > 0
            if real(downsampledSig_Xpol(i)) > division
                if imag(downsampledSig_Xpol(i)) > division
                    demappedBits_Xpol(i, :) = [1 0 0 0];
                    demappedSymb_Xpol(i) = 8;
                else
                    demappedBits_Xpol(i, :) = [1 0 0 1];
                    demappedSymb_Xpol(i) = 9;
                end
            else
                if imag(downsampledSig_Xpol(i)) > division
                    demappedBits_Xpol(i, :) = [1 1 0 0];
                    demappedSymb_Xpol(i) = 12;
                else
                    demappedBits_Xpol(i, :) = [1 1 0 1];
                    demappedSymb_Xpol(i) = 13;
                end
            end
        else
            if real(downsampledSig_Xpol(i)) > division
                if imag(downsampledSig_Xpol(i)) > -division
                    demappedBits_Xpol(i, :) = [1 0 1 1];
                    demappedSymb_Xpol(i) = 11;
                else
                    demappedBits_Xpol(i, :) = [1 0 1 0];
                    demappedSymb_Xpol(i) = 10;
                end
            else
                if imag(downsampledSig_Xpol(i)) > -division
                    demappedBits_Xpol(i, :) = [1 1 1 1];
                    demappedSymb_Xpol(i) = 15;
                else
                    demappedBits_Xpol(i, :) = [1 1 1 0];
                    demappedSymb_Xpol(i) = 14;
                end
            end
        end
    else
         if imag(downsampledSig_Xpol(i)) > 0
            if real(downsampledSig_Xpol(i)) > -division
                if imag(downsampledSig_Xpol(i)) > division
                    demappedBits_Xpol(i, :) = [0 1 0 0];
                    demappedSymb_Xpol(i) = 4;
                else
                    demappedBits_Xpol(i, :) = [0 1 0 1];
                    demappedSymb_Xpol(i) = 5;
                end
            else
                if imag(downsampledSig_Xpol(i)) > division
                    demappedBits_Xpol(i, :) = [0 0 0 0];
                    demappedSymb_Xpol(i) = 0;
                else
                    demappedBits_Xpol(i, :) = [0 0 0 1];
                    demappedSymb_Xpol(i) = 1;
                end
            end
        else
            if real(downsampledSig_Xpol(i)) > -division
                if imag(downsampledSig_Xpol(i)) > -division
                    demappedBits_Xpol(i, :) = [0 1 1 1];
                    demappedSymb_Xpol(i) = 7;
                else
                    demappedBits_Xpol(i, :) = [0 1 1 0];
                    demappedSymb_Xpol(i) = 6;
                end
            else
                if imag(downsampledSig_Xpol(i)) > -division
                    demappedBits_Xpol(i, :) = [0 0 1 1];
                    demappedSymb_Xpol(i) = 3;
                else
                    demappedBits_Xpol(i, :) = [0 0 1 0];
                    demappedSymb_Xpol(i) = 2;
                end
            end
         end 
    end
end

for i = 1:length(downsampledSig_Ypol)
    % This is a simplistic approach; real demapping would consider noise, etc.
    if real(downsampledSig_Ypol(i)) > 0
        if imag(downsampledSig_Ypol(i)) > 0
            if real(downsampledSig_Ypol(i)) > division
                if imag(downsampledSig_Ypol(i)) > division
                    demappedBits_Ypol(i, :) = [1 0 0 0];
                    demappedSymb_Ypol(i) = 8;
                else
                    demappedBits_Ypol(i, :) = [1 0 0 1];
                    demappedSymb_Ypol(i) = 9;
                end
            else
                if imag(downsampledSig_Ypol(i)) > division
                    demappedBits_Ypol(i, :) = [1 1 0 0];
                    demappedSymb_Ypol(i) = 12;
                else
                    demappedBits_Ypol(i, :) = [1 1 0 1];
                    demappedSymb_Ypol(i) = 13;
                end
            end
        else
            if real(downsampledSig_Ypol(i)) > division
                if imag(downsampledSig_Ypol(i)) > -division
                    demappedBits_Ypol(i, :) = [1 0 1 1];
                    demappedSymb_Ypol(i) = 11;
                else
                    demappedBits_Ypol(i, :) = [1 0 1 0];
                    demappedSymb_Ypol(i) = 10;
                end
            else
                if imag(downsampledSig_Ypol(i)) > -division
                    demappedBits_Ypol(i, :) = [1 1 1 1];
                    demappedSymb_Ypol(i) = 15;
                else
                    demappedBits_Ypol(i, :) = [1 1 1 0];
                    demappedSymb_Ypol(i) = 14;
                end
            end
        end
    else
         if imag(downsampledSig_Ypol(i)) > 0
            if real(downsampledSig_Ypol(i)) > -division
                if imag(downsampledSig_Ypol(i)) > division
                    demappedBits_Ypol(i, :) = [0 1 0 0];
                    demappedSymb_Ypol(i) = 4;
                else
                    demappedBits_Ypol(i, :) = [0 1 0 1];
                    demappedSymb_Ypol(i) = 5;
                end
            else
                if imag(downsampledSig_Ypol(i)) > division
                    demappedBits_Ypol(i, :) = [0 0 0 0];
                    demappedSymb_Ypol(i) = 0;
                else
                    demappedBits_Ypol(i, :) = [0 0 0 1];
                    demappedSymb_Ypol(i) = 1;
                end
            end
        else
            if real(downsampledSig_Ypol(i)) > -division
                if imag(downsampledSig_Ypol(i)) > -division
                    demappedBits_Ypol(i, :) = [0 1 1 1];
                    demappedSymb_Ypol(i) = 7;
                else
                    demappedBits_Ypol(i, :) = [0 1 1 0];
                    demappedSymb_Ypol(i) = 6;
                end
            else
                if imag(downsampledSig_Ypol(i)) > -division
                    demappedBits_Ypol(i, :) = [0 0 1 1];
                    demappedSymb_Ypol(i) = 3;
                else
                    demappedBits_Ypol(i, :) = [0 0 1 0];
                    demappedSymb_Ypol(i) = 2;
                end
            end
         end 
    end
end
end