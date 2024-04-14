function [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = QPSK_demapping(downsampledSig_Xpol,downsampledSig_Ypol)
downsampledSig_Xpol = downsampledSig_Xpol(2:2:end); %downsample to 1SpS
downsampledSig_Ypol = downsampledSig_Ypol(2:2:end);

% Symbol Demapping
demappedBits_Xpol = zeros(length(downsampledSig_Xpol),2); % Adjust size accordingly for bit pairs, etc.
demappedSymb_Xpol = zeros(length(downsampledSig_Xpol),1);
demappedBits_Ypol = zeros(length(downsampledSig_Ypol),2); % Adjust size accordingly for bit pairs, etc.
demappedSymb_Ypol = zeros(length(downsampledSig_Ypol),1);

% Assuming QPSK and not accounting for noise, just a direct mapping
for i = 1:length(downsampledSig_Xpol)
    % This is a simplistic approach; real demapping would consider noise, etc.
    if real(downsampledSig_Xpol(i)) > 0
        if imag(downsampledSig_Xpol(i)) > 0
            demappedBits_Xpol(i, :) = [1 0];
            demappedSymb_Xpol(i) = 2;
        else
            demappedBits_Xpol(i, :) = [1 1];
            demappedSymb_Xpol(i) = 3;
        end
    else
        if imag(downsampledSig_Xpol(i)) > 0
            demappedBits_Xpol(i, :) = [0 0];
            demappedSymb_Xpol(i) = 0;
        else
            demappedBits_Xpol(i, :) = [0 1];
            demappedSymb_Xpol(i) = 1;
        end
    end
end

for i = 1:length(downsampledSig_Ypol)
    % This is a simplistic approach; real demapping would consider noise, etc.
    if real(downsampledSig_Ypol(i)) > 0
        if imag(downsampledSig_Ypol(i)) > 0
            demappedBits_Ypol(i, :) = [1 0];
            demappedSymb_Ypol(i) = 2;
        else
            demappedBits_Ypol(i, :) = [1 1];
            demappedSymb_Ypol(i) = 3;
        end
    else
        if imag(downsampledSig_Ypol(i)) > 0
            demappedBits_Ypol(i, :) = [0 0];
            demappedSymb_Ypol(i) = 0;
        else
            demappedBits_Ypol(i, :) = [0 1];
            demappedSymb_Ypol(i) = 1;
        end
    end
end
end