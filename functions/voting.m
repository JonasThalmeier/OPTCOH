function [consolidatedBits_Xpol,consolidatedBits_Ypol] = voting(N, M, demappedBits_Xpol, demappedBits_Ypol)
% Majority voting over Npp bits
% Let's assume demappedBits_Xpol is your bit outcomes with size [N * Npp, 2]
% Where N is the number of unique symbols and Npp is the number of repetitions

consolidatedBits_Xpol = zeros(N, M);
consolidatedBits_Ypol = zeros(N, M);

for i = 1:N
    consolidatedBits_Xpol(i, :) = mode(demappedBits_Xpol(i:N:end,:),1);
    consolidatedBits_Ypol(i, :) = mode(demappedBits_Ypol(i:N:end,:),1);
end
end

