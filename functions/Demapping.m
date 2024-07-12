function [demappedBits_Xpol, demappedSymb_Xpol, demappedBits_Ypol, demappedSymb_Ypol] = Demapping(Xpol, Ypol, Sig, M)
% Demapping Demaps the received symbols to bits for X and Y polarizations.
%
% Inputs:
%   Xpol - Received signal for X polarization
%   Ypol - Received signal for Y polarization
%   Sig - Structure containing the signal mapping information:
%       .decSymbols - Decimal representation of the transmitted symbols
%       .txSymb - Transmitted symbol values
%       .bits - Corresponding bits for each symbol
%   M - Modulation order (e.g., 16 for 16-QAM)
%
% Outputs:
%   demappedBits_Xpol - Demapped bits for X polarization
%   demappedSymb_Xpol - Demapped symbol indices for X polarization
%   demappedBits_Ypol - Demapped bits for Y polarization
%   demappedSymb_Ypol - Demapped symbol indices for Y polarization

% Initialize variables
dec = 0:M-1; % Decimal values for all possible symbols
symb = zeros(length(dec), 1); % Array to hold the transmitted symbols
bits = zeros(length(dec), log2(M)); % Array to hold the corresponding bits

% Populate symb and bits arrays with the transmitted symbol values and corresponding bits
for idx = 1:length(dec)
    index = find(Sig.decSymbols == dec(idx), 1); % Find the index of the current symbol
    symb(idx) = Sig.txSymb(index); % Store the transmitted symbol value
    bits(idx, :) = Sig.bits(index, :); % Store the corresponding bits
end

% Initialize output variables for X polarization
demappedBits_Xpol = zeros(length(Xpol), log2(M));
demappedSymb_Xpol = zeros(length(Xpol), 1);

% Demap the received symbols for X polarization
for idx = 1:length(Xpol)
    [~, index] = min(abs(symb - Xpol(idx))); % Find the closest transmitted symbol
    demappedSymb_Xpol(idx) = index - 1; % Store the symbol index
    demappedBits_Xpol(idx, :) = bits(index, :); % Store the corresponding bits
end

% Initialize output variables for Y polarization
demappedBits_Ypol = zeros(length(Ypol), log2(M));
demappedSymb_Ypol = zeros(length(Ypol), 1);

% Demap the received symbols for Y polarization
for idx = 1:length(Ypol)
    [~, index] = min(abs(symb - Ypol(idx))); % Find the closest transmitted symbol
    demappedSymb_Ypol(idx) = index - 1; % Store the symbol index
    demappedBits_Ypol(idx, :) = bits(index, :); % Store the corresponding bits
end
end