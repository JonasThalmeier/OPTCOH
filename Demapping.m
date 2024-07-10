function [demappedBits_Xpol,demappedSymb_Xpol,demappedBits_Ypol, demappedSymb_Ypol] = Demapping(Xpol,Ypol,Sig,M)
dec = 0:M-1;
symb = zeros(length(dec),1);
bits = zeros(length(dec),log2(M));
for idx=1:length(dec)
    index = find(Sig.decSymbols  == dec(idx), 1);
    symb(idx) = Sig.txSymb(index);
    bits(idx,:) = Sig.bits(index,:);
end
demappedBits_Xpol = zeros(length(Xpol), log2(M));
demappedSymb_Xpol = zeros(length(Xpol),1);
for idx = 1:length(Xpol)
    [~,index] = min(abs(symb-Xpol(idx)));
    demappedSymb_Xpol(idx) = index-1;
    demappedBits_Xpol(idx,:) = bits(index,:);
end

demappedBits_Ypol = zeros(length(Ypol), log2(M));
demappedSymb_Ypol = zeros(length(Ypol),1);
for idx = 1:length(Ypol)
    [~,index] = min(abs(symb-Ypol(idx)));
    demappedSymb_Ypol(idx) = index-1;
    demappedBits_Ypol(idx,:) = bits(index,:);
end
end