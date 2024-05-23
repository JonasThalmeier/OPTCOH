close all;
clc;

N = 1e6;
seq1 = pskmod(randi([0 3],N,1),4);
seq2 = pskmod(randi([0 3],N,1),4);

rng(42);
h = randi([0 1e3],5,1)./1e3;
h = h./sum(h);

dist1 = conv(seq1,h);
dist2 = conv(seq2,h);
dist1 = dist1./mean(abs(dist1));
dist2 = dist2./mean(abs(dist2));

% mu=1e-3.25 seems to work fine for CMA
for idx = -5:.5:-2
    test=EQ_func(dist1,dist2,10^idx,7,"LMS",seq1(1:round(3*N/4)),seq2(1:round(3*N/4)));
    scatterplot(test(round(N/2):end,1)); title(sprintf('Output EQ with mu=10^{%d}',idx));
end
scatterplot(dist1(round(N/2):end)); title("EQ input");
