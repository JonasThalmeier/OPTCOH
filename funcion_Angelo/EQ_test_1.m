close all;
clear all;
clc;

N = 1e6;
seq1 = pskmod(randi([0 3],N,1),4);
seq2 = pskmod(randi([0 3],N,1),4);
seq1_up=upsample(seq1,2);
seq2_up=upsample(seq2,2);
seq1_up(2:2:end)=seq1_up(1:2:end);
seq2_up(2:2:end)=seq2_up(1:2:end);

rng(42);
h = randi([0 1e3],5,1)./1e3;
h = h./sum(h);

dist1 = conv(seq1_up,h);
dist2 = conv(seq2_up,h);
dist1 = dist1./mean(abs(dist1));
dist2 = dist2./mean(abs(dist2));

% mu=1e-3.25 seems to work fine for CMA
for idx = -5:.5:-2
    test=EQ_func(dist1,dist2,10^idx,7,"CMA",seq1(1:round(3*N/4)),seq2(1:round(3*N/4)));
    scatterplot(test(round(N/2):end,1)); title(sprintf('Output EQ with mu=10^{%d}',idx));
end
scatterplot(dist1(round(N/2):end)); title("EQ input");