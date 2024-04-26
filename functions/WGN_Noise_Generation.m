function [delay_phase_noise_distorted_RX_Apol, Noise] = WGN_Noise_Generation(delay_phase_distorted_RX_Apol,SpS, M, EbN0)

Nbit = log2(M);
%EbN0 = 30; %dB
EbN0_lin = 10^(EbN0/10);
Nsamples = length(delay_phase_distorted_RX_Apol);
% Signal Power evaluation
SignalPower=mean(abs(delay_phase_distorted_RX_Apol).^2);
% Noise variance evaluation
NoiseVariance=SignalPower/EbN0_lin*SpS/Nbit;
% Normalized Complex WGN generation (Noise Power set to 1)
NoiseNormalized=(randn(1,Nsamples)+1i*randn(1,Nsamples))/sqrt(2);
% Noise scaling to the given Eb/N0
Noise=NoiseNormalized*sqrt(NoiseVariance); %1e-3 pratically no noise it works, 1e-2, 1e-1 too

fprintf('The power noise generated is %.4f dBm\n ', 10*log10((mean(abs(Noise).^2))));

delay_phase_noise_distorted_RX_Apol = delay_phase_distorted_RX_Apol + Noise';

end