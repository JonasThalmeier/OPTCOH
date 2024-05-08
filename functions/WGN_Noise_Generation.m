function [delay_phase_noise_distorted_RX_Apol, Noise] = WGN_Noise_Generation(delay_phase_distorted_RX_Apol,SpS, M, EbN0, Rs)

Nbit = log2(M);
%EbN0 = 30; %dB
EbN0_lin = 10^(EbN0/10);
Nsamples = length(delay_phase_distorted_RX_Apol);
% Signal Power evaluation
SignalPower=mean(abs(delay_phase_distorted_RX_Apol).^2);
% Noise variance evaluation (E_s = SignaPower*SpS; E_b = E_s/Nbit;
% N_0=E_b/SNR; sigma^2 =N_0/2)
NoisePower = SignalPower*SpS/(Nbit*EbN0_lin);
% Normalized Complex WGN generation (Noise Power set to 1)
% NoiseNormalized=(randn(1,Nsamples)+1i*randn(1,Nsamples)/sqrt(2));
NoiseNormalized=(randn(1,Nsamples)+1i*randn(1,Nsamples));
% Noise scaling to the given Eb/N0
Noise=NoiseNormalized*sqrt(NoisePower/4); %1e-3 pratically no noise it works, 1e-2, 1e-1 too

% sigma = 1/sqrt(2)*10^(-snr/20);

% Nsamples = length(delay_phase_distorted_RX_Apol);
% EbN0_lin = 10^(EbN0/10);
% SignalPower=mean(abs(delay_phase_distorted_RX_Apol).^2);
% Bsim = SpS * Rs;
% N0 = SignalPower / (2*EbN0_lin*Rs);
% Tot_power_noise = 1/2*N0*Bsim;
% noise_power_perquad_perpol = Tot_power_noise/4;
% sigma = sqrt(noise_power_perquad_perpol);

% Noise = sigma*randn(Nsamples,1)/sqrt(2)+1i*sigma*randn(Nsamples,1)/sqrt(2);


% fprintf('The power noise generated is %.4f dBm\n', 10*log10(noise_power_perquad_perpol));

delay_phase_noise_distorted_RX_Apol = delay_phase_distorted_RX_Apol + Noise';

end