function [downsampledSig_Xpol,downsampledSig_Ypol] = samp_time_recovery(Xpol,Ypol,Sps)
% SAMP_PHASE_RECOVERY This function performs sample phase recovery on 
% dual-polarization signals by determining the optimal sampling phase.
% The function first interpolates the input signals for both polarizations
% (Xpol and Ypol) by a factor of Sps/2. It then calculates the timing error
% for each sample in the interpolated signals using a derivative-based
% method. The error is computed using the difference between late and early
% samples, multiplied by the complex conjugate of the mid sample.
%
% Inputs:
% Xpol - Array containing the signal of the X polarization.
% Ypol - Array containing the signal of the Y polarization.
% Sps - Samples per symbol in the original signal before interpolation.
%
% Outputs:
% downsampledSig_Xpol - Downsampled signal of X polarization, sampled at the
%                       optimal phase.
% downsampledSig_Ypol - Downsampled signal of Y polarization, sampled at the
%                       optimal phase.
%
% The function uses interpolation to double the sampling points, aiding in the
% fine-tuning of the sampling phase. By evaluating the error terms across a range
% of potential sampling points (from the start to the original Sps value), the 
% function identifies the minimum cumulative error phase for resampling. This
% optimal sampling phase is used to downsample the interpolated signals to half
% the interpolated rate, effectively aligning the sampling points with the 
% estimated optimal phase for symbol decision-making.
Xpol = interp(Xpol,Sps/2);
Ypol = interp(Ypol,Sps/2);
earlyX = Xpol(1:end-2);
lateX = Xpol(3:end);
midX = Xpol(2:end-1);
earlyY = Ypol(1:end-2);
lateY = Ypol(3:end);
midY = Ypol(2:end-1);
erX = real(conj(midX).*(lateX-earlyX));
erY = real(conj(midY).*(lateY-earlyY));
N = length(erX)/Sps;
cumEr = zeros(Sps,1);
for k=1:Sps+1
    cumEr(k) = (sum(erX(k:Sps:end))+sum(erY(k:Sps:end)))/N;
end
[M IND] = min(cumEr);
% figure;
% for k=0:Sps
%     subplot(Sps+1, 1, k+1);
%     scatter(real(Xpol(IND+k:Sps:end)), imag(Xpol(IND+k:Sps:end)));
%     title(k);
% end
IND
downsampledSig_Xpol = Xpol(IND-1:Sps/2:end);
downsampledSig_Ypol = Ypol(IND-1:Sps/2:end);
end

