function [X_out,Y_out] = vit_n_vit(X_in,Y_in, rate, init)
% Viterbi & Viterbi Algorithm for Phase Estimation
M = 4; % 4th power for 16QAM, treating it similar to QPSK in phase recovery
% Raise received symbols to the Mth power
X_raised = (X_in(1:2:end)./abs(X_in(1:2:end))).^M;
Y_raised = (Y_in(1:2:end)./abs(Y_in(1:2:end))).^M;

N_X = length(X_raised);
N_Y = length(Y_raised);
% Calculate the phase using the moving average of the raised signal
% phase_X = zeros(length(X_in),1);
% phase_Y = zeros(length(Y_in),1);
% initphase_X = mean(X_raised(1:init));
% initphase_Y = mean(Y_raised(1:init));
% X_raised = X_raised*exp(-1i*initphase_X);
% Y_raised = Y_raised*exp(-1i*initphase_Y);
% for k = 2:N_X
%     phase_X(2*k) = ((1-rate)*phase_X(2*(k-1)))+rate*angle(X_raised(k-1));
%     X_raised(k) = X_raised(k)*exp(-1i*phase_X(2*k));
% end
% for k = 2:N_Y
%     phase_Y(2*k) = ((1-rate)*phase_Y(2*(k-1)))+rate*angle(Y_raised(k-1));
%     Y_raised(k) = Y_raised(k)*exp(-1i*phase_Y(2*k));
% end


windowSize = 21; % Window size for moving average

% Calculate the phase using the moving average of the raised signal
avgPhase_X = zeros(N_X, 1);
for k = (windowSize+1)/2:N_X-(windowSize-1)/2
    windowedSegment = X_raised(k-(windowSize-1)/2:k+(windowSize-1)/2);
    avgPhase_X(k) = angle(mean(windowedSegment));
end
avgPhase_Y = zeros(N_Y, 1);
for k = (windowSize+1)/2:N_Y-(windowSize-1)/2
    windowedSegment = Y_raised(k-(windowSize-1)/2:k+(windowSize-1)/2);
    avgPhase_Y(k) = angle(mean(windowedSegment));
end

% Normalize the phase estimate and compensate
X_phase1 = repelem(avgPhase_X, 2);
X_phase2 = X_phase1(1:length(X_in));
Y_phase1 = repelem(avgPhase_Y, 2);
Y_phase2 = Y_phase1(1:length(Y_in));
X_out = X_in.*exp(-1i.*X_phase2./M);
Y_out = Y_in.*exp(-1i.*Y_phase2./M);
end

