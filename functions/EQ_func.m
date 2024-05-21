% function [out] = EQ(Xpol_in,mu)
% num_taps = 9;  % Number of filter taps
% fir = zeros(1, num_taps);  % Initialize the FIR filter coefficients
% N = length(Xpol_in);  % Length of the input signal
% XPol_out = zeros(N - num_taps + 1, 1);  % Output signal
% fir(ceil(num_taps / 2)) = 1;  % Initialize the center tap to 1
% XPol_out(1) = fir*Xpol_in(1:num_taps);
% for idx=1:N-num_taps
%     % Calculate the error term for the CMA update
%     er = conj(XPol_out(idx))*(1-abs(XPol_out(idx))^2);
%     % Update the FIR filter coefficients using CMA
%     fir = fir +  mu*er.*Xpol_in(idx:idx+num_taps-1)';
%     XPol_out(idx+1) = fir*Xpol_in(idx+1:idx+num_taps);
% end

function [out] = EQ_func(Xpol_in,Ypol_in,mu,NTaps,alg,Xorg,Yorg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADAPTIVEEQUALIZER [y] = AdaptiveEqualizer(x,SpS,ParamDE) %
% %
% This function performs equalization in a dual polarization signal 'y' %
% using the CMA or the RDE algorithms. The supported modulation formats %
% are QPSK and 16-QAM. The CMA can be used for pre-convergence of the RDE %
% algorithm. %
% %
% Input: %
% x = Input signal for two pol. orientation (matrix with two %
% column vectors, where each column corresponds to the signal %
% of one pol. orientation). 'x' must be normalized to unitary %
% power and obtained at 2 Sa/Symbol; %
% SpS = Number of samples per symbol in the input signal; %
% ParamDE =Struct that specifies parameters for the adaptive equalization%
% - ParamDE.Eq: Defines the algorithm to be used: %
% 'CMA' - CMA only; %
% 'RDE' - RDE only; %
% 'CMA+RDE' - CMA is used to initialize the RDE; %
% - ParamDE.NTaps: Number of taps for the filters in the butterfly %
% configuration; %
% - ParamDE.Mu: Step-size for coefficients calculation; %
% - ParamDE.SingleSpike: 'true': Single spike initialization; %
% 'false': All taps are initialized with zeros%
% - ParamDE.N1: Number of coefficient calculations to perform prior %
% to proper initialization of the filters w2H and w2V; %
% ('ParamDE.N1' is related to the single spike initialization) %
% - ParamDE.N2: Number of coefficient calculations to perform prior %
% to swicth from CMA to RDE (Note that N2 must only be %
% defined if CMA is used for RDE initialization); %
% - ParamDE.NOut: Number of samples to discard after equalization; %
% %
% Output: %
% y = Output signal (at 1 Sa/Symbol) after adaptive equalization; %
% %
% This function is part of the book Digital Coherent Optical Systems; %
% Darli A. A. Mello and Fabio A. Barbosa; %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_CMA = 1;
SpS = 1;
NOut = 0;
train_len = length(Xorg);
% Input blocks:
x = [Xpol_in, Ypol_in];
x = [x(end-floor(NTaps/2)+1:end,:) ; x ; x(1:floor(NTaps/2),:)];
xV = convmtx(x(:,1).',NTaps) ; xH = convmtx(x(:,2).',NTaps);
% xV = xV(:,NTaps:SpS:end-NTaps+1) ; xH = xH(:,NTaps:SpS:end-NTaps+1); %
% Not needed if used wit 1SpS
% Output length:
OutLength = floor((size(x,1)-NTaps+1)) ; clearvars x % removed the /2
% Initializing the outputs
y1 = zeros(OutLength,1) ; y2 = zeros(OutLength,1);
% Initial filter coefficients:
w1V = zeros(NTaps,1); w1H = zeros(NTaps,1); w2V = zeros(NTaps,1);w2H = zeros(NTaps,1);
w1V(floor(NTaps/2)+1) = 1;w1H(floor(NTaps/2)+1) = 1;w2V(floor(NTaps/2)+1) = 1;w2H(floor(NTaps/2)+1) = 1;
for i = 1:OutLength
    % Calculating the outputs:
    % y1(i) = w1V'*xV(:,i) + w1H'*xH(:,i);
    % y2(i) = w2V'*xV(:,i) + w2H'*xH(:,i);
    y1(i) = w1V'*xV(:,i);
    y2(i) = w2H'*xH(:,i);
    % Updating the filter coefficients:
    if alg=="CMA"
        [w1V,w1H,w2V,w2H] = CMA(xV(:,i),xH(:,i),y1(i),y2(i),w1V,w1H,w2V,w2H,R_CMA,mu);
    elseif alg=="LMS"
        [w1V,w1H,w2V,w2H] = LMS(xV(:,i),xH(:,i),y1(i),y2(i),w1V,w1H,w2V,w2H,mu,Xorg(i),Yorg(i));
        if i>=train_len
            alg = "";
        end
    end
end
% Output samples:
out = [y1 y2] ; out = out(1+NOut:end,:);
end

function [w1V,w1H,w2V,w2H] = CMA(xV,xH,y1,y2,w1V,w1H,w2V,w2H,R,Mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CMA [w1V,w1H,w2V,w2H] = CMA(xV,xH,y1,y2,w1V,w1H,w2V,w2H,R,Mu) %
% %
% This function performs the update of the filters of the MIMO butterfly %
% equalizer using the CMA algorithm. %
% %
% Input: %
% xV = Column vector that represents the complex samples at the MIMO %
% butterfly equalizer input for the vertical pol. orientation; %
% xH = Column vector that represents the complex samples at the MIMO %
% butterfly equalizer input for the horizontal pol. orientation; %
% y1 = Sample at the output 1 of the MIMO butterfly equalizer; %
% y2 = Sample at the output 2 of the MIMO butterfly equalizer; %
% w1V, w1H, w2V, w2H = N-coefficient FIR filters that compose the MIMO %
% butterfly equalizer; %
% R = radius used as reference for coefficients calculation; %
% Mu = Step-size for coefficients calculation; %
% Note: xV and xH must have the same length as the FIR filters w1V, w1H,%
% w2V, w2H; %
% %
% Output: %
% w1V, w1H, w2V, w2H = Updated N-coefficient FIR filters that compose %
% the MIMO butterfly equalizer; %
% %
% This function is part of the book Digital Coherent Optical Systems; %
% Darli A. A. Mello and Fabio A. Barbosa; %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Updating the filters:
w1V = w1V + Mu*xV*(R-abs(y1).^2)*conj(y1);
w1H = w1H + Mu*xH*(R-abs(y1).^2)*conj(y1);
w2V = w2V + Mu*xV*(R-abs(y2).^2)*conj(y2);
w2H = w2H + Mu*xH*(R-abs(y2).^2)*conj(y2);
end

function [w1V,w1H,w2V,w2H] = LMS(xV,xH,y1,y2,w1V,w1H,w2V,w2H,Mu,Xorg,Yorg)
%Updating the filters:
w1V = w1V + Mu*xV*y1*conj(Xorg-y1);
w1H = w1H + Mu*xH*y1*conj(Yorg-y1);
w2V = w2V + Mu*xV*y2*conj(Xorg-y2);
w2H = w2H + Mu*xH*y2*conj(Yorg-y2);
end