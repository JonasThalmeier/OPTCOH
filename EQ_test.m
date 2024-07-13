close all;clear;clc;    
r=1;
Rs=64;
OSNR_dB=100;
delta_nu=0e4;
rad_sec=0e3;
f_offset = 0e6;
EQ_N_tap=31;
EQ_mu=1e-5;
EQ_mu2=1e-5;
EQ_N1=3e5;
EQ_N2=1e5;
CarSync_DampFac=50;
EQ_mode = 'LMS';
% OSNR_dB = 5:5:20;
MODULATIONS = ["QPSK","16QAM","64QAM"];
modulation = ["QPSK","QAM","QAM"];
mod_mat = [1,1,2,2,3,3;64,128,64,128,64,128];


sweep_vec = logspace(4,5.5,10);
Ber_Tot = zeros(length(mod_mat), length(sweep_vec));

for mod = 1:length(mod_mat)
    r = mod_mat(1,mod);
    Rs = mod_mat(2,mod);
    Baud_rate = num2str(Rs);
    % load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));
    % Construct the file name dynamically
    fileName = sprintf('TXsequence_%s_%sGBaud.mat', MODULATIONS{r}, Baud_rate);
    % Construct the full path to the .mat file
    matFilePath = fullfile(fileparts(mfilename('fullpath')), 'TXsequences', fileName);
    % Load the .mat file
    load(matFilePath);

    if r == 1
        M = 4;
        power_norm = 2;
    elseif r == 2
        M = 16;
        power_norm = 10;
    else
        M = 64;
        power_norm = 42;
    end
    if size(SIG.Xpol.bits,2) ~=log2(M)
        fprintf('Right file not loaded\n');
        pause;
    end
    TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
    TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
    TX_SYMB= [repmat(SIG.Xpol.txSymb,10,1), repmat(SIG.Ypol.txSymb,10,1)]; % repeat the bits 10 times to simulate the original transmission in symb for LMS

    for idx=1:length(sweep_vec)
        % Create delay and phase convolved signals
        [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, sweep_vec(idx), rad_sec, SIG.symbolRate,f_offset);
        %add chromatic dispersion
        [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1,SIG.symbolRate);

        %% SIUMULATION
        Ber_Tot(mod,idx) = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac, 0);
    end
    figure(1);
    hold on;
    loglog(sweep_vec,Ber_Tot(mod,:),'o-','DisplayName', sprintf('r=%d, Rs=%d', mod_mat(1,mod),mod_mat(2,mod)));
    hold off;
end     
legend;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('delta nu [Hz]');
ylabel('BER');
grid on;
savefig('deltanu_sweep')








%--------Pol roation------------------------------------------------------
sweep_vec = logspace(4,6,10);
Ber_Tot = zeros(length(mod_mat), length(sweep_vec));

for mod = 1:length(mod_mat)
    r = mod_mat(1,mod);
    Rs = mod_mat(2,mod);
    Baud_rate = num2str(Rs);
    % load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));
    % Construct the file name dynamically
    fileName = sprintf('TXsequence_%s_%sGBaud.mat', MODULATIONS{r}, Baud_rate);
    % Construct the full path to the .mat file
    matFilePath = fullfile(fileparts(mfilename('fullpath')), 'TXsequences', fileName);
    % Load the .mat file
    load(matFilePath);

    if r == 1
        M = 4;
        power_norm = 2;
    elseif r == 2
        M = 16;
        power_norm = 10;
    else
        M = 64;
        power_norm = 42;
    end

    TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
    TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
    TX_SYMB= [repmat(SIG.Xpol.txSymb,10,1), repmat(SIG.Ypol.txSymb,10,1)]; % repeat the bits 10 times to simulate the original transmission in symb for LMS

    for idx=1:length(sweep_vec)
        % Create delay and phase convolved signals
        [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, sweep_vec(idx), SIG.symbolRate,f_offset);
        %add chromatic dispersion
        [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1,SIG.symbolRate);

        %% SIUMULATION
        Ber_Tot(mod,idx) = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac, 0);
    end
    figure(2);
    hold on;
    loglog(sweep_vec,Ber_Tot(mod,:),'o-','DisplayName', sprintf('r=%d, Rs=%d', mod_mat(1,mod),mod_mat(2,mod)));
    hold off;
end     
legend;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('polarisation rotation [rad/s]');
ylabel('BER');
grid on;
savefig('polrot_sweep')



%--------Freq offset------------------------------------------------------
sweep_vec = logspace(7,15,10);
Ber_Tot = zeros(length(mod_mat), length(sweep_vec));

for mod = 1:length(mod_mat)
    r = mod_mat(1,mod);
    Rs = mod_mat(2,mod);
    Baud_rate = num2str(Rs);
    % load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));
    % Construct the file name dynamically
    fileName = sprintf('TXsequence_%s_%sGBaud.mat', MODULATIONS{r}, Baud_rate);
    % Construct the full path to the .mat file
    matFilePath = fullfile(fileparts(mfilename('fullpath')), 'TXsequences', fileName);
    % Load the .mat file
    load(matFilePath);

    if r == 1
        M = 4;
        power_norm = 2;
    elseif r == 2
        M = 16;
        power_norm = 10;
    else
        M = 64;
        power_norm = 42;
    end

    TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
    TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
    TX_SYMB= [repmat(SIG.Xpol.txSymb,10,1), repmat(SIG.Ypol.txSymb,10,1)]; % repeat the bits 10 times to simulate the original transmission in symb for LMS

    for idx=1:length(sweep_vec)
        % Create delay and phase convolved signals
        [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate,sweep_vec(idx));
        %add chromatic dispersion
        [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1,SIG.symbolRate);

        %% SIUMULATION
        Ber_Tot(mod,idx) = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB, EQ_mode, EQ_N_tap, EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac, 0);
    end
    figure(3);
    hold on;
    loglog(sweep_vec,Ber_Tot(mod,:),'o-','DisplayName', sprintf('r=%d, Rs=%d', mod_mat(1,mod),mod_mat(2,mod)));
    hold off;
end     
legend;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('freq offset [Hz]');
ylabel('BER');
grid on;
savefig('foffset_sweep')





%--------mu------------------------------------------------------
sweep_vec = logspace(-6,-4,10);
Ber_Tot = zeros(length(mod_mat), length(sweep_vec));

for mod = 1:length(mod_mat)
    r = mod_mat(1,mod);
    Rs = mod_mat(2,mod);
    Baud_rate = num2str(Rs);
    % load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));
    % Construct the file name dynamically
    fileName = sprintf('TXsequence_%s_%sGBaud.mat', MODULATIONS{r}, Baud_rate);
    % Construct the full path to the .mat file
    matFilePath = fullfile(fileparts(mfilename('fullpath')), 'TXsequences', fileName);
    % Load the .mat file
    load(matFilePath);

    if r == 1
        M = 4;
        power_norm = 2;
    elseif r == 2
        M = 16;
        power_norm = 10;
    else
        M = 64;
        power_norm = 42;
    end

    TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
    TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
    TX_SYMB= [repmat(SIG.Xpol.txSymb,10,1), repmat(SIG.Ypol.txSymb,10,1)]; % repeat the bits 10 times to simulate the original transmission in symb for LMS

    for idx=1:length(sweep_vec)
        % Create delay and phase convolved signals
        [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate,f_offset);
        %add chromatic dispersion
        [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1,SIG.symbolRate);

        %% SIUMULATION
        Ber_Tot(mod,idx) = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB, EQ_mode, EQ_N_tap, sweep_vec(idx), EQ_mu2, EQ_N1, CarSync_DampFac, 0);
    end
    figure(4);
    hold on;
    loglog(sweep_vec,Ber_Tot(mod,:),'o-','DisplayName', sprintf('r=%d, Rs=%d', mod_mat(1,mod),mod_mat(2,mod)));
    hold off;
end     
legend;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('mu');
ylabel('BER');
grid on;
savefig('mu_sweep')



%--------NTap------------------------------------------------------
sweep_vec = linspace(11,101,10);
Ber_Tot = zeros(length(mod_mat), length(sweep_vec));

for mod = 1:length(mod_mat)
    r = mod_mat(1,mod);
    Rs = mod_mat(2,mod);
    Baud_rate = num2str(Rs);
    % load(strcat('TXsequences/TXsequence_', MODULATIONS(r) , '_',Baud_rate,'GBaud.mat'));
    % Construct the file name dynamically
    fileName = sprintf('TXsequence_%s_%sGBaud.mat', MODULATIONS{r}, Baud_rate);
    % Construct the full path to the .mat file
    matFilePath = fullfile(fileparts(mfilename('fullpath')), 'TXsequences', fileName);
    % Load the .mat file
    load(matFilePath);

    if r == 1
        M = 4;
        power_norm = 2;
    elseif r == 2
        M = 16;
        power_norm = 10;
    else
        M = 64;
        power_norm = 42;
    end

    TX_BITS_Xpol = repmat(SIG.Xpol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
    TX_BITS_Ypol = repmat(SIG.Ypol.bits,10,1); %repeat the bits 10 times to simulate the original transmission
    TX_SYMB= [repmat(SIG.Xpol.txSymb,10,1), repmat(SIG.Ypol.txSymb,10,1)]; % repeat the bits 10 times to simulate the original transmission in symb for LMS

    for idx=1:length(sweep_vec)
        % Create delay and phase convolved signals
        [X_distorted, Y_distorted] = DP_Distortion(SIG.Xpol.txSig, SIG.Ypol.txSig, delta_nu, rad_sec, SIG.symbolRate,f_offset);
        %add chromatic dispersion
        [X_CD,Y_CD]=Chromatic_Dispersion(X_distorted, Y_distorted, SIG.Sps, 1,SIG.symbolRate);

        %% SIUMULATION
        Ber_Tot(mod,idx) = core_simulation(X_CD, Y_CD, r, Rs, OSNR_dB, EQ_mode, sweep_vec(idx), EQ_mu, EQ_mu2, EQ_N1, CarSync_DampFac, 0);
    end
    figure(5);
    hold on;
    plot(sweep_vec,Ber_Tot(mod,:),'o-','DisplayName', sprintf('r=%d, Rs=%d', mod_mat(1,mod),mod_mat(2,mod)));
    hold off;
end     
legend;
xlabel('NTaps');
ylabel('BER');
grid on;
savefig('Ntaps')