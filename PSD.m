clc; clear; close all;

%% parameters
N = 64;          % Number of subcarriers
M = 4;           % Modulation order (QPSK)
numSymbols = 100; % Number of symbols
K = 4;          % Overlapping factor for FBMC

%% generate QPSK symbols
data = randi([0 M-1], N, numSymbols);
qpskMod = pskmod(data, M, pi/4, 'gray');

%% OFDM Modulation
ofdmSymbols = ifft(qpskMod, N, 1); % IFFT across columns
ofdmSignal = reshape(ofdmSymbols, [], 1); % Serializing

%% compute PSD for OFDM
[pxx_ofdm, f_ofdm] = pwelch(ofdmSignal, hamming(512), 256, 1024, 'centered');

%% FBMC Modulation (Using PHYDYAS Filter)
h = phydyas_filter(K, N); % Compute prototype filter
filterLen = length(h);

% allocate FBMC signal with sufficient length
fbmcSignal = zeros((numSymbols + K - 1) * N, 1);

for symbolIdx = 1:numSymbols
    for k = 1:N
        % generate FBMC symbol with offset
        tempSignal = conv(qpskMod(k, symbolIdx), h); % Apply prototype filter
        
        % Compute correct start index (each symbol is spaced N samples apart)
        startIdx = (symbolIdx - 1) * N + 1;
        endIdx = startIdx + length(tempSignal) - 1;
        
        % Ensure indexing is within bounds
        if endIdx > length(fbmcSignal)
            endIdx = length(fbmcSignal);
            tempSignal = tempSignal(1:endIdx - startIdx + 1);
        end
        
        % add to FBMC signal with proper summation
        fbmcSignal(startIdx:endIdx) = fbmcSignal(startIdx:endIdx) + tempSignal;
    end
end

%% Compute PSD for FBMC
[pxx_fbmc, f_fbmc] = pwelch(fbmcSignal, hamming(512), 256, 1024, 'centered');

%% plot PSD
figure;
plot(f_ofdm, 10*log10(pxx_ofdm), 'r', 'LineWidth', 1.5); hold on;
plot(f_fbmc, 10*log10(pxx_fbmc), 'b', 'LineWidth', 1.5);
grid on;
legend('OFDM PSD', 'FBMC PSD');
xlabel('Frequency (normalized)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density of OFDM and FBMC');

%% function: PHYDYAS Prototype Filter
function h = phydyas_filter(K, N)
    % PHYDYAS prototype filter coefficients for FBMC
    a = [1 0.97196 1/sqrt(2) 0.23515]; % Coefficients for K = 4
    h = zeros(K*N, 1); % Initialize filter

    % Construct the filter
    for n = 0:(K*N-1)
        h(n+1) = a(1) * cos(pi * (n - (K*N)/2) / (K*N)) ...
               + a(2) * cos(3*pi * (n - (K*N)/2) / (K*N)) ...
               + a(3) * cos(5*pi * (n - (K*N)/2) / (K*N)) ...
               + a(4) * cos(7*pi * (n - (K*N)/2) / (K*N));
    end

    % Normalize the filter energy
    h = h / sum(h);
end
