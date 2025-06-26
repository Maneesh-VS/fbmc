clc; clear; close all;

% Parameters
M = 4; % QPSK
k = log2(M);
numSymbols = 1e4;
numSubcarriers = 64;
K = 4; % FBMC overlapping factor

% QPSK Modulator
mod = comm.QPSKModulator('BitInput', true);
bits = randi([0 1], numSymbols * k, 1);
symbols = mod(bits);

% --- Trim symbols to match OFDM block size ---
numSymbolsTrimmed = floor(length(symbols)/numSubcarriers) * numSubcarriers;
symbols = symbols(1:numSymbolsTrimmed);

%% ---------------- OFDM ---------------- %%
ofdmSymbols = reshape(symbols, numSubcarriers, []);
tx_ofdm = ifft(ofdmSymbols, numSubcarriers);
tx_ofdm = tx_ofdm(:); % Serializing

% --- PAPR (OFDM) ---
papr_ofdm = 10*log10(max(abs(tx_ofdm).^2) / mean(abs(tx_ofdm).^2));
fprintf('PAPR (OFDM): %.2f dB\n', papr_ofdm);

%% ---------------- FBMC (Approximate) ---------------- %%
L = K * numSubcarriers;                 % Prototype filter length
h = sqrt(hamming(L));                   % Prototype filter

symbols_fbmc = reshape(symbols, numSubcarriers, []);
numBlocks = size(symbols_fbmc, 2);

% Preallocate FBMC signal length
fbmc_block_len = length(conv(h, upsample(ones(numSubcarriers,1), K)));
tx_fbmc_len = (numBlocks-1)*numSubcarriers + fbmc_block_len;
tx_fbmc = zeros(tx_fbmc_len, 1);

for n = 1:numBlocks
    s = symbols_fbmc(:,n);
    upsampled = upsample(s, K);               % Upsample by K
    filtered = conv(h, upsampled);            % Apply prototype filter

    idx_start = (n-1)*numSubcarriers + 1;
    idx_end = idx_start + length(filtered) - 1;

    tx_fbmc(idx_start:idx_end) = tx_fbmc(idx_start:idx_end) + filtered;
end

% --- PAPR (FBMC) ---
papr_fbmc = 10*log10(max(abs(tx_fbmc).^2) / mean(abs(tx_fbmc).^2));
fprintf('PAPR (FBMC): %.2f dB\n', papr_fbmc);

%% ---------------- CCDF for PAPR ---------------- %%
N = 1000; % Trials for CCDF
paprValsOFDM = zeros(N,1);
paprValsFBMC = zeros(N,1);

for i = 1:N
    % Generate random symbols
    bits = randi([0 1], numSubcarriers * k, 1);
    sym = mod(bits);
    
    % OFDM block
    ofdm_block = ifft(sym);
    paprValsOFDM(i) = 10*log10(max(abs(ofdm_block).^2) / mean(abs(ofdm_block).^2));
    
    % FBMC block (approximate)
    up_fbmc = upsample(sym, K);
    filtered_fbmc = conv(h, up_fbmc);
    paprValsFBMC(i) = 10*log10(max(abs(filtered_fbmc).^2) / mean(abs(filtered_fbmc).^2));
end

% CCDF (1 - CDF)
[p_ofdm, x_ofdm] = ecdf(paprValsOFDM);
ccdf_ofdm = 1 - p_ofdm;

[p_fbmc, x_fbmc] = ecdf(paprValsFBMC);
ccdf_fbmc = 1 - p_fbmc;

%% ---------------- Plot CCDF ---------------- %%
figure;
semilogy(x_ofdm, ccdf_ofdm, 'r', 'LineWidth', 2); hold on;
semilogy(x_fbmc, ccdf_fbmc, 'b', 'LineWidth', 2);
grid on;
xlabel('PAPR (dB)');
ylabel('CCDF (P[PAPR > PAPR₀])');
legend('OFDM', 'FBMC', 'Location', 'southwest');
xlim([2 11]);
ylim([1e-3 1]);
set(gca, 'FontSize', 12);
title('PAPR CCDF of OFDM and FBMC', 'FontWeight', 'bold');
