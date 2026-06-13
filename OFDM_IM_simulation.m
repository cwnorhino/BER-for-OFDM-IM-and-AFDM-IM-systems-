%% OFDM-IM BER vs Eb/N0 Simulation
clear; clc; close all;

%% System parameters
N = 4;           % subcarriers per subblock
K = 2;           % active subcarriers per subblock
M = 4;           % QPSK modulation
cpLen = 16;      % cyclic prefix length
N_F = 128;       % total subcarriers (must be multiple of N)
G = N_F / N;     % number of subblocks

% Derived parameters
p1 = floor(log2(nchoosek(N, K)));    % index bits per subblock
p2 = K * log2(M);                    % modulation bits per subblock
p = p1 + p2;                         % total bits per subblock
bits_per_frame = G * p;              % total bits per OFDM-IM symbol
Eb = K / p;                          % Energy per bit 

%% Channel parameters 
L = 10;                  

%% Eb/N0 range (dB) 
EbN0_dB = 0:2:24;
numEbN0 = length(EbN0_dB);

%% Simulation settings
numFrames = 500;         
maxBitErrors = 1000;     
minFrames = 50;          

% Preallocate BER array
BER = zeros(1, numEbN0);

%% simulation loop
for snrIdx = 1:numEbN0
    EbN0 = EbN0_dB(snrIdx);
    totalBits = 0;
    bitErrors = 0;
    frameCount = 0;

    
    EbN0_linear = 10^(EbN0 / 10);
    N0 = Eb / EbN0_linear; 

    while (frameCount < numFrames) && (bitErrors < maxBitErrors || frameCount < minFrames)
        % Generate random bits
        txBits = randi([0,1], bits_per_frame, 1);

        % Transmitter
        txSig = transmitter(txBits, N, K, M, cpLen, N_F);

        
        % Generate new random channel realization for each frame (Rayleigh fading)
        h = sqrt(1/L) * (randn(L,1) + 1i*randn(L,1)) / sqrt(2);

        % Convolve tx signal with channel
        rxSig_filtered = conv(txSig, h);

        % Truncate to original length 
        rxSig_filtered = rxSig_filtered(1:length(txSig));

        % Add AWGN
        noise = sqrt(N0/2) * (randn(size(rxSig_filtered)) + 1i*randn(size(rxSig_filtered)));
        rxSig = rxSig_filtered + noise;

        % Compute the frequency-domain channel coefficients needed by the receiver
        H_freq = fft(h, N_F);   

        % Receiver
        rxBits = receiver(rxSig, N, K, M, cpLen, N_F, H_freq, N0);

        % Count errors
        errors = sum(rxBits ~= txBits);
        bitErrors = bitErrors + errors;
        totalBits = totalBits + length(txBits);
        frameCount = frameCount + 1;
    end

    BER(snrIdx) = bitErrors / totalBits;
    fprintf('Eb/N0 = %d dB, BER = %.2e (%d errors / %d bits)\n', ...
        EbN0, BER(snrIdx), bitErrors, totalBits);
end

%% Plot results
figure;
semilogy(EbN0_dB, BER, 'b-o', 'LineWidth', 1.5);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title(sprintf('OFDM-IM Performance (N=%d, K=%d, M=%d, N_F=%d)', N, K, M, N_F));
legend('OFDM-IM LLR Detector');
hold on;