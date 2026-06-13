function rxBits = receiver(rxSignal, N, K, M, cpLen, N_F, H_freq, N0)
% OFDM-IM Receiver using LLR detector
G = N_F / N;
p1 = floor(log2(nchoosek(N, K)));
p2 = K * log2(M);
p = p1 + p2;

% 1. Remove Cyclic Prefix
rx_time = rxSignal(cpLen+1:end);

% 2. FFT
Y_interleaved = fft(rx_time) / sqrt(N_F);

% 3. Block deinterleaving
Y_matrix = reshape(Y_interleaved, G, N);          
Y_deinterleaved = reshape(Y_matrix.', N_F, 1);

H_matrix = reshape(H_freq, G, N);
H_deinterleaved = reshape(H_matrix.', N_F, 1);

% 4. Precompute constellation
const = qammod(0:M-1, M, 'UnitAveragePower', true);

rxBits = zeros(G * p, 1);

for g = 1:G
    y = Y_deinterleaved((g-1)*N+1 : g*N);
    h = H_deinterleaved((g-1)*N+1 : g*N);

    % Compute LLR for each subcarrier
    LLR = zeros(N,1);
    for i = 1:N
        diff = y(i) - h(i)*const;
        sum_exp = sum(exp(-(abs(diff).^2)/N0));
        LLR_h1 = log(sum_exp) - log(M);
        LLR_h0 = -(abs(y(i)).^2)/N0;
        LLR(i) = LLR_h1 - LLR_h0;
    end

    % Select K active indices
    [~, sortedIdx] = sort(LLR, 'descend');
    active_local = sort(sortedIdx(1:K)) - 1;  

    % Demodulate symbols on active subcarriers
    sym_ints = zeros(K,1);
    for idx = 1:K
        pos = active_local(idx) + 1;
        metric = abs(y(pos) - h(pos)*const).^2;
        [~, minIdx] = min(metric);
        sym_ints(idx) = minIdx - 1;
    end

   
    idx_val = index_from_combination(active_local, K);


    if idx_val >= 2^p1
        
        idx_val = 0; 
    end

    idx_bits = de2bi(idx_val, p1, 'left-msb');

    % Modulation bits
    sym_bits = de2bi(sym_ints, log2(M), 'left-msb');
    sym_bits = reshape(sym_bits.', 1, []);

    % Assemble subblock bits
    subblock_bits = [idx_bits, sym_bits];
    rxBits((g-1)*p+1 : g*p) = subblock_bits.';
end
end

% Combinatorial Inverse Function

function idx = index_from_combination(indices, K)

idx = 0;
for k = 1:K
    c = indices(k);
    if c >= k
        idx = idx + nchoosek(c, k);
    end
end
end