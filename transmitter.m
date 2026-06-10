function txSignal = transmitter(bits, N, K, M, cpLen, N_F)
% OFDM-IM Transmitter Parameters and defintions
%   bits    - input bit vector (length = G * (p1 + p2))
%   N       - subcarriers per subblock
%   K       - active subcarriers per subblock (1 <= K < N)
%   G       - subblock 
%   M       - QAM order (power of 2)
%   cpLen   - cyclic prefix length
%   N_F     - total subcarriers (multiple of N)
%   txSignal - time-domain transmitted signal (column vector

assert(mod(N_F, N) == 0, 'N_F must be multiple of N');
G = N_F / N;  
p1 = floor(log2(nchoosek(N, K)));
p2 = K * log2(M);
p = p1 + p2;
total_bits = G * p;
assert(length(bits) == total_bits, 'Bit length mismatch');

% Reshape bits into subblocks (G rows, p columns)
bits_matrix = reshape(bits, p, G).';

X_freq = zeros(N_F, 1);

for g = 1:G
    idx_bits = bits_matrix(g, 1:p1);
    sym_bits = bits_matrix(g, p1+1:end);

    % Index bits -> active subcarrier bits
    % Convert p1 bits to integer using bit2int (MSB first)
    idx_val = bit2int(idx_bits, p1, 'msb');   % integer 0..2^p1-1
    active_local = combination_from_index(idx_val, N, K);

    % M-ary symbols
    % Reshape sym_bits into K groups of log2(M) bits
    sym_bits_matrix = reshape(sym_bits, log2(M), []).';
    sym_ints = bit2int(sym_bits_matrix, log2(M), 'msb');
    symbols = qammod(sym_ints, M, 'UnitAveragePower', true);

    % Place symbols into global frequency vector
    global_pos = (g-1)*N + active_local + 1;
    X_freq(global_pos) = symbols;
end

% Block interleaving
X_matrix = reshape(X_freq, N, G).'; %the transpose might be an issue idk yet
interleaved = X_matrix(:);

% IFFT and Cyclic prefix
x_time = ifft(interleaved) * sqrt(N_F);
txSignal = [x_time(end-cpLen+1:end); x_time];
end

% Combinatorial Number theory method 
function indices = combination_from_index(idx, N, K)
indices = zeros(1, K);
current = idx;
n = N - 1;
for k = K:-1:1
    if k == 1
        indices(K) = current;
        break;
    end
    while n >= k && nchoosek(n, k) > current
        n = n - 1;
    end
    indices(K - k + 1) = n;
    current = current - nchoosek(n, k);
    n = n - 1;
end
indices = sort(indices);
end
