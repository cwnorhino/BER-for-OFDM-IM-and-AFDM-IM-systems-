function txSignal = transmitter(bits, N, K, M, cpLen, N_F)
% OFDM-IM Transmitter Parameters and definitions
%   bits    - input bit vector (length = G * (p1 + p2))
%   N       - subcarriers per subblock
%   K       - active subcarriers per subblock (1 <= K < N)
%   G       - subblock 
%   M       - QAM order (power of 2)
%   cpLen   - cyclic prefix length
%   N_F     - total subcarriers (multiple of N)
%   txSignal - time-domain transmitted signal (column vector)

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


    % Convert p1 bits to integer using bit2int (MSB first)
    idx_val = bit2int(idx_bits(:), p1, true);  
    active_local = combination_from_index(idx_val, N, K);

    % M-ary symbols
    sym_bits_matrix = reshape(sym_bits, log2(M), []);
    sym_ints = bit2int(sym_bits_matrix, log2(M), true);
    symbols = qammod(sym_ints, M, 'UnitAveragePower', true);

 
    global_pos = (g-1)*N + active_local + 1;
    X_freq(global_pos) = symbols;
end

% Block interleaving
X_matrix = reshape(X_freq, N, G).'; 
interleaved = X_matrix(:);

% IFFT and Cyclic prefix
x_time = ifft(interleaved) * sqrt(N_F);
txSignal = [x_time(end-cpLen+1:end); x_time];
end

% Combinatorial Number theory method 
function active_local = combination_from_index(idx_val, N, K)
active_local = zeros(1, K);
current = idx_val;

for k = K:-1:1
    n = N - 1;


    while true  
     
        if n < k
            val = 0;
        else
            val = nchoosek(n, k);
        end

        if val <= current
            break;
        end
        n = n - 1;
    end

    active_local(k) = n;

    
    if n >= k
        current = current - nchoosek(n, k);
    end
end
end