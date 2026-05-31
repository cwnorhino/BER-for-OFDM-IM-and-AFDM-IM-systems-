function txSignal = transmitter(bits,N,K,M,cpLen)

%splitting of bits into index and modulation bits

indexBits = bits(:,1:2);
symbolBits = bits(:,3:6);


%Symbol activation pattern lookup table 
patterns = [
    1 2;
    1 3;
    2 4;
    3 4
    ];

%converting bits into decimal values 
%correction:Changed bi2de to bit2int as bi2de is outdated acc to docs

symbolInts = bit2int(reshape(symbolBits.',2,[]).',2);

%normalization 
symbols = qammod(symbolInts,M,...
    'UnitAveragePower',true);

symbols = reshape(symbols,[],K);

%Frequency division vector
X = zeros(size(bits,1),N);

for b = 1:size(bits,1)

%convert index bits into Pattern number 

    idx = bit2int(indexBits(b,:),2) + 1;

%active subcarrier 

    active = patterns(idx,:);
%Place symbols on active carriers 
    X(b,active(1)) = symbols(b,1);
    X(b,active(2)) = symbols(b,2);

end

%Inverse Fourier transform 

x_time = ifft(X,N,2);

%Adding Cyclic Prefix

x_cp = [x_time(:,end-cpLen+1:end) x_time];
%serialize 

txSignal = x_cp.';
txSignal = txSignal(:);

end