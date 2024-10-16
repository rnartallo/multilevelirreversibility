function P = DegreeDist(K)
%This calculates the degree dist matrix from degree matrix K
[T,N]=size(K);
maxDegree = max(max(K));
P=zeros(maxDegree,N);
for k=1:maxDegree
    [~,c] = find(K==k);
    for n=1:N
        P(k,n) = sum(c==n)/T;
    end
end