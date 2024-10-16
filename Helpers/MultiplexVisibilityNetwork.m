function [AdjMatrix] = MultiplexVisibilityNetwork(X)
tic
[N,T] = size(X);
%This function generates a horizontal multilayer visibility network for a time-series
AdjMatrix = zeros(T,T,N);
for n=1:N
    for t1=1:T-1
        for t2=t1+1:T
            connected = true;
            for k=t1+1:t2-1
                if X(n,k)>=X(n,t2)+(X(n,t1)-X(n,t2))*((t2-k)/(t2-t1))
                    connected=false;
                end
            end
            if connected
                AdjMatrix(t1,t2,n)=1;
                AdjMatrix(t2,t1,n)=1;
            end
        end
    end
end
toc
end