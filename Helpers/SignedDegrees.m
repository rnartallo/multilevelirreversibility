function [K_in,K_out] = SignedDegrees(A)
%A - multilayer adjacency matrix
%K_in - in-degree matrix
%K_out - out-degree matrix
[~,T,N] = size(A);
K_in = zeros([T,N]);
K_out = zeros([T,N]);
for m=1:N
    for t1=1:T
        for t2=1:T
            if t1<t2
                K_in(t1,m)=K_in(t1,m)+A(t1,t2,m);
            elseif t1>t2
                K_out(t1,m)=K_out(t1,m)+A(t1,t2,m);
            end
        end
    end
end
end