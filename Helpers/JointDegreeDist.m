function P = JointDegreeDist(K,alpha,beta,a,b)
%Returns the probability of having a node with degree a in layer alpha and
%degree b in layer beta
%This has to be run case by case as computing the whole matrix takes too
%much memory
[T,~]=size(K);
Ba = K==a;
Bb = K==b;
J = 0;
for node=1:T
    if Ba(node,alpha) && Bb(node,beta)
        J=J+1;
    end
end
P=J/T;
end