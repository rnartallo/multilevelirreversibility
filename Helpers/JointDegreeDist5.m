function P = JointDegreeDist5(K,alpha,beta,theta,delta,phi,a,b,c,d,f)
%Returns the probability of having a node with degree a in layer alpha and
%degree b in layer beta and degree c in layer theta
%This has to be run case by case as computing the whole matrix takes too
%much memory
[T,~]=size(K);
Ba = K==a;
Bb = K==b;
Bc = K==c;
Bd = K==d;
Bf = K==f;
J = 0;
for node=1:T
    if Ba(node,alpha) && Bb(node,beta) && Bc(node,theta) && Bd(node,delta) && Bf(node,phi)
        J=J+1;
    end
end
P=J/T;
end