function [Dist] = DirichletPrior1D(Dist)
%Dist - distribution containing zeros
[L,N]= size(Dist);
Dist = Dist*L;
m = sum(Dist);
for n=1:N
    for d=1:L
        Dist(d,n) = (Dist(d,n)+1)/(m(n) + L);
    end
end