function [Dist] = DirichletPrior5D(Dist)
%Dist - distribution containing zeros
[L,~]= size(Dist);
Dist = Dist*L^5;
m = sum(sum(sum(sum(sum(Dist)))));
Dist = (Dist+1)/(m + L^5);
end