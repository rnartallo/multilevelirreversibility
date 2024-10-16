function [Dist] = DirichletPrior4D(Dist)
%Dist - distribution containing zeros
[L,~]= size(Dist);
Dist = Dist*L^4;
m = sum(sum(sum(sum(Dist))));
Dist = (Dist+1)/(m + L^4);
end