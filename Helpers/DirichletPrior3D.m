function [Dist] = DirichletPrior3D(Dist)
%Dist - distribution containing zeros
[L,~]= size(Dist);
Dist = Dist*L^3;
m = sum(sum(sum(Dist)));
Dist = (Dist+1)/(m + L^3);
end