function [Dist] = DirichletPrior2D(Dist)
%Dist - distribution containing zeros
[L,~]= size(Dist);
Dist = Dist*L^2;
m = sum(sum(Dist));
Dist = (Dist+1)/(m + L^2);
end