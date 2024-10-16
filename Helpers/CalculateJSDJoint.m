function [JSD] = CalculateJSDJoint(DistIn,DistOut)
[MaxDegree,~] = size(DistIn);
%Smooth using Dirichlet prior
DistInSmoothed = DirichletPrior2D(DistIn);
DistOutSmoothed = DirichletPrior2D(DistOut);
DistMixed = 0.5*(DistOutSmoothed+DistInSmoothed);
KL_InMixed=0;
KL_OutMixed=0;
for x =1:MaxDegree
    for y=1:MaxDegree
        KL_InMixed=KL_InMixed+DistInSmoothed(x,y)*log(DistInSmoothed(x,y)/DistMixed(x,y));
        KL_OutMixed=KL_OutMixed+DistOutSmoothed(x,y)*log(DistOutSmoothed(x,y)/DistMixed(x,y));
        JSD = 0.5*(KL_InMixed + KL_OutMixed);
    end
end
