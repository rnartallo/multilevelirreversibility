function [JSD] = CalculateJSD3(DistIn,DistOut)
[MaxDegree,~] = size(DistIn);
%Smooth using Dirichlet prior
DistInSmoothed = DirichletPrior3D(DistIn);
DistOutSmoothed = DirichletPrior3D(DistOut);
DistMixed = 0.5*(DistOutSmoothed+DistInSmoothed);
KL_InMixed=0;
KL_OutMixed=0;
for x =1:MaxDegree
    for y=1:MaxDegree
        for z=1:MaxDegree
            KL_InMixed=KL_InMixed+DistInSmoothed(x,y,z)*log(DistInSmoothed(x,y,z)/DistMixed(x,y,z));
            KL_OutMixed=KL_OutMixed+DistOutSmoothed(x,y,z)*log(DistOutSmoothed(x,y,z)/DistMixed(x,y,z));
            JSD = 0.5*(KL_InMixed + KL_OutMixed);
        end
    end
end
