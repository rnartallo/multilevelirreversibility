function [JSD] = CalculateJSD4(DistIn,DistOut)
[MaxDegree,~] = size(DistIn);
%Smooth using Dirichlet prior
DistInSmoothed = DirichletPrior4D(DistIn);
DistOutSmoothed = DirichletPrior4D(DistOut);
DistMixed = 0.5*(DistOutSmoothed+DistInSmoothed);
KL_InMixed=0;
KL_OutMixed=0;
for x =1:MaxDegree
    for y=1:MaxDegree
        for z=1:MaxDegree
            for w=1:MaxDegree
                KL_InMixed=KL_InMixed+DistInSmoothed(x,y,z,w)*log(DistInSmoothed(x,y,z,w)/DistMixed(x,y,z,w));
                KL_OutMixed=KL_OutMixed+DistOutSmoothed(x,y,z,w)*log(DistOutSmoothed(x,y,z,w)/DistMixed(x,y,z,w));
                JSD = 0.5*(KL_InMixed + KL_OutMixed);
            end
        end
    end
end
