function [JSD] = CalculateJSD5(DistIn,DistOut)
[MaxDegree,~] = size(DistIn);
%Smooth using Dirichlet prior
DistInSmoothed = DirichletPrior5D(DistIn);
DistOutSmoothed = DirichletPrior5D(DistOut);
DistMixed = 0.5*(DistOutSmoothed+DistInSmoothed);
KL_InMixed=0;
KL_OutMixed=0;
for x =1:MaxDegree
    for y=1:MaxDegree
        for z=1:MaxDegree
            for w=1:MaxDegree
                for u=1:MaxDegree
                    KL_InMixed=KL_InMixed+DistInSmoothed(x,y,z,w,u)*log(DistInSmoothed(x,y,z,w,u)/DistMixed(x,y,z,w,u));
                    KL_OutMixed=KL_OutMixed+DistOutSmoothed(x,y,z,w,u)*log(DistOutSmoothed(x,y,z,w,u)/DistMixed(x,y,z,w,u));
                    JSD = 0.5*(KL_InMixed + KL_OutMixed);
                end
            end
        end
    end
end
