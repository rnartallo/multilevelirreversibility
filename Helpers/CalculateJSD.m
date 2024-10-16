function [JSD] = CalculateJSD(DistIn,DistOut)
[MaxIn,N] = size(DistIn);
[MaxOut,~] = size(DistOut);
%Pad to the same length
MaxDegree = max(MaxIn,MaxOut);
DistInPadded = zeros([MaxDegree,N]);
DistOutPadded = zeros([MaxDegree,N]);
%Rescale for same support
DistInPadded(1:MaxIn,:)=DistIn*MaxIn/MaxDegree;
DistOutPadded(1:MaxOut,:)=DistOut*MaxOut/MaxDegree;
%Smooth using Dirichlet prior
DistInSmoothed = DirichletPrior1D(DistInPadded);
DistOutSmoothed = DirichletPrior1D(DistOutPadded);
DistMixed = 0.5*(DistOutSmoothed+DistInSmoothed);
KL_InMixed=0;
KL_OutMixed=0;
for n=1:N
    for x =1:MaxDegree
        KL_InMixed=KL_InMixed+DistInSmoothed(x,n)*log(DistInSmoothed(x,n)/DistMixed(x,n));
        KL_OutMixed=KL_OutMixed+DistOutSmoothed(x,n)*log(DistOutSmoothed(x,n)/DistMixed(x,n));
        JSD(n) = 0.5*(KL_InMixed + KL_OutMixed);
    end
end
