function K = Degrees(A)
%This calculates the degrees of all nodes in the multilayer network
[~,T,N] = size(A);
K = reshape(sum(A),[T,N]);
end