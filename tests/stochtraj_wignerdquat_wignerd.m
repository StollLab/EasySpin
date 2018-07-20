function [err,data] = test(opt,olddata)

%-------------------------------------------------------------------------------
% Assert that wignerd() and stochtraj_wignerdquat() and euler2quat() are
% consistent.
%-------------------------------------------------------------------------------

rng(1);
a = rand*2*pi;
b = rand*pi;
c = rand*2*pi;

L = 4;

q = euler2quat(a,b,c);

idxM = 0;
for M = -L:L
  idxM = idxM+1;
  idxK = 0;
  for K = -L:L
    idxK = idxK+1;
    D1(idxM,idxK) = wignerd([L M K],a,b,c);    
    D2(idxM,idxK) = runprivate('stochtraj_wignerdquat',L,M,K,q);
  end
end

err = max(abs(D1(:)-D2(:))) > 1e-15;

data = [];

end
