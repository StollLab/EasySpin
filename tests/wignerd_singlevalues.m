function [err,data] = test(opt,olddata)

%============================================================================
% Compare D matrix calculated with a single call to wignerd() to D matrix
% compiled from single-element calls to wignerd()
%============================================================================

rng_(15654,'twister');

Jlist = [0:0.5:6 7:10 40];

for iJ = 1:numel(Jlist)
  
  J = Jlist(iJ);
  alpha = rand*2*pi;
  beta = rand*pi;
  gamma = rand*2*pi;
  
  % Full matrix in one call
  D1 = wignerd(J,alpha,beta,gamma);
  
  % Full matrix build up with single-matrix-element calls
  D2 = zeros(2*J+1);
  for m1 = J:-1:-J
    for m2 = J:-1:-J
      D2(J-m1+1,J-m2+1) = wignerd([J m1 m2],alpha,beta,gamma);
    end
  end
  
  err(iJ) = ~areequal(D1,D2,1e-13);
  
end

data = [];
