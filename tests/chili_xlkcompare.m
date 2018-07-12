function [err,data] = test(opt,olddata)

%===============================================================================
% Compare X^L_0K coefficients from chili_xlk and chili_xlmk
%===============================================================================
Potential.lambda = [1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8];
Potential.L = [2 2 2 4 4 4 4 4];
Potential.M = [0 0 0 0 0 0 0 0];
Potential.K = [0 1 2 0 1 2 3 4];

R = [1 2 3];
maxLx = max(Potential.L)*2;

[XLMK{1:maxLx+1}] = runprivate('chili_xlmk',Potential,R);
XLK = runprivate('chili_xlk',Potential,R);

for k = numel(XLMK):-1:1
  X_(k,1:2*(k-1)+1) = XLMK{k}(k,:);
end

err = max(max(abs(X_-XLK)))>1e-13;
data = [];
