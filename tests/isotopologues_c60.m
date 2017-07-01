function [err,data] = test(opt,olddata)

% Isotopologues of C60 with natural-abundance carbon
%-------------------------------------------------------------------------------
n = 60;
threshold = 0;
isoList = isotopologues('C',n,[],threshold);

% Check numbe rof equivalent nuclei (this assumes isotopologues are not
% sorted by abundance)
n_calc = reshape([isoList(2:n).n],2,[]).';
n_expected = [n-1:-1:1; 1:n-1]';
ok = all(n_calc(:)==n_expected(:));
ok = ok && isoList(1).n==n;
ok = ok && isoList(end).n==n;

% Check abundances
calculatedWeights = [isoList.weight];
ab_c13 = nucabund('13C');
ab_c12 = nucabund('12C');
precision = 1e-7;
ok = ok && ...
  areequal(calculatedWeights(1),ab_c12^n,precision) && ...
  areequal(calculatedWeights(2),ab_c12^(n-1)*ab_c13*60,precision) && ...;
  areequal(calculatedWeights(3),ab_c12^(n-2)*ab_c13^2*60*59/2,precision);
  
err = ~ok;
data = [];
