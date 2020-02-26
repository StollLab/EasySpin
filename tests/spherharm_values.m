function [ok,data] = test(opt,olddata)

% Regression test on values
%=======================================================

L = [1 1 2 2 2  3 4  5 20];
m = [0 1 -1 0 2 2 1 -3 -7];

for k = 1:numel(L)
  v(k) = spherharm(L(k),m(k),1,1);
end

data.v = v;

if ~isempty(olddata)
  ok = all(abs(olddata.v-v)<1e-10);
else
  ok = [];
end

