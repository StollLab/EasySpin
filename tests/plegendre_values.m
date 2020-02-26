function [ok,data] = test(opt,olddata)

% Regression test: Values of associated Legendre polynomials

n = [1 1 2 2 2  3 4  5 20];
m = [0 1 -1 0 2 2 1 -3 -7];

for k = 1:numel(n)
  v(k) = plegendre(n(k),m(k),0.5154);
end

data.v = v;

if ~isempty(olddata)
  ok = areequal(v,olddata.v,1e-10,'abs');
else
  ok = [];
end

