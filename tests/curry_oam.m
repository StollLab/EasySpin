function ok = test()

% orbital angular momenta can also be defined as spins, therefore the 
% results should be identical for the two cases

rng(5);
n = 1;
SysSL.S = 1/2;
SysSL.g = diag([2 2.1 2.2]);
SysSL.L = 1;
SysSL.soc = rand(1,2)*1000;
SysSL.gL = rand;
lenS = length(SysSL.S);

PureSpin.S = [SysSL.S,SysSL.L];
% build array of g matrices
PureSpin.g = [SysSL.g;zeros(3*n,3)];
for k = 1:n
  gL = diag(SysSL.gL(k)*ones(1,3));
  PureSpin.g(3*(n+k-1)+1:3*(n+k),:) = gL;
end

% distribute soc over ee and ee2
len = 2*lenS;
k = nchoosek(1:len,2);
eelen = nchoosek(len,2);
PureSpin.ee = zeros(eelen,1);
PureSpin.ee2 = zeros(eelen,1);
for m = 1:lenS
  x = (k(:,1)==m) & (k(:,2)==m+lenS);
  PureSpin.ee(x) = SysSL.soc(m,1);
  PureSpin.ee2(x) = SysSL.soc(m,2);
end

% build zero-field splitting part
for k = 2:2:8
  lfieldname = sprintf('CF%d',k);
  sfieldname = sprintf('B%d',k);
  SysSL.(sfieldname) = rand(n,2*k+1).*repmat(((k/2)<=SysSL.S).',1,2*k+1);
  SysSL.(lfieldname) = rand(n,2*k+1).*repmat(((k/2)<=SysSL.L).',1,2*k+1);
  PureSpin.(sfieldname) = [SysSL.(sfieldname);SysSL.(lfieldname)];
end

Exp.Temperature = [0.5 1 2 4 10 20 50 100 200 500];
Exp.Field = [1 10 100 1000 10000];

[m1,chi1] = curry(SysSL,Exp);
[m2,chi2] = curry(PureSpin,Exp);

thr = 1e-6;
ok = areequal(m1,m2,thr,'rel') && areequal(chi1,chi2,thr,'rel');
