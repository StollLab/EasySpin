function [err,data] = test(opt,olddata)
%orbital angular momenta can also be defined as spins, therefore the 
%results should be identical for the two cases
rng_(5,'twister');
n = randi_(2);
Sys.S = randi_(3,1,n)/2;
Sys.g = rand(3*n,3);
Sys.L = randi_(2,1,n);
Sys.soc = rand(n,2)*1000;
Sys.orf = rand(n,1);
lenS = length(Sys.S);
if n>1, Sys.ee = zeros(nchoosek(lenS,2),1);end

PureSpin.S = [Sys.S,Sys.L];
%build array of g matrices:
PureSpin.g = [Sys.g;zeros(3*n,3)];
for k=1:n
  PureSpin.g(3*(n+k-1)+1:3*(n+k),:) = -diag(Sys.orf(k)*ones(1,3));
end

%distribute soc over ee and ee2
len = 2*lenS;
k = nchoosek(1:len,2);
eelen = nchoosek(len,2);
PureSpin.ee = zeros(eelen,1);
PureSpin.ee2 = zeros(eelen,1);
for m=1:lenS
  x = logical((k(:,1)==m).*(k(:,2)==m+lenS));
  PureSpin.ee(x) = Sys.orf(m)*Sys.soc(m,1);
  PureSpin.ee2(x) = Sys.orf(m)*Sys.orf(m)*Sys.soc(m,2);
end

%build Zero-Field splitting part
for k=2:2:8
  lfieldname = sprintf('CF%d',k);
  sfieldname = sprintf('B%d',k);
  Sys.(sfieldname) = rand(n,2*k+1).*repmat(((k/2)<=Sys.S).',1,2*k+1);
  Sys.(lfieldname) = rand(n,2*k+1).*repmat(((k/2)<=Sys.L).',1,2*k+1);
  PureSpin.(sfieldname) = [Sys.(sfieldname);Sys.(lfieldname)];
end

%create random experiment
Exp.Temperature = rand(1,50)*100;
Exp.Field = rand(1,50)*100;

[m1,chi1] = curry(Sys,Exp);
[m2,chi2] = curry(PureSpin,Exp);

thr = 1e-6;
err = ~areequal(m1,m2,thr,'rel') || ~areequal(chi1,chi2,thr,'rel');
data = [];