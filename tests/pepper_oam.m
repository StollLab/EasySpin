function [err,data] = test(opt,olddata)
%orbital angular momenta can also be defined as spins, therefore the two
%Hamiltonians should be identical
rand('twister',5);

n = 1;%randi_(2);
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

%build random experiment for frequency-domain pepper
FDExp.Temperature = rand * 300;
FDExp.Field = rand *1e3;


[nu,fd1] = pepper(Sys,FDExp);
fd2 = pepper(PureSpin,FDExp);


%build experimet based on FD sim, always a transition in spectral window 
[ignore, ind] = max(fd1);
Exp.mwFreq = nu(ind);
Exp.CenterSweep = FDExp.Field*[1 0.5];
Exp.Temperature = FDExp.Temperature;

Opt = struct();

s1 = pepper(Sys,Exp,Opt);
s2 = pepper(PureSpin,Exp,Opt);

Sys.Nucs = '1H,14N';
Sys.A = rand(3*2,3*n);
PureSpin.Nucs = Sys.Nucs;
PureSpin.A = [Sys.A, zeros(3*2,3*n)];
Opt.Method = 'hybrid';

fd3 = pepper(Sys,FDExp,Opt);
fd4 = pepper(PureSpin,FDExp,Opt);

s3 = pepper(Sys,Exp,Opt);
s4 = pepper(PureSpin,Exp,Opt);

err = ~all([areequal(fd1,fd2,1e-10),areequal(fd3,fd4,1e-6),...
  areequal(s1,s2,1e-10),areequal(s3,s4,1e-6)]);
data = [];
