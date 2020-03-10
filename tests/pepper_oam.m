function ok = test()
%orbital angular momenta can also be defined as spins, therefore the two
%Hamiltonians should be identical

rng(5);

% Build spin system with spin in Sys.S and orbital angular momenta in Sys.L
n = 1;
SysSL.S = randi(3,1,n)/2;
SysSL.g = rand(3*n,3);
SysSL.L = randi(2,1,n);
SysSL.soc = rand(n,2)*1000;
SysSL.orf = rand(n,1);
nSpins = length(SysSL.S);
if n > 1
  SysSL.ee = zeros(nchoosek(nSpins,2),1);
end

% Build spin system with both spin and orbital angular momenta in Sys.S
SysS.S = [SysSL.S,SysSL.L];
% Build array of g matrices:
SysS.g = [SysSL.g;zeros(3*n,3)];
for k = 1:n
  SysS.g(3*(n+k-1)+1:3*(n+k),:) = -diag(SysSL.orf(k)*ones(1,3));
end

% Distribute soc over ee and ee2
len = 2*nSpins;
k = nchoosek(1:len,2);
eelen = nchoosek(len,2);
SysS.ee = zeros(eelen,1);
SysS.ee2 = zeros(eelen,1);
for m = 1:nSpins
  x = logical((k(:,1)==m).*(k(:,2)==m+nSpins));
  SysS.ee(x) = SysSL.orf(m)*SysSL.soc(m,1);
  SysS.ee2(x) = SysSL.orf(m)*SysSL.orf(m)*SysSL.soc(m,2);
end

% Build zero-field splitting part
for k = 2:2:8
  lfieldname = sprintf('CF%d',k);
  sfieldname = sprintf('B%d',k);
  SysSL.(sfieldname) = rand(n,2*k+1).*repmat(((k/2)<=SysSL.S).',1,2*k+1);
  SysSL.(lfieldname) = rand(n,2*k+1).*repmat(((k/2)<=SysSL.L).',1,2*k+1);
  SysS.(sfieldname) = [SysSL.(sfieldname);SysSL.(lfieldname)];
end

% Build random experiment for frequency-domain pepper
FDExp.Temperature = rand * 300;
FDExp.Field = rand *1e3;

% Compare S&L with S-only spin system
[nu,fd1] = pepper(SysSL,FDExp);
fd2 = pepper(SysS,FDExp);
ok(1) = areequal(fd1,fd2,1e-8,'abs');

% Build field-sweep experimet based on FD sim, always a transition in spectral window 
[~, ind] = max(fd1);
Exp.mwFreq = nu(ind);
Exp.CenterSweep = FDExp.Field*[1 0.5];
Exp.Temperature = FDExp.Temperature;

Opt = struct;

s1 = pepper(SysSL,Exp,Opt);
s2 = pepper(SysS,Exp,Opt);
ok(2) = areequal(s1,s2,1e-10,'abs');

% Test with an added nucleus
%-------------------------------------------------------------------------------
SysSL.Nucs = '1H';
SysSL.A = rand(3,3*n);
SysS.Nucs = SysSL.Nucs;
SysS.A = [SysSL.A, zeros(3,3*n)];
Opt.Method = 'hybrid';

% Frequency sweep
fd3 = pepper(SysSL,FDExp,Opt);
fd4 = pepper(SysS,FDExp,Opt);
ok(3) = areequal(fd3,fd4,1e-6,'abs');

% Field sweep
s3 = pepper(SysSL,Exp,Opt);
s4 = pepper(SysS,Exp,Opt);
ok(4) = areequal(s3,s4,1e-6,'abs');
