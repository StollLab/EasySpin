function [err,data] = test(opt,olddata)
%orbital angular momenta can also be defined as spins, therefore the two
%Hamiltonians should be identical
rng_(5,'twister');

n = 3;%randi_(2);
Sys.S = 1/2*ones(1,3);%randi_(3,1,n)/2;
Sys.g = rand(3*n,3);
Sys.L = 1*ones(1,3);%randi_(2,1,n);
Sys.soc = rand(n,2)*1000;
Sys.orf = rand(n,1);
lenS = length(Sys.S);
kS =0;
if n>1
  kS = nchoosek(lenS,2);
  coup = nchoosek(1:lenS,2);
  Sys.ee = rand(kS,1);
  Sys.ee2 = rand(kS,1);
end

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
  if m<=kS
    y = logical((k(:,1)==coup(m,1).*(k(:,2)==coup(m,2))));
    PureSpin.ee(y) = Sys.ee(m);
    PureSpin.ee2(y) = Sys.ee2(m);
  end
end

%build Zero-Field splitting part
for k=2:2:8
  lfieldname = sprintf('CF%d',k);
  sfieldname = sprintf('B%d',k);
  Sys.(sfieldname) = rand(n,2*k+1).*repmat(((k/2)<=Sys.S).',1,2*k+1);
  Sys.(lfieldname) = rand(n,2*k+1).*repmat(((k/2)<=Sys.L).',1,2*k+1);
  PureSpin.(sfieldname) = [Sys.(sfieldname);Sys.(lfieldname)];
end

[H1,GX1,GY1,GZ1] = sham(PureSpin);
[H2,GX2,GY2,GZ2] = sham(Sys);

% comparison including nuclear spins will be done on basis of the
% eigenvalues, as the order in the spin-vector differ between the two
% cases, in Sys we have |m_S,m_I,m_L> and in PureSpin |m_S, m_L, m_I>
Sys.Nucs = '1H,1H';
Sys.A = rand(3*2,3*n);
PureSpin.Nucs = Sys.Nucs;
PureSpin.A = [Sys.A, zeros(3*2,3*n)];
field = rand(1,3)*1e3;
E1 = eig(sham(PureSpin,field));
E2 = eig(sham(Sys,field));

err = ~all([areequal(H1,H2,1e-10),areequal(GX1,GX2,1e-10),...
  areequal(GY1,GY2,1e-10),areequal(GZ1,GZ2,1e-10),areequal(E1,E2,1e-6)]);
data = [];