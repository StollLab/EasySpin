function ok = test()

% Orbital angular momenta can also be defined as spins, therefore the two
% Hamiltonians should be identical

rng(5);

n = 2; % number of centers (each has spin and orbital angular momentum)

Sys.S = 1/2*ones(1,n);
Sys.g = rand(3*n,3);
Sys.L = 1*ones(1,n);
Sys.soc = rand(n,2)*1000;
Sys.gL = rand(n,1);
lenS = length(Sys.S);
kS = 0;
if n>1
  kS = nchoosek(lenS,2);
  coup = nchoosek(1:lenS,2);
  Sys.ee = rand(kS,1);
  Sys.ee2 = rand(kS,1);
end

% Build zero-field splitting part
for k=2:2:8
  lfieldname = sprintf('CF%d',k);
  sfieldname = sprintf('B%d',k);
  Sys.(sfieldname) = rand(n,2*k+1).*repmat(((k/2)<=Sys.S).',1,2*k+1);
  Sys.(lfieldname) = rand(n,2*k+1).*repmat(((k/2)<=Sys.L).',1,2*k+1);
  pureSpin.(sfieldname) = [Sys.(sfieldname);Sys.(lfieldname)];
end

% Build pure-spin Hamiltonian (i.e. model the OAMs at electron spins)
pureSpin.S = [Sys.S,Sys.L];
pureSpin.g = [Sys.g; zeros(3*n,3)];
for k = 1:n
  pureSpin.g(3*(n+k-1)+1:3*(n+k),:) = Sys.gL(k)*eye(3);
end

% Add Sys.soc as pureSpin.ee and pureSpin.ee2
len = 2*lenS;
k = nchoosek(1:len,2);
eelen = nchoosek(len,2);
pureSpin.ee = zeros(eelen,1);
pureSpin.ee2 = zeros(eelen,1);
for m = 1:lenS
  x = logical((k(:,1)==m).*(k(:,2)==m+lenS));
  pureSpin.ee(x) = Sys.soc(m,1);
  pureSpin.ee2(x) = Sys.soc(m,2);
  if m<=kS
    y = logical((k(:,1)==coup(m,1).*(k(:,2)==coup(m,2))));
    pureSpin.ee(y) = Sys.ee(m);
    pureSpin.ee2(y) = Sys.ee2(m);
  end
end

[H01,muX1,muY1,muZ1] = ham(pureSpin);
[H02,muX2,muY2,muZ2] = ham(Sys);

% comparison including nuclear spins will be done on basis of the
% eigenvalues, as the order in the spin-vector differ between the two
% cases, in Sys we have |m_S,m_I,m_L> and in pureSpin |m_S, m_L, m_I>
Sys.Nucs = '1H,1H';
Sys.A = rand(3*2,3*n);
pureSpin.Nucs = Sys.Nucs;
pureSpin.A = [Sys.A, zeros(3*2,3*n)];
field = rand(1,3)*1e3;
E1 = eig(ham(pureSpin,field));
E2 = eig(ham(Sys,field));

ok(1) = areequal(H01,H02,1e-10,'abs');
ok(2) = areequal(muX1,muX2,1e-10,'abs');
ok(3) = areequal(muY1,muY2,1e-10,'abs');
ok(4) = areequal(muZ1,muZ2,1e-10,'abs');
ok(5) = areequal(E1,E2,1e-6,'abs');
