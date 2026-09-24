function ok = test()

% Compare principal values read from the main ORCA 6 output files with
% reference values copied manually from these files:
%   triplet formaldehyde: g(tot) and diagonalized D matrix (cm^-1)
%   nitroxide, 14N: A(Tot) (MHz), e**2qQ (MHz), and eta

folder = 'orca/v6.1.0/';

g_ref = [2.0033261 2.0046919 2.0057848];
D_ref = [-0.319688 -0.182797 0.304358];  % cm^-1
A_ref = [-0.8057 -1.2281 73.4439];  % MHz
e2qQ_ref = -4.738579;  % MHz
eta_ref = 0.298034;
Qtol = 1e-4;  % MHz

iN = 4;  % index of nitrogen atom in nitroxide
cm2MHz = 100*clight/1e6;

Sys = orca2easyspin([folder 'tripletformaldehyde.out']);
ok(1) = areequal(sort(Sys.g),g_ref,1e-6,'abs');
ok(2) = areequal(sort(Sys.D),D_ref*cm2MHz,0.1,'abs');

Sys = orca2easyspin([folder 'nitroxide.out']);
n = find(Sys.NucsIdx==iN);
ok(3) = areequal(sort(Sys.A(n,:)),sort(A_ref),1e-3,'abs');
K = e2qQ_ref/4;  % I = 1
Q_ref = K*[-(1-eta_ref),-(1+eta_ref),2];
ok(4) = areequal(sort(Sys.Q(n,:)),sort(Q_ref),Qtol,'abs');
