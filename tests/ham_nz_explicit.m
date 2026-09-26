function ok = test()

% Assert that the nuclear Zeeman Hamiltonian is H = -mu.B with
% mu = +gn*muN*I (in MHz), by comparing it to an explicitly constructed matrix

B0 = 350;  % mT
isotope = '14N';
phi = deg2rad(47);
theta = deg2rad(22);

zL = ang2vec(phi,theta);
B = B0*zL;

% Reference Hamiltonian
gn = nucgval(isotope);
I = nucspin(isotope);
[Ix,Iy,Iz] = sop([1/2 I],'x2','y2','z2');
mux = +gn*nmagn*Ix;  % J/T
muy = +gn*nmagn*Iy;  % J/T
muz = +gn*nmagn*Iz;  % J/T
Href = -(B(1)*mux+B(2)*muy+B(3)*muz)*1e-3;  % J, includes mT -> T
Href = Href/planck/1e6;  % J -> MHz

% Hamiltonian from ham_nz
Sys.S = 1/2;
Sys.Nucs = isotope;
Sys.A = 1;  % MHz
H = ham_nz(Sys,B);

ok = areequal(H,Href,1e-12,'abs');
