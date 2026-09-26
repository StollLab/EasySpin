function ok = test()

% Assert that the magnetic-moment operators returned by ham for a spin system
% with higher-order Zeeman terms are consistent with the full Hamiltonian

rng(5);

B = [1 2 3];
Sys.S = 1/2;
Sys.Ham110 = rand;
Sys.Ham112 = rand(1,5);

% Highest order in B is 1
[H0,mx,my,mz] = ham(Sys);
Hfull = ham(Sys,B);
ok(1) = areequal(Hfull,H0-B(1)*mx-B(2)*my-B(3)*mz,1e-10,'abs');

[H0,muz] = ham(Sys,[0 0 B(3)]);
ok(2) = areequal(ham(Sys,[0 0 B(3)]),H0-B(3)*muz,1e-10,'abs');

% Highest order in B is 3
Sys.Ham312 = rand(1,5);
[H0,mu,~,~] = ham(Sys);
Hfull = ham(Sys,B);
Href = H0 - B(1)*mu{1} - B(2)*mu{2} - B(3)*mu{3} + ham_ezho(Sys,B,[],'',2:3);
ok(3) = areequal(Hfull,Href,1e-10,'abs');
