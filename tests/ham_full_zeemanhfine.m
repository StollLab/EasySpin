function ok = test()

% ham vs ham_ez/ham_nz/ham_hf

Sys = struct('S',1/2,'g',[2 3 4],'Nucs','63Cu','A',[50 50 350]);

B = rand(1,3)*400;

% Construct Hamiltonian using ham()
[F1,mu1] = ham(Sys,B);
H1 = F1 - norm(B)*mu1;

% Construct Hamiltonian using ham_*()
F = ham_nq(Sys) + ham_hf(Sys);
[muxe,muye,muze] = ham_ez(Sys);
[muxn,muyn,muzn] = ham_nz(Sys);
H = F - B(1)*(muxe+muxn) - B(2)*(muye+muyn) - B(3)*(muze+muzn);

ok = areequal(H1,H,1e-10,'rel');
