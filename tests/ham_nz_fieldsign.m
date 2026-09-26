function ok = test()

% Assert that the nuclear Zeeman Hamiltonian is H = -mu.B

Sys.Nucs = '1H';
Sys.A = 1;
B = [100 -200 300];

[mux,muy,muz] = ham_nz(Sys);
H = ham_nz(Sys,B);

ok = areequal(H,-(mux*B(1)+muy*B(2)+muz*B(3)),1e-10,'abs');
