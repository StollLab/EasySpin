function ok = test()

% Test whether hyperfine hamiltonian derivatives are correct,
% for a system with one electron and one nucleus

Sys.S = 1/2;
Sys.Nucs = '14N';
Sys.A = [1.435 2.765 3.9876];

[Hhf,dHhf] = ham_hf(Sys);

% Calculate hamiltonian from derivatives
Hhf_2 = Sys.A(1)*dHhf{1}{1} + Sys.A(2)*dHhf{1}{2} + Sys.A(3)*dHhf{1}{3};

ok = areequal(Hhf_2,Hhf,1e-10,'abs');

