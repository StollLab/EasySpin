function ok = ham_hf_deriv_value()

% Test whether hyperfine hamiltonian derivatives are correct,
% for a system with one electron and one nucleus

Sys.S = 1/2;
Sys.Nucs = '14N';
Sys.A = [1.435 2.765 3.9876];

[Hhf,dHhf] = ham_hf(Sys);

% Calculate hamiltonian from derivatives
Hhf_2 = Sys.A(1)*dHhf{1}{1} + Sys.A(2)*dHhf{1}{2} + Sys.A(3)*dHhf{1}{3};

ok(1) = areequal(Hhf_2,Hhf,1e-10,'abs');

% test for a larger spin system and different euler angles
Sys.S = [1/2 1];
Sys.Nucs = '14N,1H';
Sys.A = [1.435 2.765 3.9876, 1 2 3 ; 4 5 6, 7 8 9];
Sys.AFrame = [35 26 0, 87 63 42; 81 23 60, 42 36 47]*pi/180;
Sys.J=1;

[Hhf,dHhf] = ham_hf(Sys);

% Calculate hamiltonian from derivatives
Hhf_2 = Sys.A(1,1)*dHhf{1,1}{1} + Sys.A(1,2)*dHhf{1,1}{2} + Sys.A(1,3)*dHhf{1,1}{3}+...
        Sys.A(2,1)*dHhf{1,2}{1} + Sys.A(2,2)*dHhf{1,2}{2} + Sys.A(2,3)*dHhf{1,2}{3}+...
        Sys.A(1,4)*dHhf{2,1}{1} + Sys.A(1,5)*dHhf{2,1}{2} + Sys.A(1,6)*dHhf{2,1}{3}+...
        Sys.A(2,4)*dHhf{2,2}{1} + Sys.A(2,5)*dHhf{2,2}{2} + Sys.A(2,6)*dHhf{2,2}{3};

ok(2) = areequal(Hhf_2,Hhf,1e-10,'abs');

ok=all(ok);