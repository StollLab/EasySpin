function ok = ham_zf_deriv_value()

% Test whether ZFS hamiltonian derivatives are correct,
% for a system with 1 electron S = 1

Sys.S = 1;
Sys.D = 100*[-1 -1 2];

[Hzf,dHzf] = ham_zf(Sys);

% Calculate hamiltonian from derivatives
Hzf_2 = Sys.D(1)*dHzf{1}{1} + Sys.D(2)*dHzf{1}{2} + Sys.D(3)*dHzf{1}{3};

ok(1) = areequal(Hzf_2,Hzf,1e-10,'abs');


% test for a larger spin system and asymmetry and euler angles
clear Sys
Sys.S = [3/2 5/2];
D = [1 2]';
E = [0.5 0.3]';
Sys.D = D*[-1,-1,2]/3 + E*[1,-1,0];
Sys.DFrame = [23 96 30 ; 10 6 81]*pi/180;
Sys.J=0;

[Hzf,dHzf] = ham_zf(Sys);

% Calculate hamiltonian from derivatives
Hzf_2 = Sys.D(1,1)*dHzf{1}{1} + Sys.D(1,2)*dHzf{1}{2} + Sys.D(1,3)*dHzf{1}{3}+...
        Sys.D(2,1)*dHzf{2}{1} + Sys.D(2,2)*dHzf{2}{2} + Sys.D(2,3)*dHzf{2}{3};

ok(2) = areequal(Hzf_2,Hzf,1e-10,'abs');
