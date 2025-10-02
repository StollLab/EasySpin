function ok = ham_ee_deriv_value()

% Test whether electron-electron hamiltonian derivatives are correct,
% for a system with 2 electrons

Sys.S = [1/2 1/2];
J=10;
Sys.ee = [1 1 -2]+ J;

[Hee,dHee] = ham_ee(Sys);

% Calculate hamiltonian from derivatives
Hee_2 = Sys.ee(1)*dHee{1}{1} + Sys.ee(2)*dHee{1}{2} + Sys.ee(3)*dHee{1}{3};

ok(1) = areequal(Hee_2,Hee,1e-10,'abs');

clear Sys
% test for a larger spin system and bilinear term
Sys.S = [1/2 1 1/2];
dip=[[1 1 -2];2*[1 1 -2];3*[1 1 -2]];
J=[1 2 3];
Sys.ee=dip+J';
Sys.eeFrame=[10 23 65; 78 42 136; -52 86 276]*pi/180;
Sys.ee2=[100 200 300];
[Hee,dHee] = ham_ee(Sys);

% Calculate hamiltonian from derivatives
Hee_2 = Sys.ee(1,1)*dHee{1}{1} + Sys.ee(1,2)*dHee{1}{2} + Sys.ee(1,3)*dHee{1}{3}+...
         Sys.ee(2,1)*dHee{2}{1} + Sys.ee(2,2)*dHee{2}{2} + Sys.ee(2,3)*dHee{2}{3}+...
         Sys.ee(3,1)*dHee{3}{1} + Sys.ee(3,2)*dHee{3}{2} + Sys.ee(3,3)*dHee{3}{3};

Hee_2 = Hee_2 + Sys.ee2(1)*dHee{1}{4} + Sys.ee2(2)*dHee{2}{4} + Sys.ee2(3)*dHee{3}{4};

ok(2) = areequal(Hee_2,Hee,1e-10,'abs');
