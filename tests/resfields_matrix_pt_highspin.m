function ok = test()

% For a high-spin system, compare line positions and intensities
% as obtained by matrix diagonalization, eigenfields and perturbation theory.

clear Sys Exp
Sys.S = 5/2;
Sys.D = 10;
Sys.lwpp = 0.3;

Exp.mwFreq = 9.5;
Exp.Range = [300 380];
Exp.SampleFrame = [pi/4 pi/4 0];

[Ba,Ia] = resfields_perturb(Sys,Exp);
[Bb,Ib] = resfields(Sys,Exp);
[Bc,Ic] = resfields_eig(Sys,Exp);

[Ba,idx] = sort(Ba); Ia = Ia(idx);
[Bb,idx] = sort(Bb); Ib = Ib(idx);
[Bc,idx] = sort(Bc); Ic = Ic(idx);

ok(1) = all(abs(Ba-Bb)<0.0001);
ok(2) = all(abs(Ia-Ib)<0.01*max(Ib));
ok(3) = numel(Bc)==numel(Bb) && all(abs(Bc(:)-Bb(:))<0.0001);
ok(4) = all(abs(Ic(:)-Ib(:))<0.01*max(Ib));
