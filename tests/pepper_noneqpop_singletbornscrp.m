function ok = test(opt)

% Spin-correlated radical pair: different types of input
%==========================================================

% Spin system
Sys = struct('S',[1/2 1/2],'g',[2.00 1.99],'HStrain',1);
Sys.J = -3; % exchange coupling in MHz
Sys.dip = [1 1 -2]*5; % electron-electron coupling in MHz

% Experimental parameters
Exp = struct('mwFreq',9.5,'Range',[337 343],'Harmonic',0);

% Simulation options
Opt = struct;

% Density matrix in uncoupled basis
S = 1/sqrt(2)*[0; 1; -1; 0];
Sys.initState = S*S';
[~,spc1] = pepper(Sys,Exp,Opt);

% Density matrix in coupled basis
S = [0 0 0 1].';
Sys.initState = {S*S','coupled'};
[~,spc2] = pepper(Sys,Exp,Opt);

% Population vector in coupled basis
Sys.initState = {[0 0 0 1],'coupled'};
[~,spc3] = pepper(Sys,Exp,Opt);

% Initial state defined through shortcut
Sys.initState = 'singlet';
[x,spc4] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(x,spc1,x,spc2,x,spc3,x,spc4);
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: Spin-correlated radical pair');
end

ok = areequal(spc1,spc2,1e-12,'abs') && areequal(spc1,spc3,1e-12,'abs') && areequal(spc1,spc4,1e-12,'abs');
