function ok = test(opt)

% Spin-correlated radical pair: shortcut input
%===============================================

% Spin system
Sys = struct('S',[1/2 1/2],'g',[2.00 1.99],'HStrain',1);
Sys.J = -3; % exchange coupling in MHz
Sys.dip = [1 1 -2]*5; % electron-electron coupling in MHz

% Experimental parameters
Exp = struct('mwFreq',9.5,'Range',[337 343],'Harmonic',0);

% User-specified population vector
S = 1/sqrt(2)*[0; 1; -1; 0];
Sys.initState = S*S';

% Simulation options
Opt = struct;

[x,spc] = pepper(Sys,Exp,Opt);

% Initial state defined through shortcut
Sys.initState = 'singlet';

[x,spcshortcut] = pepper(Sys,Exp,Opt);

if opt.Display
    plot(x,spc,x,spcshortcut);
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: Spin-correlated radical pair');
end

ok = areequal(spc,spcshortcut,1e-4,'rel');
