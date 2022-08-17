function ok = test(opt)

% Photoexcited triplet state coupled to a nuclear spin
% Check that input of electron state populations matches
% result for input of electron and nuclear state
% populations
%=======================================================

% Spin system and experimental parameters
Sys = struct('S',1,'g',2,'lw',0.3,'D',[-800 100],'Nucs','14N','A',[10 10 60]);
Exp = struct('mwFreq',9.5,'Range',[290 390],'Harmonic',0);

% User-specified population vector for each electron manifold
Sys.initState = {[0.1 0.5 0.4],'zerofield'};

% Simulation options
Opt = struct;

[x,spc] = pepper(Sys,Exp,Opt);

% User-specified population vector for all sublevels
Sys.initState = {[0.1 0.1 0.1 0.5 0.5 0.5 0.4 0.4 0.4],'zerofield'};

[x,spcshortcut] = pepper(Sys,Exp,Opt);

if opt.Display
    plot(x,spc,x,spcshortcut);
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: Spin-correlated radical pair');
end

ok = areequal(spc,spcshortcut,1e-4,'rel');
