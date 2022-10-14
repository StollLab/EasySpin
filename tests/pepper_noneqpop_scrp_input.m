function ok = test(opt)

% Check equivalence of different inputs for a spin-correlated radical pair coupled
% to a 14N nuclear spin formed from a thermalised triplet state

% Spin system
Sys.S = [1/2 1/2];
Sys.g = [2.0027; 2.0000];
Sys.J = -6; % MHz

Sys.Nucs = '14N';
Sys.A = [0 0 0 5 5 20]; % MHz

Sys.lwpp = 0.1; % mT

% Experimental parameters
Exp.mwFreq = 9.75; % GHz
Exp.Range = [346 350]; % mT
Exp.Harmonic = 0;

% Thermalised triplet precursor

% Population vector in coupled basis (electron states only)
Sys.initState = {[1/3 1/3 1/3 0],'coupled'};
[x,spc1] = pepper(Sys,Exp);

% Population vector in coupled basis (electron and nuclear states)
Sys.initState = {[1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 0 0 0],'coupled'};
[~,spc2] = pepper(Sys,Exp);

% Density matrix in uncoupled basis (electron states only)
V = cgmatrix(Sys.S(1),Sys.S(2));
Tp = V(1,:)';
T0 = V(2,:)';
Tm = V(3,:)';
Sys.initState = 1/3*(Tp*Tp' + T0*T0' + Tm*Tm');
[~,spc3] = pepper(Sys,Exp);

% Plot
if opt.Display
  plot(x,spc1,x,spc2,x,spc3);
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: Spin-correlated radical pair coupled to 14N');
end

ok = areequal(spc1,spc2,1e-12,'abs') && areequal(spc1,spc3,1e-12,'abs');
