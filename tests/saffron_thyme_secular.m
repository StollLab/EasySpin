function ok = test()

% In a shifted simulation frame, saffron uses the secular Hamiltonian:
% block-diagonal in the total electron spin projection, with the exact
% lab-frame energies. Results then do not depend on the simulation frame
% or on the rotation of the molecule about the field.

% Secular Hamiltonian: commutes with electron Sz, exact lab-frame energies
Sys.S = 1;
Sys.g = [2.00 2.01 2.02];
Sys.gFrame = [0.2 0.7 0.4];
Sys.D = [600 100];
Sys.DFrame = [0.5 1.0 0.2];
Sys.Nucs = '14N';
Sys.A = [5 20];
Sys.AFrame = [1.1 0.3 0];
Sys = runprivate('validatespinsys',Sys);
Sz = sop(Sys.Spins,'ze');
Hlab = ham(Sys,[0 0 340]);
Hlab = (Hlab+Hlab')/2;
Hsec = runprivate('s_secularham',Hlab,Sys.Spins,Sys.nElectrons);
ok(1) = norm(Hsec*Sz-Sz*Hsec)<1e-9*norm(Hsec);
ok(2) = areequal(eig(Hsec),eig(Hlab),1e-9,'rel');

% saffron: independence of simulation frame and of chi
Pulse.tp = 0.02;
Pulse.Flip = pi/4;
Exp.Sequence = {Pulse 0.05};
Exp.Field = 340;
Exp.mwFreq = 9.53;
Exp.DetWindow = [0 0.05];
Opt.SimulationMode = 'thyme';

SysT.S = 1;
SysT.g = 2.0028;
SysT.D = [60 12];

tol = 5e-3;
tc = linspace(0,0.049,50);  % inside both time axes (interp1 gives NaN outside)

OptP = Opt;
OptP.GridSymmetry = 'O3';
OptP.GridFrame = [0.4 0.9 0];
SimFreq = [5 7];
for k = 1:numel(SimFreq)
  OptP.SimFreq = SimFreq(k);
  [t,y] = saffron(SysT,Exp,OptP);
  yF(k,:) = interp1(t,y,tc); %#ok
end
ok(3) = areequal(yF(1,:),yF(2,:),tol,'rel');

ExpC = Exp;
ExpC.SampleFrame = [0 0 0];
OptC = Opt;
OptC.SimFreq = 7;
chi = [0 2];
for k = 1:numel(chi)
  ExpC.MolFrame = -fliplr([0.4 0.9 chi(k)]);
  [~,yC(k,:)] = saffron(SysT,ExpC,OptC); %#ok
end
ok(4) = areequal(yC(1,:),yC(2,:),tol,'rel');

% Manifold weights: S=1/2 in a field tilted by theta from z has weight
% cos(theta/2)^2 in each mS manifold. Below 0.7, s_secularham errors.
H = @(theta) 1000*(cos(theta)*sop(1/2,'z') + sin(theta)*sop(1/2,'x'));
[~,minWeight] = runprivate('s_secularham',H(pi/3),1/2,1);
ok(5) = areequal(minWeight,cos(pi/6)^2*[1;1],1e-12,'abs');
try
  runprivate('s_secularham',H(84*pi/180),1/2,1);
  ok(6) = false;
catch err
  ok(6) = contains(err.message,'too strongly mixed');
end
