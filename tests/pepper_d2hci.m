function ok = test(opt)

% Simulations using D2h vs Ci vs C2h

% Define spin system with C2h symmetry
Sys.g = [1.9,2,2.3];
Sys.Nucs = '1H';
Sys.A = [20 80 300];
Sys.AFrame = [-pi/4 -pi/2 0];

SysPg = symm(Sys);
if ~strcmp(SysPg,'C2h')
  error('Incorrect symmetry %s of test system, should be C2h.',SysPg);
end

Exp = struct('Range',[285 365],'mwFreq',9.5,'nPoints',2e3,'Harmonic',0);
Opt = struct('nKnots',[20 4],'Method','perturb');
GridSymmetry = {'Ci','C2h','D2h'};

for k = 1:length(GridSymmetry)
  Opt.Symmetry = GridSymmetry{k};
  [x,y(k,:)] = pepper(Sys,Exp,Opt);
end

ok = areequal(y(1,:),y(2,:),5e-3,'rel');

if opt.Display
  title('D2h, C2h, Ci computation (D2h is incorrect)');
  plot(x,y);
  legend(GridSymmetry{:});
end
