function ok = test()
%orbital angular momenta can also be defined as spins, therefore the two
%Hamiltonians should be identical

rng(5);

% Build spin system with spin in Sys.S and orbital angular momentum in Sys.L
SysSL.S = randi(3)/2;
SysSL.g = rand(3,3);
SysSL.L = randi(2);
SysSL.soc = rand(1,2)*1000;
SysSL.gL = rand;

% Build spin system with both spin and orbital angular momentum in Sys.S
SysS.S = [SysSL.S SysSL.L];
SysS.g = [SysSL.g; SysSL.gL*eye(3)];
SysS.ee = SysSL.soc(1);  % spin-orbit coupling as bilinear spin-spin coupling
SysS.ee2 = SysSL.soc(2);  % and as biquadratic spin-spin coupling

% Build zero-field splitting and crystal-field part
for k = 2:2:8
  lfieldname = sprintf('CF%d',k);
  sfieldname = sprintf('B%d',k);
  SysSL.(sfieldname) = rand(1,2*k+1)*(k/2<=SysSL.S);
  SysSL.(lfieldname) = rand(1,2*k+1)*(k/2<=SysSL.L);
  SysS.(sfieldname) = [SysSL.(sfieldname); SysSL.(lfieldname)];
end

% Build random experiment for frequency-domain pepper, using a few random
% crystal orientations (much faster than a powder)
FDExp.Temperature = rand * 300;
FDExp.Field = rand *1e3;
FDExp.MolFrame = [0 0 0];
FDExp.SampleFrame = rand(3,3)*pi;

% Switch off random Hamiltonian fuzzing, so that both spin systems give
% identical spectra
Opt.FuzzLevel = 0;

% Compare S&L with S-only spin system
[nu,fd1] = pepper(SysSL,FDExp,Opt);
fd2 = pepper(SysS,FDExp,Opt);
ok(1) = areequal(fd1,fd2,1e-12,'rel');

% Build field-sweep experimet based on FD sim, always a transition in spectral window 
[~, ind] = max(fd1);
Exp.mwFreq = nu(ind);
Exp.CenterSweep = FDExp.Field*[1 0.5];
Exp.Temperature = FDExp.Temperature;
Exp.MolFrame = FDExp.MolFrame;
Exp.SampleFrame = FDExp.SampleFrame;

s1 = pepper(SysSL,Exp,Opt);
s2 = pepper(SysS,Exp,Opt);
ok(2) = areequal(s1,s2,1e-12,'rel');

% Test with an added nucleus
%-------------------------------------------------------------------------------
SysSL.Nucs = '1H';
SysSL.A = rand(3,3);
SysS.Nucs = SysSL.Nucs;
SysS.A = [SysSL.A, zeros(3,3)];
Opt.Method = 'hybrid';

% Frequency sweep
fd3 = pepper(SysSL,FDExp,Opt);
fd4 = pepper(SysS,FDExp,Opt);
ok(3) = areequal(fd3,fd4,1e-12,'rel');

% Field sweep
s3 = pepper(SysSL,Exp,Opt);
s4 = pepper(SysS,Exp,Opt);
ok(4) = areequal(s3,s4,1e-12,'rel');
