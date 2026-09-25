function ok = test(opt)

% pepper with Method='eig' should include thermal populations
% and give the same integral intensity as Method='matrix'
% (S=1 powder at low temperature)

Sys.S = 1;
Sys.g = 2;
Sys.D = 3000;  % MHz
Sys.lwpp = 20;  % mT

Exp.mwFreq = 95;  % GHz
Exp.Range = [2000 4800];  % mT
Exp.Harmonic = 0;
Exp.Temperature = 2;  % K

Opt.Method = 'matrix';
[B,spc1] = pepper(Sys,Exp,Opt);
Opt.Method = 'eig';
[B,spc2] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(B,spc1,B,spc2);
  legend('matrix','eig');
end

dB = B(2)-B(1);
ok = areequal(sum(spc1)*dB,sum(spc2)*dB,1e-3,'rel');
