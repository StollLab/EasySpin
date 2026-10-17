function ok = test(opt)

% Sys.initState with the hybrid method (perturbational nuclei, electron-only
% core system) agrees with the matrix method, for field and frequency sweeps

Sys = struct('S',1,'D',[900 160],'Nucs','1H','A',[2 3 4]);
Sys.initState = {[0.1 0.5 0.4],'zerofield'};

% Field sweep (resfields)
Sys.lwpp = 0.5;  % mT
Exp = struct('mwFreq',9.5,'Range',[250 400],'Harmonic',0,'SampleFrame',[0.2 0.7 0.1]);
Opt.Method = 'matrix';
[B,spcM] = pepper(Sys,Exp,Opt);
Opt.Method = 'hybrid';
spcH = pepper(Sys,Exp,Opt);
spcH = spcH/max(abs(spcM));
spcM = spcM/max(abs(spcM));
ok(1) = areequal(spcH,spcM,1e-2,'abs');

% Frequency sweep (resfreqs_matrix)
Sys.lwpp = 14;  % MHz
Expf = struct('Field',340,'mwRange',[8 11],'Harmonic',0,'SampleFrame',[0.2 0.7 0.1]);
Opt.Method = 'matrix';
[nu,spcfM] = pepper(Sys,Expf,Opt);
Opt.Method = 'hybrid';
spcfH = pepper(Sys,Expf,Opt);
spcfH = spcfH/max(abs(spcfM));
spcfM = spcfM/max(abs(spcfM));
ok(2) = areequal(spcfH,spcfM,1e-2,'abs');

if opt.Display
  subplot(2,1,1)
  plot(B,spcM,B,spcH);
  legend('matrix','hybrid');
  xlabel('magnetic field (mT)');
  subplot(2,1,2)
  plot(nu,spcfM,nu,spcfH);
  legend('matrix','hybrid');
  xlabel('frequency (GHz)');
end
