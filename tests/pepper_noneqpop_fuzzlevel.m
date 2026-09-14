function ok = test(opt)

% Check FuzzLevel behavior for spin-polarized spin systems

% Spin system
Sys.S = [1/2 1/2];
Sys.g = [2.0027; 2.0000];
Sys.J = -6; % MHz

Sys.Nucs = '1H,1H';
Sys.A = [0 0 0 5 5 20; 0 0 0 2 2 10]; % MHz

Sys.lwpp = 0.1; % mT

Sys.initState = 'singlet';

% Experimental parameters
Exp.mwFreq = 9.75; % GHz
Exp.Range = [346 350]; % mT
Exp.Harmonic = 0;

% Default options
[x,spcDefault] = pepper(Sys,Exp);

% Explicitly set Opt.FuzzLevel = 0
Opt.FuzzLevel = 0;
[x,spcZeroFuzz] = pepper(Sys,Exp,Opt);

% Check that non-zero FuzzLevel is applied
Opt.FuzzLevel = 1e-10;
[x,spcFuzz] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(x,spcDefault,x,spcZeroFuzz,x,spcFuzz);
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: FuzzLevel for spin-polarized spin systems');
  legend('default','FuzzLevel = 0','FuzzLevel = 1e-10');
end

ok(1) = areequal(spcDefault,spcZeroFuzz,1e-12,'abs');
ok(2) = ~areequal(spcZeroFuzz,spcFuzz,1e-12,'abs');