function ok = test(opt)

% Test whether double integral is independent of the number of hyperfine lines

Sys.g = [2 2.001];

Sys.A = [50 100];
Sys.tcorr = 1e-6;
Sys.lw = 5;

Exp.Field = 350;
Exp.mwRange = [9.4 10.2];
Exp.nPoints = 3e3;

Opt.LLMK = [30 0 0 0];

A = [200 400]*2;
Sys.Nucs = '15N';
Sys.A = A/2;
[B,spc2] = chili(Sys,Exp,Opt);
Sys.Nucs = '14N';
Sys.A = A/3;
[B,spc3] = chili(Sys,Exp,Opt);
Sys.Nucs = '63Cu';
Sys.A = A/4;
[B,spc4] = chili(Sys,Exp,Opt);

dB = B(2)-B(1);
int2 = sum(spc2)*dB;
int3 = sum(spc3)*dB;
int4 = sum(spc4)*dB;

ok = areequal(int2,int3,0.01,'rel') && areequal(int3,int4,0.01,'rel');

if opt.Display
  plot(B,spc2,B,spc3,B,spc4);
  legend('I=1/2','I=1','I=3/2');
end
