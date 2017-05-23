clear

%==============================================================
% Compare intensities for rigid limit chili and pepper (Freed)
%==============================================================

Sys.g = [2.01 2.003];
Sys.tcorr = 1000e-9;
Sys.lw = 5;

Exp.Field = 339;
Exp.Harmonic = 0;
Exp.mwRange = [9.49 9.555];

[x,y1] = pepper(Sys,Exp);

Opt.LiouvMethod = 'Freed';
%Opt.Solver = '\';
%Opt.ExplicitFieldSweep = true;
Opt.LLKM = [20 0 0 0];

[x,y2] = chili(Sys,Exp,Opt);

if 1
  plot(x,y1,x,y2);
  legend('pepper','chili');
end

%err = ~areequal(y1,y2,9e-1*max(y1));

data = [];

%%

clear
%==============================================================
% Compare intensities for rigid limit chili and pepper (Freed)
%==============================================================

Sys.g = [2.01 2.003];
Sys.tcorr = 1000e-9;
Sys.lw = 0.2;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.Range = [337 339.5];

[x,y1] = pepper(Sys,Exp);

Opt.LiouvMethod = 'Freed';
Opt.Solver = '\';
%Opt.ExplicitFieldSweep = true;
Opt.LLKM = [20 0 0 0];

[x,y2] = chili(Sys,Exp,Opt);

if 1
  plot(x,y1,x,y2);
  legend('pepper','chili');
end

%err = ~areequal(y1,y2,9e-1*max(y1));

data = [];
