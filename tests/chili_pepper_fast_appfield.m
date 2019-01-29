function [err,data] = test(opt,olddata)

%==============================================================================
% Intensity match for rigid limit chili and pepper approximate B sweep (fast)
%==============================================================================

Sys.g = [2.01 2.003];
Sys.tcorr = 1000e-9;
Sys.lw = 0.2;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.Range = [337 339.5];


[x,y1] = pepper(Sys,Exp);

Opt.LiouvMethod = 'fast';
Opt.LLMK = [20 0 0 0];

[x2,y2] = chili(Sys,Exp,Opt);

if opt.Display
  plot(x,y1,x2,y2);
  legend('pepper','chili');
end

err = ~areequal(y1,y2,1e-1*max(y1)); % 10% tolerance limit

data = [];
