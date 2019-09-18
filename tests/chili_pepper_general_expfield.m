function [err,data] = test(opt,olddata)

%=============================================================================
% Intensity match for rigid limit chili and pepper explicit B sweep (general)
%=============================================================================

Sys.g = [2.01 2.003];
Sys.tcorr = 1e-5;
Sys.lw = 0.2;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.Range = [337 339.5];


[x,y1] = pepper(Sys,Exp);

Opt.LiouvMethod = 'general';
Opt.ExplicitFieldSweep = true;
Opt.LLMK = [20 0 0 0];

[x2,y2] = chili(Sys,Exp,Opt);

if opt.Display
  subplot(2,1,1);
  plot(x,y1,x2,y2);
  legend('pepper','chili');
  title('unscaled');
  subplot(2,1,2);
  plot(x,y1/max(y1),x2,y2/max(y2));
  legend('pepper','chili');
  title('scaled');
end

err = ~areequal(y1/max(y1),y2/max(y2),0.02);

data = [];
