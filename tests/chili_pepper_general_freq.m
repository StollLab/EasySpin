function [err,data] = test(opt,olddata)

%============================================================================
% Intensity match for rigid limit chili and pepper frequency sweep (general)
%============================================================================

Sys.g = [2.01 2.003];
Sys.tcorr = 1e-5;
Sys.lw = 5;

Exp.Field = 339;
Exp.Harmonic = 0;
Exp.mwRange = [9.49 9.555];

[x,y1] = pepper(Sys,Exp);

Opt.LiouvMethod = 'general';
Opt.LLMK = [20 0 0 0];

[x,y2] = chili(Sys,Exp,Opt);

if opt.Display
  subplot(2,1,1)
  plot(x,y1,x,y2);
  title('unscaled')
  legend('pepper','chili');
  subplot(2,1,2)
  plot(x,y1/max(y1),x,y2/max(y2));
  title('scaled')
  legend('pepper','chili');
end

err = ~areequal(y1/max(y1),y2/max(y2),0.01);

data = [];