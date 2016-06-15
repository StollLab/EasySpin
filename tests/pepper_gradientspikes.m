function [err,data] = test(opt,olddata)

% Check for gradient smoothing spikes

Sys = struct('S',1/2,'g',[2 2 3],'HStrain',[1 1 2]*0.001);
Exp = struct('mwFreq',9.5,'Range',[200 400],'Harmonic',0);
Opt = struct('Verbosity',opt.Verbosity,'nKnots',[10 0]);

Opt.nKnots = 10;
[x,y1] = pepper(Sys,Exp,Opt);
Opt.nKnots = 61;
[x,y2] = pepper(Sys,Exp,Opt);

if (opt.Display)
  plot(x,y1,'r',x,y2,'b');
  title('Absence of gradient smoothing spikes');
end

err = max(y1)>max(y2);
data = [];
