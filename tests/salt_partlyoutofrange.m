function [err,data] = test(opt,olddata)

% I=1/2 powder spectrum, transitions partly out of range

Sys = struct('S',1/2,'g',[2 2 2],'lwEndor',0.03);
Sys = nucspinadd(Sys,'1H',[-1 -1 2]*2+3);
Exp = struct('mwFreq',9.5,'Field',339.377,'Range',[13 19]);
Opt = struct('Verbosity',opt.Verbosity);
[x,y] = salt(Sys,Exp,Opt);
if (opt.Display)
  plot(x,y);
  title('I=1/2 powder spectrum ,transitions partly out of range');
end

err = 0;
data = [];
