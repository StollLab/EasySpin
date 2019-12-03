function [err,data] = test(opt,olddata)

Sys1 = struct('g',[2 2.2],'lwpp',1);
Sys2 = struct('g',[2.05 2.1 2.15],'lwpp',2);
Exp = struct('mwFreq',9.7,'Range',[300 360]);
Opt = struct('Verbosity',opt.Verbosity);

Sys1.weight = 0.567;
Sys2.weight = 1.734;

[x,y1] = pepper(Sys1,Exp,Opt);
[x,y2] = pepper(Sys2,Exp,Opt);
[x,y] = pepper({Sys1,Sys2},Exp,Opt);

y_sum = y1 + y2;

err = ~areequal(y,y_sum,1e-10,'rel');

if (opt.Display)
  plot(x,y_sum,x,y);
end

data = [];
