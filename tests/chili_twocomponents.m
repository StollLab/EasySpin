function [err,data] = test(opt,olddata)

Sys = struct('g',[2.008 2.0061 2.0027],'Nucs','14N','A',[16 16 86]);
Sys.lw = 0.1;

Sys1 = Sys;
Sys2 = Sys;
Sys1.tcorr = 0.2e-9;
Sys2.tcorr = 3e-9;

Exp = struct('mwFreq',9.7,'Range',[340 354]);
Opt = struct('Verbosity',opt.Verbosity);

Sys1.weight = 0.05;
Sys2.weight = 1;

[x,y1] = chili(Sys1,Exp,Opt);
[x,y2] = chili(Sys2,Exp,Opt);
[x,y] = chili({Sys1,Sys2},Exp,Opt);

y_sum = y1 + y2;

err = ~areequal(y,y_sum,1e-10,'rel');

if (opt.Display)
  plot(x,y_sum,x,y);
  legend('separate/added','multi');
end

data = [];
