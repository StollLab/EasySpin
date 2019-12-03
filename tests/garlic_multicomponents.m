function [err,data] = test(opt,olddata)

Sys1 = struct('g',2,'Nucs','1H','lw',[0,0.1],'A',20);
Sys2 = struct('g',2,'Nucs','14N','lw',[0,0.03],'A',34);
Sys3 = struct('g',2,'Nucs','C','lw',[0,0.03],'A',5);
Exp = struct('mwFreq',9.7,'Range',[344 349]);
Opt = struct('Verbosity',opt.Verbosity);

Sys1.weight = 0.567;
Sys2.weight = 1.567;

[x,y1] = garlic(Sys1,Exp,Opt);
[x,y2] = garlic(Sys2,Exp,Opt);
[x,y3] = garlic(Sys3,Exp,Opt);
[x,y] = garlic({Sys1,Sys2,Sys3},Exp,Opt);

y_sum = y1 + y2 + y3;

err = ~areequal(y,y_sum,1e-8,'rel');

if (opt.Display)
  plot(x,y_sum,x,y,x,y_sum-y,'r');
end

data = [];
