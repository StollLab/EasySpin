function [err,data] = test(opt,olddata)


Sys1.Nucs = '1H';
Sys1.A = [2 2 5];
Sys1.lwEndor = 0.3;
Sys2.Nucs = '1H';
Sys2.A = [8 10 12];
Sys2.lwEndor = 0.2;

Sys1.weight = 0.567;
Sys2.weight = 0.734;

Exp.Field = 330;
Exp.Range = [0 30];

Opt = struct;

[x,y1] = salt(Sys1,Exp,Opt);
[x,y2] = salt(Sys2,Exp,Opt);
[x,y] = salt({Sys1,Sys2},Exp,Opt);

y_sum = y1 + y2;

err = ~areequal(y,y_sum);

if (opt.Display)
  plot(x,y_sum,x,y);
end

data = [];
