function ok = test(opt)

%-------------------------------------------------------------
%S=3/2, I=1/2
%-------------------------------------------------------------

Sys.S = 3/2;
Sys.D = 200;
Sys.Nucs='1H';
Sys.A = [3 5];
Sys.lwEndor = 0.2;
Exp.Field = 350;
Exp.mwFreq = 9.5;
Exp.Range = [0 30];
Opt.Method='matrix';
[x,a]=salt(Sys,Exp,Opt);
Opt.Method='perturb2';
[x,b]=salt(Sys,Exp,Opt);
b=rescaledata(b,a,'minmax');
if opt.Display
  plot(x,a,'k',x,b,'r');
  legend('matrix','perturb2');
end

ok = areequal(b,a,0.01,'rel');
