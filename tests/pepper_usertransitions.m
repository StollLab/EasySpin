function ok = test(opt)

Sys = struct('S',1,'g',[2 2 2],'D',[-1 -1 2]*500,'lw',10);
Exp = struct('Range',[0 1000],'mwFreq',9.5,'nPoints',2e3);
Opt = struct('Verbosity',0,'Output','transitions');
Opt.Transitions = [1 2; 1 3; 2 3];

[x,y] = pepper(Sys,Exp,Opt);

ok = size(y,1)==size(Opt.Transitions,1);

if opt.Display
  title('User-specified transitions, separate output');
  plot(x,y);
end
