function ok = test(opt)

% For a multi-component or multi-isotopologue powder simulation, return
% transitions for all components separately if Opt.Output='separate'

% (1) multiple isotopologues
%----------------------------------------------------
Sys.lw = 0.1;
Sys.Nucs = 'Cu';
Sys.A = 50;

Exp.mwFreq = 9.3260;
Exp.Range = [325 340];

Opt.Output = 'separate';
[B,spc1] = pepper(Sys,Exp,Opt);

Opt.Output = 'summed';
[B,spc2] = pepper(Sys,Exp,Opt);

ok1 = size(spc1,1)==2;
ok1 = ok1 && max(abs(sum(spc1,1)-spc2))/max(abs(spc2))<1e-10;

if opt.Display
  subplot(2,1,1);
  plot(B,spc1,B,spc2)
end

% (2) multiple components
%----------------------------------------------------
Sys1 = Sys;
Sys2 = Sys;
Sys1.Nucs = '63Cu';
Sys2.Nucs = '65Cu';
Sys2.A = Sys1.A*1.4;
Sys2.weight = 2;

Opt.Output = 'separate';
[B,spc1] = pepper({Sys1,Sys2},Exp,Opt);

Opt.Output = 'summed';
[B,spc2] = pepper({Sys1,Sys2},Exp,Opt);

ok2 = size(spc1,1)==2;
ok2 = ok2 && max(abs(sum(spc1,1)-spc2))/max(abs(spc2))<1e-10;

if opt.Display
  subplot(2,1,2);
  plot(B,spc1,B,spc2);
end

ok = ok1 && ok2;
