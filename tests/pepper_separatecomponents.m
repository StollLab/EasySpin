function [err,data] = test(opt,olddata)

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

err1 = size(spc1,1)~=8;
err1 = err1 || max(abs(sum(spc1,1)-spc2))/max(abs(spc2))>1e-10;

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

err2 = size(spc1,1)~=8;
err2 = err2 || max(abs(sum(spc1,1)-spc2))/max(abs(spc2))>1e-10;

if opt.Display
  subplot(2,1,2);
  plot(B,spc1,B,spc2);
end

err = err1 || err2;
data = [];
