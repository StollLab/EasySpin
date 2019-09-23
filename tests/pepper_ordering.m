function [err,data] = test(opt,olddata)

%=======================================================================
% Partially oriented system with axial spin parameters
%=======================================================================
Sys = struct('g',[2, 2, 2.2],'lw',1);
Exp = struct('mwFreq',9.5,'Range',[300 350]);
Opt = struct('Verbosity',0);
[x,spc] = pepper(Sys,Exp,Opt);

lambda = -10:2:10;

for k = 1:numel(lambda)
  Exp.Ordering = lambda(k);
  [x,spc(k+1,:)] = pepper(Sys,Exp,Opt);
end

if (opt.Display)
  if ~isempty(olddata)
    subplot(3,1,1);
    plot(x,spc); title('current');
    subplot(3,1,2);
    plot(x,olddata.spc); title('old');
    subplot(3,1,3);
    plot(x,olddata.spc-spc); title('difference');
  else
    plot(x,spc);
    xlabel('magnetic field [mT]');
    ylabel('intensity [a.u.]');
    title('pepper: non-isotropic orientational distribution, axial system');
  end
end

data.spc = spc;

if ~isempty(olddata)
  err = ~areequal(spc/max(spc),olddata.spc/max(olddata.spc),1e-4,'abs');
else
  err = [];
end
