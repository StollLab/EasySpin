function [ok,data] = test(opt,olddata)

%=======================================================================
% Partially oriented system with axial g tensor and axial ordering
%=======================================================================
Sys.g = [2, 2, 2.2];
Sys.lw = 1;
Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.Range = [300 350];
[B,spc] = pepper(Sys,Exp);

lambda = [-5 -2 2 5];

for k = 1:numel(lambda)
  Exp.Ordering = lambda(k);
  [B,spc(k+1,:)] = pepper(Sys,Exp);
end

if opt.Display
  if ~isempty(olddata)
    subplot(3,1,1);
    plot(B,spc); title('current');
    subplot(3,1,2);
    plot(B,olddata.spc); title('old');
    subplot(3,1,3);
    plot(B,olddata.spc-spc); title('difference');
  else
    plot(B,spc);
    xlabel('magnetic field (mT)');
    ylabel('intensity (arb.u.)');
    title('pepper: axial orientational distribution');
  end
end

data.spc = spc;

if ~isempty(olddata)
  ok = areequal(spc/max(spc),olddata.spc/max(olddata.spc),1e-4,'abs');
else
  ok = [];
end
