function [err,data] = test(opt,olddata)

%===============================================================================
% Compare Freed and general code with various potentials
%===============================================================================

Sys.g = [2.05 2.02 2.00];
Sys.logDiff = 7;

Exp.Field = 340;
Exp.mwRange = [9.4 9.9];
Exp.nPoints = 512;

Opt.LLKM = [10 5 2 2];

for p = 1:4
  Sys.lambda = [0 0 0 0];
  Sys.lambda(p) = +2;
  
  Opt.LiouvMethod = 'Freed';
  [B,y1] = chili(Sys,Exp,Opt);
  Opt.LiouvMethod = 'general';
  [B,y2] = chili(Sys,Exp,Opt);
  
  err(p) = ~areequal(y1,y2,1e-6*max(abs(y1)));
  
  if opt.Display
    subplot(4,1,p);
    plot(B,y1,B,y2);
    legend('Freed','general');
    title(sprintf('lambda(%d) = %g',p,Sys.lambda(p)));
  end
  
end

err = any(err);
data = [];
