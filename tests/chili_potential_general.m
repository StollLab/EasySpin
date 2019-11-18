function [err,data] = test(opt,olddata)

%===============================================================================
% Test different types of potentials (M=K=0, M=0, K=0, both nonzero)
%===============================================================================

Sys.g = [2.05 2.00];
Sys.logDiff = 7;

Exp.Field = 340;
Exp.mwRange = [9.4 9.9];
Exp.nPoints = 200;
Exp.CrystalOrientation = [0 0 0];

Opt.LLMK = [4 0 2 2];
Opt.highField = true;

lam = +2;
LMKlam = [0 0 0 0; 2 0 0 lam; 2 1 0 lam; 2 0 1 lam; 2 1 1 lam];

for p = 1:5
  Sys.Potential = LMKlam(p,:);
  
  Opt.LiouvMethod = 'general';
  [nu,spc{p}] = chili(Sys,Exp,Opt);
end

if opt.Display
  for p = 1:5
    subplot(5,1,p);
    plot(nu,spc{p});
    axis tight
    title(sprintf('lambda(%d,%d,%d) = %g',...
      LMKlam(p,1),LMKlam(p,2),LMKlam(p,3),LMKlam(p,4)));
  end
end

err = false;
data = [];
