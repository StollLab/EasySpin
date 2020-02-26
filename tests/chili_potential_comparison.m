function ok = test(opt)

% Compare fast and general code with various potentials

Sys.g = [2.05 2.00];
Sys.logDiff = 7;

Exp.Field = 340;
Exp.mwRange = [9.4 9.9];
Exp.nPoints = 512;
Exp.CrystalOrientation = [0 0 0];

Opt.LLMK = [10 0 2 4];

lam = +2;
LMK = [2 0 0; 2 0 2; 4 0 0; 4 0 2];

for p = 1:4
  Sys.Potential = [LMK(p,:) lam];
  
  Opt.LiouvMethod = 'fast';
  [~,y_fast] = chili(Sys,Exp,Opt);
  Opt.LiouvMethod = 'general';
  [nu,y_general] = chili(Sys,Exp,Opt);
  
  ok(p) = areequal(y_fast,y_general,1e-4,'rel');
  
  if opt.Display
    subplot(4,1,p);
    plot(nu,y_fast,nu,y_general);
    axis tight
    legend('fast','general');
    title(sprintf('lambda(%d,%d,%d) = %g',LMK(p,1),LMK(p,2),LMK(p,3),lam));
  end
  
end
