function err = test(opt)

%=======================================================
% Compare all solvers
%=======================================================
Sys.g = [2 2 2.1];
Sys.tcorr = 1e-9;

Exp.Field = 350;
Exp.mwRange = [9 11];
Exp.nPoints = 100;

Opt.LLMK = [14 0 2 0];

Solvers = {'L','\','E','C','B'};
for s = 1:numel(Solvers)
  Opt.Solver = Solvers{s};
  [nu,spc(s,:)] = chili(Sys,Exp,Opt);
end

if opt.Display
  plot(nu,spc)
  legend(Solvers{:},'Interpreter','none');
  legend boxoff
end

% Assort that all spectra are close
spread = max(max(spc)-min(spc))/max(spc(:));

err = spread>1e-3;
