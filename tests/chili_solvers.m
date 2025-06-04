function ok = test(opt)

% Compare all solvers

Sys.g = [2 2 2.1];
Sys.tcorr = 1e-9;

Exp.Field = 350;
Exp.mwRange = [9 11];
Exp.nPoints = 100;

Opt.LLMK = [8 0 2 0];

Solvers = {'L','\','E','C','B'};
for s = 1:numel(Solvers)
  Opt.Solver = Solvers{s};
  [nu,spc(s,:)] = chili(Sys,Exp,Opt);
end

% Assort that all spectra are close
spread = max(max(spc)-min(spc))/max(spc(:));

ok = spread<5e-3;

if opt.Display
  subplot(3,1,[1 2]);
  plot(nu,spc)
  legend(Solvers{:},'Interpreter','none');
  legend boxoff
  subplot(3,1,3);
  plot(nu,spc-spc(1,:));
end

end
