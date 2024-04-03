function ok = test(opt)

%=======================================================================
% "Partially ordered" system with isotropic orientational distribution
%=======================================================================
Sys.g = [2, 2.05, 2.2];
Sys.lw = 1;
Exp.mwFreq = 9.5;
Exp.Range = [300 350];
Exp.Harmonic = 0;

Opt.GridSymmetry = 'Ci';

% Isotropic distribution
Exp.Ordering = [];
[B,spc0] = pepper(Sys,Exp,Opt);

% Explicitly specified isotropic distribution
Exp.Ordering = @(a,b,c) ones(size(b));
[B,spc1] = pepper(Sys,Exp,Opt);

if opt.Display
  subplot(2,1,1)
  plot(B,spc0,B,spc1);
  legend('no ordering function','isotropic ordering function');
  legend boxoff
  subplot(2,1,2)
  plot(B,spc0-spc1);
  legend('difference');
end

ok = areequal(spc0,spc1,1e-10,'rel');

end
