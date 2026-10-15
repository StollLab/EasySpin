function ok = test(opt)

% Test whether fast-motion spectrum is invariant with respect to the sign
% of the hyperfine tensor (linewidths must be assigned to the correct lines).

Sys.g = [2.008 2.006 2.003];
Sys.Nucs = '14N';
Sys.tcorr = 3e-10;

Exp.mwFreq = 9.5;
Exp.Range = [334 344];

Sys.A = +[20 20 100]; [B,spc1] = garlic(Sys,Exp);
Sys.A = -[20 20 100]; [B,spc2] = garlic(Sys,Exp);

ok = areequal(spc1,spc2,1e-6,'rel');

if opt.Display
  plot(B,spc1,B,spc2,'--');
  legend('positive A','negative A');
  legend boxoff
end
