function ok = test(opt)

% Comparison fast and general code for Liouvillian

Sys.g = [2.008 2.0061 2.0027];
Sys.Nucs = '14N';
Sys.A = [20 20 100];
Sys.tcorr = 1e-9;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;

Opt.LiouvMethod = 'fast';
[x,y1] = chili(Sys,Exp,Opt);
Opt.LiouvMethod = 'general';
[x,y2] = chili(Sys,Exp,Opt);

if opt.Display
  plot(x,y1,x,y2);
  legend('fast','general');
end

ok = areequal(y1,y2,1e-3,'rel');
