function ok = test(opt)

% Test whether the double integral is independent of the number of nuclei

Exp.mwFreq = 9.5;
Exp.Harmonic = 1;
Exp.CenterSweep = [339 10];

Sys1.S = 1/2;
Sys1.lw = 0.3;
Sys2 = nucspinadd(Sys1,'14N',50);
Sys3 = nucspinadd(Sys2,'1H',20);

Opt.Method = 'matrix';
[B,spec1] = pepper(Sys1,Exp,Opt);
spec2 = pepper(Sys2,Exp,Opt);
spec3 = pepper(Sys3,Exp,Opt);
Opt.Method = 'perturb';
spec4 = pepper(Sys1,Exp,Opt);
spec5 = pepper(Sys2,Exp,Opt);
spec6 = pepper(Sys3,Exp,Opt);

% Calculate double integrals
dint = @(f) trapz(cumtrapz(f));
integrals1 = [dint(spec1) dint(spec2) dint(spec3)]/dint(spec1);
integrals2 = [dint(spec4) dint(spec5) dint(spec6)]/dint(spec4);

ok(1) = areequal(integrals1,[1 1 1],1e-6,'abs');
ok(2) = areequal(integrals2,[1 1 1],1e-6,'abs');

if opt.Display
  subplot(2,1,1)
  plot(B,spec1,B,spec2,B,spec3)
  title('Opt.Method=''matrix''');
  legend('no nuclei','1H','14N')
  subplot(2,1,2)
  plot(B,spec4,B,spec5,B,spec6)
  title('Opt.Method=''perturb''');
  legend('no nuclei','1H','14N')
  integrals1
  integrals2
end
