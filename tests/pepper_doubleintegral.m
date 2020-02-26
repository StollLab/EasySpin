function ok = test()

% Test whether the double integral is independent of the
% number of nuclei

Exp.mwFreq = 9.5;
Exp.Harmonic = 1;
Exp.CenterSweep = [339 10];

Sys1.S = 1/2;
Sys1.lw = 0.2;
Sys2 = nucspinadd(Sys1,'14N',50);
Sys3 = nucspinadd(Sys2,'1H',20);

Opt.Method = 'matrix';
spec1 = pepper(Sys1,Exp,Opt);
spec2 = pepper(Sys2,Exp,Opt);
spec3 = pepper(Sys3,Exp,Opt);
Opt.Method = 'perturb';
spec4 = pepper(Sys1,Exp,Opt);
spec5 = pepper(Sys2,Exp,Opt);
spec6 = pepper(Sys3,Exp,Opt);

integrals1 = [trapz(cumtrapz(spec1)) trapz(cumtrapz(spec2)) trapz(cumtrapz(spec3))]/trapz(cumtrapz(spec1));
integrals2 = [trapz(cumtrapz(spec4)) trapz(cumtrapz(spec5)) trapz(cumtrapz(spec6))]/trapz(cumtrapz(spec4));

ok = areequal(integrals1,[1 1 1],1e-6,'abs') && areequal(integrals2,[1 1 1],1e-6,'abs');
