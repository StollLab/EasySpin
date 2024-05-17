function ok = test()

t_us = 1.1;  % µs
conc_uM = 100;  % µM
lambda = 0.7;

c_mM = conc_uM/1e3;  % µM -> mM=mol/m^3
c = c_mM*avogadro;  % mol/m^3 -> m^-3

g = gfree;
D = mu0/(4*pi)*g^2*bmagn^2/hbar;  % rad/s
k = 8*pi^2/(9*sqrt(3))*D*lambda*c;
Vinter_ref = exp(-k*t_us*1e-6);
Vinter = dipbackground(t_us,conc_uM,lambda);

ok = areequal(Vinter,Vinter_ref,1e-10,'rel');

end
