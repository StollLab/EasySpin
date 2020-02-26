function ok = test(opt,olddata)

eps0_val = eps0;
eps0_ref = 8.8541878128e-12; % = 1/mu0/clight^2

ok = areequal(eps0_val,eps0_ref,1e-10,'rel');
