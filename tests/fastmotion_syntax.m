function ok = test()

%======================================================
% Syntax
%======================================================

System.g = [2.0088 2.0064 2.0027];
System.Nucs = '14N';
System.A = mt2mhz([7.59 5.95 31.76]/10); % MHz
tcorr = 1e-10; % s
field = 350; % mT

lw = fastmotion(System,field,tcorr);
[lw,mI] = fastmotion(System,field,tcorr);
[lw,mI,coeffs] = fastmotion(System,field,tcorr);

ok = true;
