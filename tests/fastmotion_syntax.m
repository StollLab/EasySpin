function [err,data] = test(opt,olddata)

%======================================================
% Syntax
%======================================================

System = struct('g',[2.0088 2.0064 2.0027],'Nucs','14N');
System.A = mt2mhz([7.59 5.95 31.76]/10);
tcorr = 1e-10;
field = 350;
lw = fastmotion(System,field,tcorr);
[lw,mI] = fastmotion(System,field,tcorr);
[lw,mI,coeffs] = fastmotion(System,field,tcorr);

err = 0;
data = [];
