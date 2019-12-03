function [err,data] = test(opt,olddata)

%===================================================================
% Check consistency of frequency-domain and field-domain linewidths
%===================================================================

g0 = 2;
System.g = g0;
System.Nucs = '15N';
System.A = [1 6]*10; % MHz
tcorr = 1e-10; % s
field = 350; % mT

[lw_mT,mI] = fastmotion(System,field,tcorr,'field');
[lw_MHz,mI] = fastmotion(System,field,tcorr,'freq');

err = ~areequal(lw_mT,mhz2mt(lw_MHz,g0),1e-10,'rel');

data = [];
