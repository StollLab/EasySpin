function [err,data] = test(opt,olddata)

%======================================================
% Linewdiths of a nitroxide radical
%======================================================

System.g = 2;
System.Nucs = '14N';
System.A = 10;
System.Q = [-1 -1 2]*5;
tcorr = 1e-10;
field = 350;

[lw,mI] = fastmotion(System,field,tcorr);

lw0 = [0.001010068477786; 0.000673378985191; 0.001010068477786];
err = ~areequal(lw,lw0,1e-8,'abs');
data = [];
