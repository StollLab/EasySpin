function [err,data] = test(opt,olddata)

% Direct value check

a = gfree;
%b = 2.00231930437180;  % 2002 CODATA
%b = 1.00115965218085*2;  % 2006 Gabrielse
b = 1.00115965218073*2;  % 2008 Gabrielse
err = abs(a-b)/b > 1e-12;
data = [];
