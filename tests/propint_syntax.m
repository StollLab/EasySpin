function [err,data] = test(opt,olddata)

% Test 1: syntax
%======================================================
H0 = complex(rand(4),rand(4)); H0 = H0+H0';
H1 = complex(rand(4),rand(4)); H1 = H1+H1';
t = [0 10];
freq = 0.5768;
phase = 0.1;
n = 100;
U = propint(H0,H1,t,freq);
U = propint(H0,H1,t,freq,phase);
U = propint(H0,H1,t,freq,phase,n);
[U,S] = propint(H0,H1,t,freq);
U = propint(S,t,freq);
err = 0;
data = [];
