function [err,data] = test(opt,olddata)

%======================================================
% Comparison of wignerd output with explicit expression
%======================================================

beta = rand*pi;

d1 = wignerd(1/2,beta,'-');

% Varshalovich, Quantum Theory of Angular Momentum, p.199, Table 4.3
% Brink/Satcher, Angular Momentum, p.24, Table I
% Richard N. Zare, Angular Momentum, p.89, Table 3.1
d2(1,1) = cos(beta/2);
d2(2,2) = d2(1,1);
d2(1,2) = -sin(beta/2);
d2(2,1) = -d2(1,2);

err = ~areequal(d1,d2,1e-10,'abs');
data = [];