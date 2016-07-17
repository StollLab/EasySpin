function [err,data] = test(opt,olddata)

%======================================================
% Comparison of wignerd output with explicit expression
%======================================================

beta = rand*pi;

d1 = wignerd(1,beta,'-');

% Varshalovich, Quantum Theory of Angular Momentum, p.199, Table 4.4
% Brink/Satcher, Angular Momentum, p.24, Table I
% Richard N. Zare, Angular Momentum, p.89, Table 3.1
d2(1,1) = cos(beta/2)^2;
d2(3,3) = d2(1,1);
d2(1,3) = sin(beta/2)^2;
d2(3,1) = d2(1,3);
d2(2,1) = sqrt(1/2)*sin(beta);
d2(3,2) = d2(2,1);
d2(2,3) = -d2(2,1);
d2(1,2) = -d2(2,1);
d2(2,2) = cos(beta);

err = ~areequal(d1,d2,1e-10);
data = [];