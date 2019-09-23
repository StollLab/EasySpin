function [err,data] = test(opt,olddata)

%======================================================
% Comparison of wignerd output with explicit expression
%======================================================

beta = rand*pi;

d1 = wignerd(3/2,beta,'-');

% '-' definition:
% Varshalovich, Quantum Theory of Angular Momentum, p.199, Table 4.5
% Brink/Satcher, Angular Momentum, p.24, Table I
% Richard N. Zare, Angular Momentum, p.89, Table 3.1
d2(1,1) = cos(beta/2)^3;
d2(4,4) = d2(1,1);
d2(1,2) = -sqrt(3)*cos(beta/2)^2*sin(beta/2);
d2(3,4) = d2(1,2);
d2(2,1) = -d2(1,2);
d2(4,3) = -d2(1,2);
d2(1,3) = sqrt(3)*cos(beta/2)*sin(beta/2)^2;
d2(3,1) = d2(1,3);
d2(2,4) = d2(1,3);
d2(4,2) = d2(1,3);
d2(1,4) = -sin(beta/2)^3;
d2(4,1) = -d2(1,4);
d2(2,2) = cos(beta/2)*(3*cos(beta/2)^2-2);
d2(3,3) = d2(2,2);
d2(2,3) = sin(beta/2)*(3*sin(beta/2)^2-2);
d2(3,2) = -d2(2,3); % wrong in Zare 1st printing, correct in 2nd

err = ~areequal(d1,d2,1e-10,'abs');
data = [];