function [err,data] = test(opt,olddata)

%======================================================
% Comparison of wignerd output with explicit expression
%======================================================

beta = rand*pi;

d1 = wignerd(2,beta,'-');

% '-' definition:
% Varshalovich, Quantum Theory of Angular Momentum, p.199, Table 4.6
% Brink/Satcher, Angular Momentum, p.24, Table I
% Richard N. Zare, Angular Momentum, p.89, Table 3.1
d2(1,1) = cos(beta/2)^4;
d2(5,5) = d2(1,1);
d2(1,2) = -1/2*sin(beta)*(1+cos(beta));
d2(2,1) = -d2(1,2);
d2(5,4) = -d2(1,2);
d2(4,5) = d2(1,2); % wrong sign in Zare 1st printing, correct in 2nd
d2(1,3) = sqrt(3/8)*sin(beta)^2;
d2(3,1) = d2(1,3);
d2(5,3) = d2(1,3);
d2(3,5) = d2(1,3);
d2(1,4) = 1/2*sin(beta)*(cos(beta)-1);
d2(2,5) = d2(1,4);
d2(5,2) = -d2(1,4);
d2(4,1) = -d2(1,4);
d2(1,5) = sin(beta/2)^4;
d2(5,1) = d2(1,5);
d2(2,2) = 1/2*(2*cos(beta)-1)*(cos(beta)+1);
d2(4,4) = d2(2,2);
d2(2,4) = 1/2*(2*cos(beta)+1)*(1-cos(beta));
d2(4,2) = d2(2,4);
d2(2,3) = -sqrt(3/2)*sin(beta)*cos(beta); % wrong sign in Zare 1st printing, correct in 2nd
d2(3,4) = d2(2,3);
d2(3,2) = -d2(2,3);
d2(4,3) = -d2(2,3);
d2(3,3) = 1/2*(3*cos(beta)^2-1);

err = ~areequal(d1,d2,1e-10);
data = [];