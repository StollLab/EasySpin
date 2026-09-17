function ok = test()

% The decay rate k is linear in the modulation depth lambda, so
% doubling lambda must square the decay factor.

t = linspace(0.1,3,15);  % µs
conc = 120;  % µM
lambda = 0.3;

V1 = dipbackground(t,conc,lambda);
V2 = dipbackground(t,conc,2*lambda);

ok = areequal(V2,V1.^2,1e-10,'rel');

end
