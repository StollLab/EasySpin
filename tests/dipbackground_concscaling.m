function ok = test()

% The decay rate k is linear in concentration, so doubling the
% concentration must square the decay factor.

t = linspace(0.1,3,15);  % µs
conc = 80;  % µM
lambda = 0.4;

V1 = dipbackground(t,conc,lambda);
V2 = dipbackground(t,2*conc,lambda);

ok = areequal(V2,V1.^2,1e-10,'rel');

end
