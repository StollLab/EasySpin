function ok = test()

% The decay factor is exp(-k*abs(t)), so doubling t must square the
% decay factor.

t = linspace(0.1,2,15);  % µs
conc = 200;  % µM
lambda = 0.6;

V1 = dipbackground(t,conc,lambda);
V2 = dipbackground(2*t,conc,lambda);

ok = areequal(V2,V1.^2,1e-10,'rel');

end
