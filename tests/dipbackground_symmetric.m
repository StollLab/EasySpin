function ok = test()

% dipbackground uses abs(t) internally, so Vinter must be an even
% function of t.

t = linspace(-3,3,25);  % µs
conc = 150;  % µM
lambda = 0.5;

V1 = dipbackground(t,conc,lambda);
V2 = dipbackground(-t,conc,lambda);

ok = areequal(V1,V2,1e-12,'abs');

end
