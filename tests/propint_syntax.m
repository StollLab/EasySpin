function ok = test()

% Test whether all input argument combinations are handled without crash

n = 4;
H0 = complex(rand(n),rand(n)); H0 = H0+H0';
H1 = complex(rand(n),rand(n)); H1 = H1+H1';

t = [0 10];  % µs
freq = 0.5768;  % MHz
phase = deg2rad(20);  % rad

Opt.n = 400;
Opt.Method = 'simple';

U = propint(H0,H1,t,freq);
U = propint(H0,H1,t,freq,phase);
U = propint(H0,H1,t,freq,phase,Opt);
[U,Ustore] = propint(H0,H1,t,freq);
U = propint(Ustore,t,freq);

ok = true;
