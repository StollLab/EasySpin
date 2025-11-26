function ok = test()

% S=1 system
S = 1;
[Sz,Sy] = sop(S,'z','y');
spinFreq = 10e3;  % MHz

mwFreq = 10e3;  % MHz
tp = 0.010;  % µs

% Calculate propagator
H0 = spinFreq*Sz;
H1 = 1/2/tp*Sy;
U = propint(H0,H1,tp,mwFreq);

Ureference = [
   0.5000 + 0.0006i  -0.7071 - 0.0004i   0.5000 + 0.0000i
   0.7071 + 0.0004i   0.0000 - 0.0000i  -0.7071 + 0.0004i
   0.5000 + 0.0000i   0.7071 - 0.0004i   0.5000 - 0.0006i
];

ok = areequal(U,Ureference,1e-4,'abs');
