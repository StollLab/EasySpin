function ok = test()

% Compare the results of all methods against ground truth

S = 1/2;
[Sz,Sy] = sop(S,'z','y');
nu0 = 10e3;  % MHz
theta = pi/2;  % pulse flip angle
nu_mw = 10e3;  % MHz
tp = 0.010;  % µs

nu1 = theta/pi/tp;
H0 = nu0*Sz;
H1 = nu1*Sy;

Uref = [0.707107 + 0.000442, -0.707106 - 0.000000i
        0.707106 + 0.000000i, 0.707107 - 0.000442i];

% Calculate propagator
methods = ["simple", "RK4"];
Opt.n = 800;
for k = 1:numel(methods)
  Opt.Method = methods(k);
  U = propint(H0,H1,tp,nu_mw,0,Opt);
  ok(k) = areequal(U,Uref,1e-3,'abs'); 
end
