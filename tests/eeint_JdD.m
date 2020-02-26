function ok = test()

% Compare ee and J/dip/dvec for a two-spin system

S = [1/2 1/2];
g = [2 2];

J = rand*1e5; % isotropic exchange, MHz
d = rand(1,3)*1e5; % antisymmetric exchange, MHz
D = [-1 -1 2]*1e5; % dipolar coupling, MHz

% 1) ee
Sys1.S = S;
Sys1.g = g;
Sys1.ee = J*eye(3) +  diag(D) + ...
  [0 d(3) -d(2); -d(3) 0 d(1); d(2) -d(1) 0];
Hee1 = eeint(Sys1);

% 2) J, dip, dvec
Sys2.S = S;
Sys2.g = g;
Sys2.J = J;
Sys2.dip = D;
Sys2.dvec = d;
Hee2 = eeint(Sys2);

% The two Hamiltonians should be _numerically_ identical
ok = all(Hee1(:)==Hee2(:));
