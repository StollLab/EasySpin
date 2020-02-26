function ok = test()

% Cl values

q = nucqmom('35Cl,37Cl');
q1 = [-0.0817  -0.0644 ];

ok = areequal(q,q1,1e-10,'rel');
