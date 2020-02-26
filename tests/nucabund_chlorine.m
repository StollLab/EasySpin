function ok = test()

% Cl isotopes

w1 = nucabund('35Cl,37Cl');
w2 = [0.7578 0.2422];
ok = areequal(w1,w2,1e-10,'rel');
