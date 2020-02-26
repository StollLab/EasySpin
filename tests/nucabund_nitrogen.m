function ok = test()

% Nitrogen isotopes

w1 = nucabund('14N,15N');
w2 = [0.99632 0.00368];
ok = areequal(w1,w2,1e-10,'rel');
