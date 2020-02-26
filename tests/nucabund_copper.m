function ok = test()

% Cu isotopes

w1 = nucabund('63Cu,65Cu');
w2 = [0.6917 0.3083];
ok = areequal(w1,w2,1e-10,'rel');
