function ok = test()

% baseline + corrected data = data

N = 6345;
Y = rand(1,N);
[c,b] = basecorr(Y,2,3);

ok = all(abs(c+b-Y)<=1e-10);
