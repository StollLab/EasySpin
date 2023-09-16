function ok = test()

% baseline + corrected data = data

N = 121;

Y = rand(1,N);
[Ycorr,baseline] = basecorr(Y,2,3);
ok(1) = all(abs(Ycorr+baseline-Y)<=1e-10);

Y = rand(N,1);
[Ycorr,baseline] = basecorr(Y,1,3);
ok(2) = all(abs(Ycorr+baseline-Y)<=1e-10);
