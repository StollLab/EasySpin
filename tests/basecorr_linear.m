function ok = test()

% Linear baseline correction of linear data

N = 201;
Y = 1.32 + 0.234*linspace(0,1,N);
Ycorr = basecorr(Y,2,1);

ok = all(abs(Ycorr)<=1e-10);
