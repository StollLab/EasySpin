function ok = test()

% Linear baseline correction of linear data

N = 1023;
Y = 1 + 0.234*linspace(0,1,N);
[c,~] = basecorr(Y,2,1);

ok = all(abs(c)<=1e-10);
