function ok = test()

% Quadratic baseline correction of quadratic data

N = 634;
X = linspace(0,1,N);
Y = 1 - 0.234*X + 0.245*X.^2;
[Ycorr,~] = basecorr(Y,2,2);

ok = all(abs(Ycorr)<=1e-10);
