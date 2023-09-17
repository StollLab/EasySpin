function ok = test()

x = linspace(0,1,201);
y0 = gaussian(x,0.3,0.1);
y0 = y0/max(y0);

N = 10;
y0 = repmat(y0,N,1);

baseline = rand(N,1) + (2*rand(N,1)-1)*(x-0.5) + (2*rand(N,1)-1)*(x-0.5).^2;

y = y0 + baseline;
region = x<0.15 | x>0.5;

yc_rows = basecorr(baseline,2,2,region);
yc_cols = basecorr(baseline.',1,2,region);
yc_rows_region = basecorr(y,2,2,region);
yc_cols_region = basecorr(y.',1,2,region);

threshold = 1e-15;
ok(1) = areequal(yc_rows,0*y,threshold,'abs');
ok(2) = areequal(yc_cols.',0*y,threshold,'abs');
threshold = 1e-3;
ok(3) = areequal(yc_rows_region,y0,threshold,'abs');
ok(4) = areequal(yc_cols_region.',y0,threshold,'abs');
