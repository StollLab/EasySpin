function ok = test()

% Interface test

rng(2);
y = rand(1,500);
m = 5;
yy = datasmooth(y,m);
yy = datasmooth(y,m,'savgol');
yy = datasmooth(y,m,'binom');
yy = datasmooth(y,m,'flat');
yy = datasmooth(y,m,'savgol',3);
yy = datasmooth(y,m,'savgol',3,1);

ok = true;
