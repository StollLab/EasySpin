function ok = test()

% non negative check

y = ones(1,100);
y(1) = -1;
y(2) = 2;
z = -2*y;

[y_,scale] = rescaledata(y,z,'lsq');

ok = all(scale>0);
