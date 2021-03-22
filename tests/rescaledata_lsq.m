function ok = test()

% check whether rescale works with various least-square fitting option

thr = 1e-12;

N = 100;
x = linspace(0,1,N);
yr = rand(1,N);

ys = 2.1*yr;
y0 = ys + 0.5;
y1 = y0 + 1.2*x;
y2 = y1 - 2.1*x.^2;
y3 = y2 + 0.7*x.^3;

ys_ = rescaledata(ys,yr,'lsq');
y0_ = rescaledata(y0,yr,'lsq0');
y1_ = rescaledata(y1,yr,'lsq1');
y2_ = rescaledata(y2,yr,'lsq2');
y3_ = rescaledata(y3,yr,'lsq3');

ok(1) = areequal(ys_,yr,thr,'abs');
ok(2) = areequal(y0_,yr,thr,'abs');
ok(3) = areequal(y1_,yr,thr,'abs');
ok(4) = areequal(y2_,yr,thr,'abs');
ok(5) = areequal(y3_,yr,thr,'abs');
