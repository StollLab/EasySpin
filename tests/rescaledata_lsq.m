function ok = test()

% check whether rescale works with various least-square fitting option

thr = 1e-12;

N = 100;
x = linspace(0,1,N);
yr = rand(1,N);

ys = 2.1*yr;

ys_ = rescaledata(ys,yr,'lsq');

ok(1) = areequal(ys_,yr,thr,'abs');
