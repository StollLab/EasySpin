function ok = test()

% check whether rescale works with 'area' option

N = 100;
y = 2 + rand(1,N);
y_ = rescaledata(y,'area');
yref = y/sum(y);

thr = 1e-12;
ok = areequal(y_,yref,thr,'abs');
ok = ok && sign(y(1))==sign(y_(1));
