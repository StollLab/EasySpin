function ok = test()

% check whether rescale works with minmax option

N = 100;
y = rand(1,N);
y_ = rescale(y,'minmax');

thr = 1e-12;
ok = areequal(min(y_),0) && areequal(max(y_),1,thr,'abs');
