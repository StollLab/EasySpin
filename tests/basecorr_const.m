function ok = test()

% Linear baseline correction of linear data

N = 1023;
rng(3434);
data = rand(N,N);

[~,baseline] = basecorr(data,0,2);
mean2 = mean(data,2);

ok = areequal(mean2,baseline(:,1),1e-10,'abs');
