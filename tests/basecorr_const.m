function ok = test()

% Linear baseline correction of linear data

N = 1023;
rng(3434);
data = rand(N,N);

[~,baseline] = basecorr(data,2,0);
mean2 = mean(data,2);

ok = areequal(mean2,baseline(:,1),1e-10,'abs');
