function ok = test()

% Linear baseline correction of linear data

N = 201;
rng(3434);
data = rand(N,N);

dim = 2;
order = 0;
[~,baseline] = basecorr(data,dim,order);
mean_ref = mean(data,dim);

ok = areequal(mean_ref,baseline(:,1),1e-10,'abs');
