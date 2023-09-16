function ok = test()

% Check offset correction along rows and data against mean

M = 10;
N = 201;
rng(3434);
data = rand(M,N);

n = 0;  % polynomial order

% Offset correction for each column
dim = 1;
[~,baseline] = basecorr(data,dim,n);
mean_ref = mean(data,dim);
ok(1) = areequal(mean_ref,baseline(1,:),1e-10,'abs');

% Offset correction for each row
dim = 2;
[~,baseline] = basecorr(data,dim,n);
mean_ref = mean(data,dim);
ok(2) = areequal(mean_ref,baseline(:,1),1e-10,'abs');
