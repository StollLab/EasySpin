function ok = test()

% Assert that baseline + corrected data = data

N = 6345;
data = rand(1,N);
[corrData,bline] = basecorr(data,3,2);

ok = areequal(corrData+bline,data,1e-10,'abs');
