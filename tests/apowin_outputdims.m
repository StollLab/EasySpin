function ok = test()

% test dimension of output

N = 2345;
w = apowin('ham',N);
ok = all(size(w)==[N 1]);
