function ok = test()

% Bin a lot of peaks

N = 20e3;
pos = rand(1,N); amp = rand(1,N);
[x,s] = makespec([-1,2],2000,pos,amp);

ok = abs(sum(amp)-sum(s))<N*1e-10;
