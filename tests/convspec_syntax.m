function ok = test()

% Test syntax

data1d = zeros(1,100);
data1d(40) = 1;
data2d = rand(100,100);

o = convspec(data1d,1,10);
o = convspec(data1d,1,10,1);
o = convspec(data1d,1,10,2);
o = convspec(data1d,1,10,0,0.3);
o = convspec(data2d,1,[3 6],[0 1],[0 1]);

ok = true;
