function ok = test()

% Save and load 1D dataset with linear X axis

fn = tempname;

N = 1001;
x = linspace(300,400,N).';
y = rand(N,1);
eprsave(fn,x,y);
[xo,yo] = eprload(fn);
delete([fn '.*']);

ok = areequal(x,xo,1e-5,'rel') && areequal(y,yo,1e-5,'rel');
