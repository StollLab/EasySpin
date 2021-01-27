function ok = test()

% Save and load 1D dataset with nonlinear X axis

fn = tempname;

% Generate non-linear axis and 1D data
N = 11;
x = linspace(1,5,N).^2;
x = x(:);
data = rand(N,1);

% Save date
eprsave(fn,x,data);

% Load data
[x_,data_] = eprload(fn);
delete([fn '.*']);

ok = areequal(x,x_,1e-5,'rel') && areequal(data,data_,1e-5,'rel');
