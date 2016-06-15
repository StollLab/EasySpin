function [err,data] = test(opt,olddata)

% Test 2: result
%====================================================
N = 1001;
x = 1:N; x0 = 400; w = 10;
y = zeros(1,N);
y(x0) = 1;
z1 = convspec(y,1,w);
z2 = gaussian(x,x0,w);

err =  ~areequal(z1,z2);
data = [];
