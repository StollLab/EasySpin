function ok = test()

%====================================================
% check whether rescale works with vectors of different size
%====================================================
y = zeros(1,100);
y(1) = -1;
y(2) = 2;
z1 = zeros(1,90);
z1(5) = 10;
z2 = zeros(90,1);
z2(5) = 10;

y_ = rescaledata(y,z1,'minmax');
y_ = rescaledata(y,z2,'minmax');

ok = true;
