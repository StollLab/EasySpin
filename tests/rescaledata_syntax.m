function ok = test()

% syntax check

y = zeros(1,100);
y(1) = -1;
y(2) = 2;
z = 2*y;

y_ = rescaledata(y,'maxabs');
y_ = rescaledata(y, 'none');

[y_,scale] = rescaledata(y, 'maxabs');
[y_,scale] = rescaledata(y, 'none');

y_ = rescaledata(y,z,'maxabs');
y_ = rescaledata(y,z,'lsq');
y_ = rescaledata(y,z,'none');

[y_,scale] = rescaledata(y,z,'maxabs');
[y_,scale] = rescaledata(y,z,'lsq');
[y_,scale] = rescaledata(y,z, 'none');

ok = true;
