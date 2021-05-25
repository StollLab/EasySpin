function ok = test()

% syntax check

y = zeros(1,100);
y(1) = -1;
y(2) = 2;
z = 2*y;

y_ = rescaledata(y,'minmax');
y_ = rescaledata(y,'maxabs');
y_ = rescaledata(y, 'none');

[y_,factors] = rescaledata(y, 'minmax');
[y_,factors] = rescaledata(y, 'maxabs');
[y_,factors] = rescaledata(y, 'none');

y_ = rescaledata(y,z,'minmax');
y_ = rescaledata(y,z,'maxabs');
y_ = rescaledata(y,z,'lsq');
y_ = rescaledata(y,z,'lsq0');
y_ = rescaledata(y,z,'lsq1');
y_ = rescaledata(y,z,'lsq2');
y_ = rescaledata(y,z,'none');

[y_,factors] = rescaledata(y,z,'minmax');
[y_,factors] = rescaledata(y,z,'maxabs');
[y_,factors] = rescaledata(y,z,'lsq');
[y_,factors] = rescaledata(y,z,'lsq0');
[y_,factors] = rescaledata(y,z,'lsq1');
[y_,factors] = rescaledata(y,z,'lsq2');
[y_,factors] = rescaledata(y,z, 'none');

ok = true;
