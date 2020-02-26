function ok = test()

% syntax check

y = zeros(1,100);
y(1) = -1;
y(2) = 2;
z = 2*y;

y_ = rescale(y,'minmax');
y_ = rescale(y,'maxabs');
y_ = rescale(y, 'none');

[y_,factors] = rescale(y, 'minmax');
[y_,factors] = rescale(y, 'maxabs');
[y_,factors] = rescale(y, 'none');

y_ = rescale(y,z,'minmax');
y_ = rescale(y,z,'maxabs');
y_ = rescale(y,z,'lsq');
y_ = rescale(y,z,'lsq0');
y_ = rescale(y,z,'lsq1');
y_ = rescale(y,z,'lsq2');
y_ = rescale(y,z,'none');

[y_,factors] = rescale(y,z,'minmax');
[y_,factors] = rescale(y,z,'maxabs');
[y_,factors] = rescale(y,z,'lsq');
[y_,factors] = rescale(y,z,'lsq0');
[y_,factors] = rescale(y,z,'lsq1');
[y_,factors] = rescale(y,z,'lsq2');
[y_,factors] = rescale(y,z, 'none');

ok = true;
