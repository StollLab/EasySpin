function [err,data] = test(opt,olddata)

% Test 1: syntax
%====================================================
y = zeros(1,100);
y(40) = 1;
y2 = rand(100,100);

o = convspec(y,1,10);
o = convspec(y,1,10,1);
o = convspec(y,1,10,2);
o = convspec(y,1,10,0,0.3);
o = convspec(y2,1,[3 6],[0 1],[0 1]);

err = 0;
data = [];
