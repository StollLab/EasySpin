function [err,data] = test(opt,olddata)

% Test 5: input is a matrix
%================================================================
y = rand(100,5);
a = deriv(y);
err = ~areequal(size(a),size(y),0);
data = [];
