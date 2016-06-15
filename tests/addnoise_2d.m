function [err,data] = test(opt,olddata)

%======================================================
% 2D data sets
%======================================================

SNR = 20;

y = rand(14,1653);

model = 'fun';
for k=1:3
  y_ = addnoise(y,SNR,model(k));
  ok(k) = all(size(y_)==size(y));
end

err = any(~ok);

data = [];
