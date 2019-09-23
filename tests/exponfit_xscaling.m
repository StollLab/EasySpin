function [err,data] = test(opt,olddata)

% Test whether fitting procedure is independent of x axis scaling
%-------------------------------------------------------------------------------
% Data contains x = linspace(0,300,300)*1e-9

load('exponfit_xscaling.mat');

[k,c,yFit0] = exponfit(x, yData);
yFit = c(1) + c(2)*exp(-k*x);

x2 = x/max(x);
[k2,c,yFit0] = exponfit(x2, yData);
yFit2 = c(1) + c(2)*exp(-k2*x2);

if opt.Display
  subplot(3,1,[1 2]);
  plot(x,yData,'.',x,yFit,x,yFit2);
  legend('data','fit orig x-scale','fit rescaled x');
  subplot(3,1,3);
  plot(x,yFit-yFit2);
  title('Residuals');
end

err = ~areequal(yFit,yFit2,1e-10,'rel');

data = [];
