function [err,data] = test(opt,olddata)

% Compare plegendre results to explicit expressions
%============================================================

x = 2*rand(1,100)-1;
x2 = x.^2;

y0 = plegendre(5,3,x);
y1 = 105/2*(1-x2).^(3/2).*(9*x2-1);
err(1) = ~areequal(y0,y1,1e-10);

y0 = plegendre(3,1,x);
y1 = 3/2*sqrt(1-x2).*(5*x2-1);
err(2) = ~areequal(y0,y1,1e-10);

y0 = plegendre(9,2,x);
y1 = -495/16*(x-1).*x.*(x+1).*(221*x2.^3 -273*x2.^2 +91*x2-7);
err(3) = ~areequal(y0,y1,1e-10);

err = any(err);

data = [];
