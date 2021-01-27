function [err,data] = test(opt,olddata)

%======================================================
% Save and load 2D dataset with nonlinear X and Y axes
%======================================================

fn = tempname;

nx = 50;
ny = 10;
x_lin = 1:nx;
y_lin = 1:ny;
x_nonlin = cumsum(rand(1,nx));
y_nonlin = cumsum(rand(1,ny));

data = rand(nx,ny) + 1i*rand(nx,ny);

xy{1} = {x_lin,y_nonlin};
xy{2} = {x_nonlin,y_lin};
xy{3} = {x_nonlin,y_nonlin};

for k = 1:3
  xy_ = xy{k};
  eprsave(fn,{xy_{1},xy_{2}},data);
  [x0,data0] = eprload(fn);
  err = ~areequal(x0{1},xy_{1},1e-5) || ...
        ~areequal(x0{2},xy_{2},1e-5) || ...
        ~areequal(data,data0,1e-5);
end

delete([fn '.*']);

data = [];
