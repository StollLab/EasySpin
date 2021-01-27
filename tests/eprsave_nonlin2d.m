function ok = test()

%======================================================
% Save and load 2D dataset with nonlinear X and Y axes
%======================================================

fn = tempname;

% Generate linear and nonlinear X and Y axes
nx = 50;
ny = 10;
x_lin = linspace(330,350,nx);
y_lin = linspace(0,3,ny);
x_nonlin = cumsum(rand(1,nx));
y_nonlin = cumsum(rand(1,ny));

% Generate data
data = rand(nx,ny) + 1i*rand(nx,ny);

% Three different combinations of linear and nonlinear axes
xy{1} = {x_lin,y_nonlin};
xy{2} = {x_nonlin,y_lin};
xy{3} = {x_nonlin,y_nonlin};

for k = 1:3
  x = xy{k}{1};
  y = xy{k}{2};
  eprsave(fn,{x,y},data);
  [x0,data0] = eprload(fn);
  
  ok = areequal(x0{1}(:),x(:),1e-5,'rel') && ...
       areequal(x0{2}(:),y(:),1e-5,'rel') && ...
       areequal(data,data0,1e-5','rel');
end

delete([fn '.*']);
