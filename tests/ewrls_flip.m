function ok = test(opt)

%=======================================================
% Test foward and backward filtering
%=======================================================
nPoints = 100;
nScans = 10;
dataarray = rand(nPoints,nScans);

p = 10;
lambda = 0.98;
delta = 10;

yf = ewrls(dataarray,p,lambda,0,delta,'f');
yb = ewrls(dataarray(end:-1:1,:),p,lambda,0,delta,'b');
yb = yb(end:-1:1);

ok = areequal(yf,yb,1e-7,'rel');

if opt.Display
  subplot(2,1,1);
  plot(1:nPoints,yf,1:nPoints,yb,'r');
  legend('forward','backward'); legend boxoff
  subplot(2,1,2);
  plot(1:nPoints,yb-yf);
  legend('backward-forward'); legend boxoff
end
