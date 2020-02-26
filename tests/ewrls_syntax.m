function ok = test()

%=======================================================
% Syntax tests
%=======================================================
nPoints = 100;
nScans = 10;
dataarray = rand(nPoints,nScans);

p = 10;
lambda = 0.98;
delta = 10;
nPreAvg = 5;

y = ewrls(dataarray,p,lambda);
y = ewrls(dataarray,p,lambda,0);
y = ewrls(dataarray,p,lambda,nPreAvg);
y = ewrls(dataarray,p,lambda,0,delta);
y = ewrls(dataarray,p,lambda,0,delta,'f');
y = ewrls(dataarray,p,lambda,0,delta,'b');
y = ewrls(dataarray,p,lambda,0,delta,'fb');

ok = true;
