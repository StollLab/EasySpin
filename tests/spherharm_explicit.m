function ok = test()

% Comparison of spherharm to explicit expressions

N = 100;
theta = linspace(0,pi,N);
phi = linspace(2*pi,N);

q = 1;
for L = 0:2
  for M = -L:L
    Y = explicitY(L,M,theta,phi);
    Z = spherharm(L,M,theta,phi);
    ok(q) = areequal(Y,Z,1e-10,'abs');
    q = q + 1;
  end
end

ok = all(ok);

return

function v = explicitY(L,M,theta,phi)

v = inf;

switch L
case 0
  switch M
  case 0, v = 1/2/sqrt(pi)*ones(size(theta));
  end
case 1
  switch M
  case -1, v = +1/2*sqrt(3/2/pi)*sin(theta).*exp(-1i*phi);
  case 0,  v = 1/2*sqrt(3/pi)*cos(theta);
  case 1,  v = -1/2*sqrt(3/2/pi)*sin(theta).*exp(+1i*phi);
  end
case 2
  switch M
  case -2, v = 1/4*sqrt(15/2/pi)*sin(theta).^2.*exp(-2i*phi);
  case -1, v = +1/2*sqrt(15/2/pi)*cos(theta).*sin(theta).*exp(-1i*phi);
  case 0,  v = 1/8*sqrt(5/pi)*(3*cos(2*theta)+1);
  case 1,  v = -1/2*sqrt(15/2/pi)*cos(theta).*sin(theta).*exp(+1i*phi);
  case 2,  v = 1/4*sqrt(15/2/pi)*sin(theta).^2.*exp(+2i*phi);
  end
end
