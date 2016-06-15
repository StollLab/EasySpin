function [err,data] = test(opt,olddata)

% consistency with erot()
%======================================================

N = 500;
phi = rand(N,1)*2*pi;
theta = rand(N,1)*2*pi;
for k = 1:N
  n1 = ang2vec(phi(k),theta(k));
  R = erot(0,-theta(k),-phi(k));
  n2 = R(:,3);
  neq(k) = ~areequal(n1,n2);
end

err = any(neq);

data = [];
