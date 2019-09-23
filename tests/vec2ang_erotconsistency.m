function [err,data] = test(opt,olddata)

%======================================================
% Test 4: compatibility with erot()
%======================================================
N = 500;
phi = rand(1,N);
theta = rand(1,N);
neq = [];
for k=1:N
  R = erot(phi(k),theta(k),0);
  [p(k),t(k)] = vec2ang(R(3,:).');
end
err = ~areequal(phi,p,1e-12,'abs') || ~areequal(theta,t,1e-12,'abs');
data = [];
