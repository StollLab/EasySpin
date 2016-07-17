function [err,data] = test(opt,olddata)

%======================================================
% rotaxi2mat should the exact inverse of rotmat2axi
%======================================================
Angles = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; ...
  1 1 1; 0.5 0 0; 0 0.5 0; 0 0 0.5; 0.5 0.5 0; 0.5 0 0.5; ...
  0 0.5 0.5; 0.5 0.5 0.5; randn(200,3)]*pi;
for k = 1:size(Angles,1)
  R1 = erot(Angles(k,:));
  [n,phi] = rotmat2axi(R1);
  R2 = rotaxi2mat(n,phi);
  ok(k) = areequal(R1,R2,1e-8);
end

err = any(~ok);
data = [];
