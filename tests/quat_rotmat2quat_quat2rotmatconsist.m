function ok = test()

% Check that rotmat2quat and quat2rotmat yield identical results

% Generate random axis-angle representation
v = rand(3,4,5);
v = bsxfun(@rdivide,v,sqrt(sum(v.*v,1)));
theta = 2*pi*rand(1,4,5);

matIn = zeros(3,3,4,5);

for k=1:5
  for n=1:4
    iv = v(:,n,k);
    itheta = theta(:,n,k);
   ex = [     0, -iv(3),  iv(2);
          iv(3),      0, -iv(1);
         -iv(2),  iv(1),      0];
   matIn(:,:,n,k) = eye(3)*cos(itheta) ...
                    + (1-cos(itheta))*(iv*iv.') + ex*sin(itheta);
  end
end

% Convert to quaternion
q = rotmat2quat(matIn);

% Convert back to Euler angles
matOut = quat2rotmat(q);

diff = abs(matOut-matIn);

% Check for consistency
ok = all(diff(:)<1e-14);
