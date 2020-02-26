function ok = test()

% Check that rotating vectors by quaternions is equivalent to using
% rotation matrices

% Generate some uniformly random, normalized quaternions
q = rand(4,5,10);
q = bsxfun(@rdivide,q,sqrt(sum(q.*q,1)));

% Generate some uniformly random, normalized vectors
v = rand(3,5,10);
v = bsxfun(@rdivide,v,sqrt(sum(v.*v,1)));


R = quat2rotmat(q);

% Rotate the vectors using rotation matrices generated from the quaternions
for i=1:5
  for j=1:10
    vrot(:,i,j) = R(:,:,i,j)*v(:,i,j);
  end
end

% Compare methods of rotating vectors
diff = vrot - quatvecmult(q,v);

ok = all(abs(diff(:)) < 1e-10);
