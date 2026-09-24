function ok = test()

% Check that tensor frames are read correctly: Compare two ORCA calculations
% of triplet formaldehyde, with the coordinates rigidly rotated between the
% two. The rotation matrix is determined from the atom coordinates, and all
% tensors (g, D, A, Q) must transform with this rotation.

mol = 'tripletformaldehyde';
exts = {'.out','.property.txt'};

% Tolerances
tol.rmsd = 1e-5;  % Angstrom
tol.g = 1e-5;
tol.D = 2;  % MHz
tol.A = 0.2;  % MHz
tol.Q = 2e-3;  % MHz

T = @(v,ang) erot(ang).'*diag(v)*erot(ang);
maxdiff = @(X,Y) max(abs(X(:)-Y(:)));

for e = numel(exts):-1:1
  Sys0 = orca2easyspin(['orca/v6.1.0/' mol exts{e}]);
  Sys1 = orca2easyspin(['orca/v6.1.0/' mol '_tilted' exts{e}]);

  % Determine rotation R that maps untilted onto tilted coordinates
  [R,rmsd] = kabsch(Sys0.xyz,Sys1.xyz);
  okk = rmsd<tol.rmsd;

  % Assert that all tensors transform with R
  rotate = @(X) R*X*R.';
  okk = okk && maxdiff(rotate(T(Sys0.g,Sys0.gFrame)),T(Sys1.g,Sys1.gFrame))<tol.g;
  okk = okk && maxdiff(rotate(T(Sys0.D,Sys0.DFrame)),T(Sys1.D,Sys1.DFrame))<tol.D;
  okk = okk && isequal(Sys0.NucsIdx,Sys1.NucsIdx);
  for n = 1:numel(Sys0.NucsIdx)
    okk = okk && maxdiff(rotate(T(Sys0.A(n,:),Sys0.AFrame(n,:))),T(Sys1.A(n,:),Sys1.AFrame(n,:)))<tol.A;
    okk = okk && maxdiff(rotate(T(Sys0.Q(n,:),Sys0.QFrame(n,:))),T(Sys1.Q(n,:),Sys1.QFrame(n,:)))<tol.Q;
  end
  ok(e) = okk;
end

end

%-------------------------------------------------------------------------------
% Kabsch algorithm: proper rotation R that minimizes the RMSD between the
% centered coordinates x1*R.' and x2 (both nAtoms x 3)
function [R,rmsd] = kabsch(x1,x2)
x1 = x1 - mean(x1,1);
x2 = x2 - mean(x2,1);
[U,~,V] = svd(x1.'*x2);
R = V*diag([1 1 sign(det(V*U.'))])*U.';
rmsd = sqrt(mean(sum((x1*R.'-x2).^2,2)));
end
