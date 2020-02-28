function ok = test()

% Check nuclear-nuclear Hamiltonian against explicit calculation

Sys.Nucs = '1H,13C';
Sys.A = [10 3];

[I1{1},I1{2},I1{3}] = sop(Sys,'x2','y2','z2');
[I2{1},I2{2},I2{3}] = sop(Sys,'x3','y3','z3');

nn{1} = 3;
nn{2} = [4 6 7];
nn{3} = [1 2 3; 4 5 6; 7 8 9];

for k = 1:3
  Sys.nn = nn{k};
  H0 = nnint(Sys);
  
  switch numel(nn{k})
    case 1
      nnMatrix = nn{k}*eye(3);
    case 3
      nnMatrix = diag(nn{k});
    case 9
      nnMatrix = nn{k};
  end
  H = 0;
  for c1 = 1:3
    for c2 = 1:3
      H = H + nnMatrix(c1,c2)*I1{c1}*I2{c2};
    end
  end
  
  ok(k) = areequal(H0,H,1e-10,'abs');
  
end
