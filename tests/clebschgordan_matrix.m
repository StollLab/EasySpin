function ok = test()

% Compute Clebsch-Gordan matrix
% for transformation from uncoupled to coupled representation.
% See Richard N. Zare, Angular Momentum, p.45

j1 = 5;
j2 = j1;

N = (2*j1+1)*(2*j2+1);

CG = [];
for m = j1+j2:-1:-(j1+j2)  
  j = j1+j2:-1:max(abs(m),abs(j1-j2));
  m1 = j1:-1:-j1;
  m1(abs(m-m1)>j2)=[];
  m2 = m-m1;
  cg_ = [];
  for q = 1:numel(m1)
    for k = 1:numel(j)
      cg_(q,k) = clebschgordan(j1,j2,j(k),m1(q),m2(q),m);
    end
  end
  CG = blkdiag(CG,cg_);
end

% Check orthogonality of Clebsch-Gordan matrix
Id = CG'*CG;
ok = areequal(Id,eye(N),1e-13,'abs');
