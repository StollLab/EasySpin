function [err,data] = test(opt,olddata)

% Compute Clebsch-Gordan matrix
% for transformation from uncoupled to coupled representation.
% See Richard N. Zare, Angular Momentum, p.45

clear
%j1 = 3; j2 = 3;
j1 = 5; j2 = j1;

N = (2*j1+1)*(2*j2+1);

%{
[m1,m2] = meshgrid(j1:-1:-j1,j2:-1:-j2);
m12 = sortrows([-m1(:)-m2(:) m1(:) m2(:)]);
m1 = m12(:,2);
m2 = m12(:,3);

k = 1;
for j_=j1+j2:-1:abs(j1-j2)
  for m_=j_:-1:-j_
    j(k) = j_;
    m(k) = m_;
    k = k+1;
  end
end
jm = sortrows([-m(:) m(:) j(:)]);
m = jm(:,2);
j = jm(:,3);

CG = zeros(N);
for r = 1:N
  for c = 1:N
    CG(r,c) = clebschgordan(j1,j2,j(r),m1(c),m2(c),m(r));
  end
end
%}

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

deviation = abs(Id-eye(N));
err = any(deviation(:)>1e-13);

data = [];

