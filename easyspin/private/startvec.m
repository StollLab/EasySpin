% startvec takes the Hilbert space spin operator Sx as input and returns as
% output the vector representation of the Sx operator in Liouville space

function StartingVector = startvec(basis,lambda,SxH)

L = basis.L;
M = basis.M;
K = basis.K;
jK = basis.jK;
nOriBasis = numel(L);

if ~any(lambda)
  idx0 = find(L==0 & M==0 & K==0);
  if numel(idx0)~=0
    error('Exactly one orientational basis function with L=M=K=0 is allowed.');
  end
  nSpinBasis = numel(SxH);
  idx = (idx0-1)*nSpinBasis + (1:nSpinBasis);
  nBasis = nOriBasis*nSpinBasis;
  StartingVector = sparse(idx,1,SxH(:),nBasis,nBasis);
  StartingVector = StartingVector/norm(StartingVector);
  return
end

c20 = lambda(1);
c22 = lambda(2);
c40 = lambda(3);
c42 = lambda(4);
c44 = lambda(5);

% set up starting vector in spatial basis
oriVector = zeros(nOriBasis,1);
for b = 1:numel(oriVector)
  L_ = L(b);
  M_ = M(b);
  K_ = K(b);
  if mod(L_,2)==0 && M_==0 && mod(K_,2)==0
    % numerically integrate
    fun = @(a,b,c) wignerd([L_ M_ K_],a,b,c).*exp(-U(a,b,c)/2);
    oriVector(b) = integral3(fun,0,2*pi,0,pi,0,2*pi);
  end
end
oriVector = oriVector/norm(oriVector);

% set up spin part of basis
SxVector = full(SxH(:)/norm(SxH(:)));

% form starting vector in direct product basis
StartingVector = real(kron(SxVector,oriVector));
StartingVector = StartingVector/norm(StartingVector);
StartingVector = sparse(StartingVector);

  function u = U(a,b,c)
    u = 0;
    if c20~=0, u = u + wignerd([2 0 0],a,b,c); end
    if c22~=0, u = u + wignerd([2 2 0],a,b,c); end
    if c40~=0, u = u + wignerd([4 0 0],a,b,c); end
    if c42~=0, u = u + wignerd([4 2 0],a,b,c); end
    if c44~=0, u = u + wignerd([4 4 0],a,b,c); end
  end

end
