% startvec takes the Hilbert space spin operator Sx as input and returns as
% output the vector representation of the Sx operator in Liouville space

function StartingVector = startvec(Basis,lambda,SxH)

c20 = lambda(1);
c22 = lambda(2);
c40 = lambda(3);
c42 = lambda(4);
c44 = lambda(5);

basisList = Basis.List;

if any(lambda)
  % set up Wigner part of basis
  wigVector = zeros(1,size(basisList,1));
  L = basisList(:,1);
  M = basisList(:,2);
  K = basisList(:,3);
  for iWig = 1:numel(wigVector)
    L_ = L(iWig);
    M_ = M(iWig);
    K_ = K(iWig);
    LKM = [L_ M_ K_];
    if mod(L_,2)==0 && M_==0 && mod(K_,2)==0
      % numerically integrate
      fun = @(a,b,c) wignerd(LKM,a,b,c).*exp(-U(a,b,c)/2);
      wigVector(iWig) = integral3(fun,0,2*pi,0,pi,0,2*pi);
    else
      wigVector(iWig) = 0;
    end
  end
  wigVector = wigVector/norm(wigVector);
  
  % set up spin part of basis
  SxVector = full(SxH(:)/norm(SxH(:)));
  
  % form starting vector in direct product basis
  %StartingVector = real(kron(wigVector,SxVector));
  StartingVector = real(kron(SxVector,wigVector));
  StartingVector = StartingVector/norm(StartingVector);
  StartingVector = sparse(StartingVector);
  
else
  % assumes L=M=K=0 is the first spatial basis function!
  SxVector = SxH(:)/norm(SxH(:));
  nSpin = numel(SxVector);
  nBasis = nSpin*size(basisList,1);
  StartingVector = sparse(1:nSpin,1,SxVector,nBasis,nBasis);
end

return

function u = U(a,b,c)
  u = 0;
  if c20~=0, u = u + wignerd([2 0 0],a,b,c); end
  if c22~=0, u = u + wignerd([2 2 0],a,b,c); end
  if c40~=0, u = u + wignerd([4 0 0],a,b,c); end
  if c42~=0, u = u + wignerd([4 2 0],a,b,c); end
  if c44~=0, u = u + wignerd([4 4 0],a,b,c); end
end

end
