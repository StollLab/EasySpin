% startvec takes the Hilbert space spin operator Sx as input and returns as
% output the vector representation of the Sx operator in Liouville space

function StartingVector = startvec(basis,Potential,SopH)

L = basis.L;
M = basis.M;
K = basis.K;
jK = basis.jK;
nOriBasis = numel(L);

lambda = Potential.lambda;
Lp = Potential.L;
Mp = Potential.M;
Kp = Potential.K;

% Detect Freed-style potential (even L, M=0, even K)
zeroM = all(Mp==0);
evenL = all(mod(Lp,2)==0);
evenK = all(mod(Kp,2)==0);

if ~any(lambda)
  idx0 = find(L==0 & M==0 & K==0);
  if numel(idx0)~=1
    error('Exactly one orientational basis function with L=M=K=0 is allowed.');
  end
  nSpinBasis = numel(SopH);
  idx = (idx0-1)*nSpinBasis + (1:nSpinBasis);
  nBasis = nOriBasis*nSpinBasis;
  StartingVector = sparse(idx,1,SopH(:),nBasis,nBasis);
  StartingVector = StartingVector/sqrt(sum(StartingVector.^2));
  return
end

% set up starting vector in orientational basis
oriVector = zeros(nOriBasis,1);
for b = 1:numel(oriVector)
  L_  = L(b);
  M_  = M(b);
  K_  = K(b);
  jK_ = jK(b);
  
  if zeroM
    if M_~=0, continue; end
    if evenL && mod(L_,2)~=0, continue; end
    if evenK && mod(K_,2)~=0, continue; end
    if jK_~=1, continue; end
    fun = @(b,c) cos(K_*c) .* wignerd([L_ 0 K_],b) .* exp(-U(0,b,c)/2) .* sin(b);
    Int = (2*pi) * integral2(fun,0,pi,0,2*pi);
  else
    fun = @(a,b,c) conj(wignerd([L_ M_ K_],a,b,c)) .* exp(-U(a,b,c)/2) .* sin(b);
    Int = integral3(fun,0,2*pi,0,pi,0,2*pi);
  end
  
  if abs(Int)<1e-10, continue; end
  
  oriVector(b) = sqrt(2/(1 + (K_==0))) * sqrt((2*L_+1)/(8*pi^2)) * Int;
end

% form starting vector in direct product basis
StartingVector = real(kron(oriVector,SopH(:)));
StartingVector = StartingVector/norm(StartingVector);
StartingVector = sparse(StartingVector);

  function u = U(a,b,c)
    u = 0;
    for p = 1:numel(lambda)
      if lambda(p)==0, continue; end
      if Kp(p)==0
        u = u - wignerd([Lp(p) 0 0],a,b,c) * lambda(p);
      else
        u = u - wignerd([Lp(p) 0 +Kp(p)],a,b,c) * lambda(p) ...
              - wignerd([Lp(p) 0 -Kp(p)],a,b,c) * lambda(p)' * (-1)^Kp(p);
      end
    end
  end

end
