% startvec takes the Hilbert space spin operator Sx as input and returns as
% output the vector representation of the Sx operator in Liouville space

function StartingVector = startvec_LMK(basis,Potential,SxH)

if isfield(basis,'jK') && ~isempty(basis.jK) && any(basis.jK)
  error('This function expects an LMK basis, without jK.');
end

L = basis.L;
M = basis.M;
K = basis.K;
nOriBasis = numel(L);

lambda = Potential.lambda;
Lp = Potential.L;
Mp = Potential.M;
Kp = Potential.K;
if any(Mp)
  error('Potential coefficients with non-zero M are not supported.');
end

if ~any(lambda)
  idx0 = find(L==0 & M==0 & K==0);
  if numel(idx0)~=1
    error('Exactly one orientational basis function with L=M=K=0 is allowed.');
  end
  nSpinBasis = numel(SxH);
  idx = (idx0-1)*nSpinBasis + (1:nSpinBasis);
  nBasis = nOriBasis*nSpinBasis;
  StartingVector = sparse(idx,1,SxH(:),nBasis,nBasis);
  StartingVector = StartingVector/sqrt(sum(StartingVector.^2));
  return
end

% set up starting vector in orientational basis
oriVector = zeros(nOriBasis,1);
for b = 1:numel(oriVector)
  L_ = L(b);
  M_ = M(b);
  K_ = K(b);
  if mod(L_,2)==0 && M_==0 && mod(K_,2)==0
    % numerically integrate
    fun = @(a,b,c) conj(wignerd([L_ M_ K_],a,b,c)) .* exp(-U(a,b,c)/2) .* sin(b);
    Int = sqrt((2*L_+1)/(8*pi^2)) * integral3(fun,0,2*pi,0,pi,0,2*pi);
    oriVector(b) = Int;
  end
end

% form starting vector in direct product basis
StartingVector = real(kron(oriVector,SxH(:)));
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
