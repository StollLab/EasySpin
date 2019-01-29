% Calculate transformation matrix for the K-symmetrization of the
% Liouvillian matrix, to render it complex symmetric for the Lanczos
% algorithm.
%
% To obtain the complex symmetric form of L and the starting vector,
% use
%    sv_symm = T'*sv
%    L_symm = T'*L*T

function T = ksymmetrizer(basis)

DebugMode = true;

if isfield(basis,'jK')
  error('This function expects an LMK basis, without jK.');
end

negjKtrim = isfield(basis,'jKmin') && basis.jKmin==1;


% Set up basis indices for LMK and LMKjK bases
%-------------------------------------------------------------------
L = basis.L;
M = basis.M;
K = basis.K;
nBasis1 = numel(L);

L_ = L;
M_ = M;
K_ = abs(K);
jK_ = sign(K);
idxKzero = K_==0;
jK_(idxKzero) = (-1).^L(idxKzero);

if negjKtrim
  jKneg = jK_<0;
  L_(jKneg) = [];
  M_(jKneg) = [];
  K_(jKneg) = [];
  jK_(jKneg) = [];
end
nBasis2 = numel(L_);


% Calculate matrix elements of T: <LMK|LMKjK>
%-------------------------------------------------------------------
idx = 0;
for b1 = 1:nBasis1
  L1 = L(b1);
  M1 = M(b1);
  K1 = K(b1);
  
  for b2 = 1:nBasis2
    L2 = L_(b2);
    if L1~=L2, continue, end
    M2 = M_(b2);
    if M1~=M2, continue, end
    K2 = abs(K_(b2));
    jK = jK_(b2);
    
    % Calculate matrix element
    v = sqrt(jK)/sqrt(2*(1+(K2==0))) * ((K1==K2) + jK*(-1)^(L2+K2)*(K1==-K2));
    
    if v==0, continue; end
        
    % Store matrix element
    idx = idx + 1;
    bra(idx) = b1;
    ket(idx) = b2;
    vals(idx) = v;
    
  end
end

% Build sparse matrix
T = sparse(bra,ket,vals,nBasis1,nBasis2);

% Assert the matrix is a transformation matrix
if DebugMode
  err = T'*T - speye(nBasis2);
  if max(abs(err(:)))>1e-10
    error('Transformation matrix is not unitary.');
  end
end
