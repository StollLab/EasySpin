% isto  Irreducible spherical tensor operators (ISTOs)
%
%   T = isto(J,kq)
%   T = isto(J,kqi)
%   T = isto(___,'sparse')
%
% Input:
%   J    array of angular momentum quantum numbers (0, 1/2, 1, etc.);
%   kq   array of ranks k and projections q for each of the angular momenta in J
%        each row contains [k q] for the corresponding J
%   kqi  array of ranks k, projections q, and ang.mom. indices i; each row
%        contains [k q i]. For any ang.mom. not listed, k=0 and q=0 are assumed.
%
% Output:
%   T    irreducible spherical tensor operator matrix
%
% If 'sparse' is given, the matrix is returned in sparse storage format,
% otherwise it is returned in full format.
%
% Examples:
%   T = isto(3/2,[2 0])            % T^2_0 for a spin 3/2
%   T = isto([1/2 1/2],[1 1 2])    % T^1_1 for the second spin

function T = isto(J,kq,sparseOpt)

% see
%   I. D. Ryabov
%   On the Generation of Operator Equivalents and the Calculation of Their Matrix Elements
%   J. Magn. Reson. 140, 141-145 (1999)
%   https://doi.org/10.1006/jmre.1999.1783

if nargin<3
  sparseOpt = 'full';
end
sparseOutput = strcmp(sparseOpt,'sparse');

% Check J
%-------------------------------------------------------------------------------
if isstruct(J)
  if ~isfield(J,'S')
    error('Sys.S is missing.');
  end
  J = J.S;
elseif ~isnumeric(J) || any(mod(2*J,1)) || any(J<0)
  error('First input must be a spin system structure, or an array of angular momentum quantum numbers (0, 1/2, 1, ...).');
end

% Check kq
%-------------------------------------------------------------------------------
if size(kq,2)==2
  if size(kq,1)~=numel(J)
    error('Provide k and q for each angular momentum.');
  end
elseif size(kq,2)==3
  kqi = kq;
  k = kqi(:,1);
  q = kqi(:,2);
  i = kqi(:,3);
  if any(i>numel(J))
    error('index in kqi is larger than number of angular momenta.');
  end
  kq = zeros(numel(J),2);
  kq(i,1) = k;
  kq(i,2) = q;
end
k = kq(:,1);
q = kq(:,2);

if any(abs(q)>k)
  error('q must satisfy -k<=q<=k.');
end
if any(mod(k,1)) || any(k<0)
  error('k must be a non-negative integer (0, 1, 2, ...).');
end
if any(mod(q,1))
  error('q must be a integer.');
end

% Calculate operator
%-------------------------------------------------------------------------------
T = sparse(1);
for j = 1:numel(J)
  Tj = isto_single(J(j),k(j),q(j));
  T = kron(T,Tj);
end

% Sparse vs full output
%-------------------------------------------------------------------------------
if ~sparseOutput
  T = full(T);
end

end
%===============================================================================

%-------------------------------------------------------------------------------
% Calculate ISTO for a single angular momentum; no input checks
function T = isto_single(J,k,q)

% Special case: no angular momentum
if J==0, T = sparse(1); return; end

% Get ladder operators
Jp = sop(J,'+','sparse');
Jm = sop(J,'-','sparse');

% Start with T^k_k (rank k, highest projection k)
Nkk = (-1)^k/2^(k/2); % overall normalization constant (Ryabov Eq.[17])
T = Nkk*Jp^k; % (Ryabov, after Eq.[2])

% Apply Racah commutation rule repeatedly to lower projection down to |q|
for q_ = k:-1:abs(q)+1
  T = (Jm*T - T*Jm)/sqrt(k*(k+1)-q_*(q_-1)); % (Ryabov Eq.[1])
end

% Apply symmetry transformation for negative q
if q<0
  T = (-1)^q*T'; % (Ryabov Eq.[3]);
end

end
