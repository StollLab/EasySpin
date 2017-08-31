% nnint  Nuclear-nuclear spin interaction Hamiltonian 
%
%   F = nnint(SpinSystem)
%   F = nnint(SpinSystem,nSpins)
%   F = nnint(SpinSystem,nSpins,'sparse')
%
%   Returns the nuclear-nuclear spin interaction (NNI)
%   Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. NNI
%       parameters are in the nn and nnFrame fields.
%   - nSpins: If given, specifies nuclear spins
%       for which the NNI should be computed. If
%       absent, all nuclei are included.
%   - 'sparse': If given, the matrix is returned in sparse format.

%   Output:
%   - F: Hamiltonian matrix containing the NNI for
%       nuclear spins specified in eSpins.

function F = nnint(System,Spins,opt)

if (nargin==0), help(mfilename); return; end

if (nargin<1) || (nargin>3), error('Wrong number of input arguments!'); end
if (nargout<0), error('Not enough output arguments.'); end
if (nargout>1), error('Too many output arguments.'); end

if nargin<3, opt = ''; end
if nargin<2, Spins = []; end
if ~ischar(opt)
  error('Third input must be a string, ''sparse''.');
end
sparseResult = strcmp(opt,'sparse');

[System,err] = validatespinsys(System);
error(err);

sys = spinvec(System);
n = prod(2*sys+1);

% Special cases: only one nuclear spin, nn not given, or all zero
F = sparse(n,n);
if System.nNuclei<2, return; end
if ~any(System.nn(:)), return; end

if isempty(Spins), Spins = 1:System.nNuclei; end

% Some error checking on the second input argument
if numel(Spins)<2
  error('Spins (2nd argument) must contain at least 2 values!'); 
end
if any(Spins<1) || any(Spins>System.nNuclei)
  error('Spins (2nd argument) must contain values between 1 and %d!',System.nNuclei);
end
if numel(unique(Spins))~=numel(Spins)
  error('Spins (2nd argument) contains double entries!');
end

F = sparse(n,n);

% Compile list of wanted interactions
Spins = sort(Spins);
[idx1,idx2] = find(tril(ones(numel(Spins)),-1));
idx = [idx2,idx1];

Pairs = Spins(idx);
nPairs = size(Pairs,1);
Coupl = Pairs(:,1) + (Pairs(:,2)-1)*System.nNuclei;

% Compile list of all nuclear spin pairs
[n2,n1] = find(tril(ones(System.nNuclei),-1));
allPairsIdx = n1 + (n2-1)*System.nNuclei;

nn = System.nn;
if ~System.fullnn
  nnFrame = System.nnFrame;
end

% Bilinear coupling term I1*nn*I2
%----------------------------------------------------------------
for iPair = 1:nPairs
  iCoupling = find(Coupl(iPair)==allPairsIdx);
  
  % Construct matrix representing coupling tensor
  if System.fullnn
    J = nn(3*(iCoupling-1)+(1:3),:);
  else
    R_M2nn = erot(nnFrame(iCoupling,:)); % mol frame -> nn frame
    R_nn2M = R_M2nn.';  % nn frame -> mol frame
    J = R_nn2M*diag(nn(iCoupling,:))*R_nn2M.';
  end
  
  % Sum up Hamiltonian terms
  for c1 = 1:3
    so1 = sop(sys,System.nElectrons+Pairs(iPair,1),c1,'sparse');
    for c2 = 1:3
      so2 = sop(sys,System.nElectrons+Pairs(iPair,2),c2,'sparse');
      F = F + so1*J(c1,c2)*so2;
    end
  end
  
end

F = (F+F')/2; % Hermitianise

if ~sparseResult
  F = full(F); % sparse -> full
end
