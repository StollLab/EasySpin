% nnint  Nuclear-nuclear spin interaction Hamiltonian 
%
%   Hnn = nnint(SpinSystem)
%   Hnn = nnint(SpinSystem,nucSpins)
%   Hnn = nnint(SpinSystem,nucSpins,'sparse')
%
%   Returns the nuclear-nuclear spin interaction (NNI)
%   Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. NNI
%       parameters are in the nn and nnFrame fields.
%   - nucSpins: Indices of nuclear spins for which the NNI should be
%       computed. E.g. [1 3] indicates the first and third nuclear spin.
%       If nucspins is omitted or empty, all nuclei are included.
%   - 'sparse': If given, the matrix is returned in sparse format.
%
%   Output:
%   - Hnn: Hamiltonian matrix containing the NNI for nuclear spins specified
%       in nucSpins.

function Hnn = nnint(System,nucSpins,opt)

if nargin==0, help(mfilename); return; end

if nargin<1 || nargin>3, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>1, error('Too many output arguments.'); end

if nargin<3, opt = ''; end
if nargin<2, nucSpins = []; end
if ~ischar(opt)
  error('Third input must be a string, ''sparse''.');
end
sparseResult = strcmp(opt,'sparse');

[System,err] = validatespinsys(System);
error(err);

sys = spinvec(System);
n = prod(2*sys+1);

Hnn = sparse(n,n);

% Special cases: only one nuclear spin, nn not given, or all zero
if System.nNuclei<2 || ~any(System.nn(:))
  if ~sparseResult
    Hnn = full(Hnn);
  end
  return
end

if isempty(nucSpins), nucSpins = 1:System.nNuclei; end

% Some error checking on the second input argument
if numel(nucSpins)<2
  error('Spins (2nd argument) must contain at least 2 values!'); 
end
if any(nucSpins<1) || any(nucSpins>System.nNuclei)
  error('Spins (2nd argument) must contain values between 1 and %d!',System.nNuclei);
end
if numel(unique(nucSpins))~=numel(nucSpins)
  error('Spins (2nd argument) contains double entries!');
end

Hnn = sparse(n,n);

% Compile list of wanted interactions
nucSpins = sort(nucSpins);
[idx1,idx2] = find(tril(ones(numel(nucSpins)),-1));
idx = [idx2,idx1];

nucPairs = nucSpins(idx);
nNucPairs = size(nucPairs,1);
Coupl = nucPairs(:,1) + (nucPairs(:,2)-1)*System.nNuclei;

% Compile list of all nuclear spin pairs
[n2,n1] = find(tril(ones(System.nNuclei),-1));
allPairsIdx = n1 + (n2-1)*System.nNuclei;

nn = System.nn;
if ~System.fullnn
  nnFrame = System.nnFrame;
end

% Bilinear coupling term I1*nn*I2
%-------------------------------------------------------------------------------
for iNucPair = 1:nNucPairs
  iCoupling = find(Coupl(iNucPair)==allPairsIdx);
  
  % Construct matrix representing coupling tensor
  if System.fullnn
    J = nn(3*(iCoupling-1)+(1:3),:);
  else
    R_M2nn = erot(nnFrame(iCoupling,:)); % mol frame -> nn frame
    R_nn2M = R_M2nn.';  % nn frame -> mol frame
    J = R_nn2M*diag(nn(iCoupling,:))*R_nn2M.';
  end
  
  % Sum up Hamiltonian terms
  spin1 = System.nElectrons+nucPairs(iNucPair,1);
  spin2 = System.nElectrons+nucPairs(iNucPair,2);
  for c1 = 1:3
    so1 = sop(sys,[spin1,c1],'sparse');
    for c2 = 1:3
      so2 = sop(sys,[spin2,c2],'sparse');
      Hnn = Hnn + so1*J(c1,c2)*so2;
    end
  end
  
end

Hnn = (Hnn+Hnn')/2; % Hermitianize

if ~sparseResult
  Hnn = full(Hnn); % sparse -> full
end
