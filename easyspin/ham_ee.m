% ham_ee  Electron-electron spin interaction Hamiltonian 
%
%   Hee = ham_ee(SpinSystem)
%   Hee = ham_ee(SpinSystem,eSpins)
%   Hee = ham_ee(SpinSystem,eSpins,'sparse')
%
%   Returns the electron-electron spin interaction (EEI)
%   Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. EEI
%       parameters are in the ee, eeFrame, and ee2 fields.
%   - eSpins: If given, specifies electron spins
%       for which the EEI should be computed. If
%       absent, all electrons are included.
%   - 'sparse': If given, the matrix is returned in sparse format.

%   Output:
%   - Hee: Hamiltonian matrix containing the EEI for
%       electron spins specified in eSpins.

function Hee = ham_ee(System,Spins,opt)

if nargin==0, help(mfilename); return; end

if nargin<1 || nargin>3, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>1, error('Too many output arguments.'); end

if nargin<3, opt = ''; end
if nargin<2, Spins = []; end
if ~ischar(opt)
  error('Third input must be a string, ''sparse''.');
end
useSparseMatrices = strcmp(opt,'sparse');

[System,err] = validatespinsys(System);
error(err);

sys = spinvec(System);
n = prod(2*sys+1);

Hee = sparse(n,n);

% Special cases: only one spins, ee not given or all zero
if (System.nElectrons==1), return; end
if ~any(System.ee(:)) && ~any(System.ee2(:)), return; end

if isempty(Spins), Spins = 1:System.nElectrons; end

% Some error checking on the second input argument
if numel(Spins)<2
  error('Spins (2nd argument) must contain at least 2 values!'); 
end
if any(Spins<1) || any(Spins>System.nElectrons)
  error('Spins (2nd argument) must contain values between 1 and %d!',System.nElectrons);
end
if numel(unique(Spins))~=numel(Spins)
  error('Spins (2nd argument) contains double entries!');
end

% Compile list of wanted pairs
idx = nchoosek(1:numel(Spins),2);
Spins = sort(Spins);
Pairs = Spins(idx);
nPairs = size(Pairs,1);
Coupl = Pairs(:,1) + (Pairs(:,2)-1)*System.nElectrons;

% Compile list of all spin pairs
idx_all = nchoosek(1:System.nElectrons,2);
allPairsIdx = idx_all(:,1) + (idx_all(:,2)-1)*System.nElectrons;

ee = System.ee;
if isfield(System,'eeFrame')
  eeFrame = System.eeFrame;
else
  eeFrame = [];
end

% Bilinear coupling term S1*ee*S2
%----------------------------------------------------------------
for iPair = 1:nPairs
  iCoupling = find(Coupl(iPair)==allPairsIdx);
  
  % Construct matrix representing coupling tensor
  if System.fullee
    J = ee(3*(iCoupling-1)+(1:3),:);
  else
    J = diag(ee(iCoupling,:));
  end
  if ~isempty(eeFrame) && any(eeFrame(iCoupling,:))
    R_M2ee = erot(eeFrame(iCoupling,:)); % mol frame -> ee frame
    R_ee2M = R_M2ee.';  % ee frame -> mol frame
    J = R_ee2M*J*R_ee2M.';
  end
  
  % Sum up Hamiltonian terms
  for c1 = 1:3
    so1 = sop(sys,[Pairs(iPair,1),c1],'sparse');
    for c2 = 1:3
      so2 = sop(sys,[Pairs(iPair,2),c2],'sparse');
      Hee = Hee + so1*J(c1,c2)*so2;
    end
  end
  
  % Isotropic biquadratic exchange coupling term +ee2*(S1.S2)^2
  %-----------------------------------------------------------------  
  if System.ee2(iCoupling)==0, continue; end
  Hee2 = 0;
  for c = 1:3
    Hee2 = Hee2 + sop(sys,[Pairs(iPair,1),c;Pairs(iPair,2),c],'sparse');
  end
  Hee = Hee + System.ee2(iCoupling)*Hee2^2;
end

Hee = (Hee+Hee')/2; % Hermitianise
if ~useSparseMatrices
  Hee = full(Hee); % sparse -> full
end

end
