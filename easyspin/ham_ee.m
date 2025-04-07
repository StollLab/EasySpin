% ham_ee  Electron-electron spin interaction Hamiltonian
%
%   [F, dF] = ham_ee(SpinSystem)
%   [F, dF] = ham_ee(SpinSystem,eSpins)
%   [F, dF] = ham_ee(SpinSystem,eSpins,'sparse')
%
%   Returns the electron-electron spin interaction (EEI)
%   Hamiltonian, in MHz and its derivatives.
%
%   Input:
%   - SpinSystem: Spin system structure. EEI
%       parameters are in the ee, eeFrame, and ee2 fields.
%   - eSpins: If given, specifies electron spins
%       for which the EEI should be computed. If
%       absent, all electrons are included.
%   - 'sparse': If given, the matrix is returned in sparse format.

%   Output:
%   - F: Hamiltonian matrix containing the EEI for
%       electron spins specified in eSpins.
%   - dF: Derivative of the Hamiltonian matrix containing the EEI for
%       electron spins specified in eSpins.

function [F,dF] = ham_ee(System,Spins,opt)

if nargin==0, help(mfilename); return; end

if nargin<1 || nargin>3, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>2, error('Too many output arguments.'); end

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

% Special cases: only one spins, ee not given or all zero
F = sparse(n,n);
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

F = sparse(n,n);

% Compile list of wanted interactions
Spins = sort(Spins);
[idx1,idx2] = find(tril(ones(numel(Spins)),-1));
idx = [idx2,idx1];

Pairs = Spins(idx);
nPairs = size(Pairs,1);
Coupl = Pairs(:,1) + (Pairs(:,2)-1)*System.nElectrons;

% Compile list of all spin pairs
[e2,e1] = find(tril(ones(System.nElectrons),-1));
allPairsIdx = e1 + (e2-1)*System.nElectrons;

ee = System.ee;
if isfield(System,'eeFrame')
  eeFrame = System.eeFrame;
else
  eeFrame = [];
end

%for conistency accross methods
nStates=System.nStates;
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
  else
    R_ee2M = eye(3);
  end

  % preparing the derivatives (specific for each electron)
  dFdeex = sparse(nStates,nStates);
  dFdeey = sparse(nStates,nStates);
  dFdeez = sparse(nStates,nStates);

  % preparing the derivatives coefficient
  deexM = R_ee2M(:,1)*R_ee2M(:,1).';  % rotate derivative wrt eex to molecular frame
  deeyM = R_ee2M(:,2)*R_ee2M(:,2).';  % rotate derivative wrt eey to molecular frame
  deezM = R_ee2M(:,3)*R_ee2M(:,3).';  % rotate derivative wrt eez to molecular frame

  % Sum up Hamiltonian terms
  for c1 = 1:3
    if ~useSparseMatrices  % sparse -> full
      so1 = sop(sys,[Pairs(iPair,1),c1]);
    else
      so1 = sop(sys,[Pairs(iPair,1),c1],'sparse');
    end
    for c2 = 1:3
      if ~useSparseMatrices  % sparse -> full
        so2 = sop(sys,[Pairs(iPair,2),c2]);
      else
        so2 = sop(sys,[Pairs(iPair,2),c2],'sparse');
      end
      tempProduct=so1*so2;
      F = F + J(c1,c2)*tempProduct;
      dFdeex = dFdeex + deexM(c1,c2)*tempProduct;
      dFdeey = dFdeey + deeyM(c1,c2)*tempProduct;
      dFdeez = dFdeez + deezM(c1,c2)*tempProduct;
    end
  end

  dF{iPair} = {dFdeex,dFdeey,dFdeez};  % derivatives

  % Isotropic biquadratic exchange coupling term +ee2*(S1.S2)^2
  %-----------------------------------------------------------------
  if System.ee2(iCoupling)==0, continue; end
  F2 = 0;
  for c = 1:3
    F2 = F2 + sop(sys,[Pairs(iPair,1),c;Pairs(iPair,2),c],'sparse');
  end
  if ~useSparseMatrices  % sparse -> full
    F2=full(F2);
  end
  dFdJ=F2^2;
  F = F + System.ee2(iCoupling)*dFdJ;
  dF{iPair} = {dFdeex,dFdeey,dFdeez,dFdJ};  % derivatives
end

F = (F+F')/2; % Hermitianise
end
