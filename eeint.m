% eeint  Electron-electron spin interaction Hamiltonian 
%
%   F = eeint(SpinSystem)
%   F = eeint(SpinSystem,eSpins)
%
%   Returns the electron-electron spin interaction (EEI)
%   Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. EEI
%       parameters are in the ee, eepa, and ee2 fields.
%   - eSpins: If given, specifies electron spins
%       for which the EEI should be computed. If
%       absent, all electrons are included.
%
%   Output:
%   - F: Hamiltonian matrix containing the EEI for
%       electron spins specified in eSpins.

function F = eeint(System,Spins)

if (nargin==0), help(mfilename); return; end

if (nargin<1) || (nargin>2), error('Wrong number of input arguments!'); end
if (nargout<0), error('Not enough output arguments.'); end
if (nargout>1), error('Too many output arguments.'); end

[System,err] = validatespinsys(System);
error(err);

sys = spinvec(System);
n = prod(2*sys+1);

% Special cases: only one spins, ee not given or all zero
F = zeros(n,n);
if (System.nElectrons==1), return; end
if ~any(System.ee(:)) && ~any(System.ee2(:)), return; end

if (nargin<2), Spins = 1:System.nElectrons; end

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
[idx1,idx2] = find(triu(ones(numel(Spins)),1));
idx = [idx1,idx2];

Pairs = Spins(idx);
nPairs = size(Pairs,1);
Coupl = Pairs(:,1) + (Pairs(:,2)-1)*System.nElectrons;

% Compile list of all spin pairs
[e2,e1] = find(tril(ones(System.nElectrons),-1));
allPairsIdx = e1 + (e2-1)*System.nElectrons;

ee = System.ee;
if ~System.fullee
  eepa = System.eepa;
end

% Isotropic biquadratic exchange
ee2 = System.ee2;

% Compute Hamiltonian matrix
for iPair = 1:nPairs
  iCoupling = find(Coupl(iPair)==allPairsIdx);

  % Construct matrix representing coupling tensor
  if System.fullee
    J = ee(3*(iCoupling-1)+(1:3),:);
  else
    Rp = erot(eepa(iCoupling,:));
    J = Rp*diag(ee(iCoupling,:))*Rp.';
  end
  
  % Sum up Hamiltonian terms
  for c1 = 1:3
    so1 = sop(sys,Pairs(iPair,1),c1,'sparse');
    for c2 = 1:3
      so2 = sop(sys,Pairs(iPair,2),c2,'sparse');
      
      % Bilinear coupling S1*ee*S2
      F = F + so1*J(c1,c2)*so2;

      % Isotropic biquadratic exchange coupling +ee2*(S1.S2)^2
      if (c1==c2)
        if ee2(iPair)~=0
          F = F + ee2(iPair)*(so1*so2)^2;
        end
      end

    end
  end
  
end

F = full(F); % sparse -> full
F = (F+F')/2; % Hermitianise
