% ham_hf  Hyperfine interaction Hamiltonian
%
%   Hhf = ham_hf(System)
%   Hhf = ham_hf(System,eSpins)
%   Hhf = ham_hf(System,eSpins,nSpins)
%   Hhf = ham_hf(System,eSpins,nSpins,'sparse')
%
%   Returns the hyperfine interaction Hamiltonian (in units of MHz) between
%   electron spins 'eSpins' and nuclear spins 'nSpins' of the system
%   'System'. eSpins=1 is the first electron spins, nSpins=1 is the first
%   nuclear spin.
%
%   If 'sparse' is given, the matrix is returned in sparse format.

function [Hhf,dHhf] = ham_hf(System,elSpins,nucSpins,opt)

if nargin==0
  help(mfilename);
  return
end

switch nargin
  case 1
    elSpins = [];
    nucSpins = [];
    opt = '';
  case 2
    nucSpins = [];
    opt = '';
  case 3
    opt = '';
  case 4
  otherwise
    error('Incorrect number of input arguments.')
end

if ~ischar(opt)
  error('Fourth input must be a string, ''sparse''.');
end
useSparseMatrices = strcmp(opt,'sparse');

% Validate spin system.
[Sys,err] = validatespinsys(System);
error(err);

% Get spin vector and space dimension
SpinVec = Sys.Spins;
nStates = Sys.nStates;
nElectrons = Sys.nElectrons;
nNuclei = Sys.nNuclei;

Hhf = sparse(nStates,nStates); % sparse zero matrix

% Special case: no nuclei present, so no hyperfine
if nNuclei==0
  if ~useSparseMatrices
    Hhf = full(Hhf);
  end
  return
end

% Get electron spin list
if isempty(elSpins)
  elSpins = 1:nElectrons;
else
  if any(elSpins<1) || any(elSpins>nElectrons)
    error('Electron spins (2nd argument) contains out-of-range values!');
  end
end

% Get electron spin list
if isempty(nucSpins)
  nucSpins = 1:nNuclei;
else
  if any(nucSpins<1) || any(nucSpins>nNuclei)
    error('Nuclear spins (3rd argument) contains out-of-range values!');
  end
end

fullAMatrix = size(System.A,1)>nNuclei;

% Generate Hamiltonian for hyperfine interaction
for eSp = elSpins
  eidx = (eSp-1)*3+(1:3);
  for nSp = nucSpins
    if Sys.I(nSp)==0, continue; end

    % Construct full hyperfine matrix
    if fullAMatrix
      A = Sys.A((nSp-1)*3+(1:3),eidx);
    else
      A = diag(Sys.A(nSp,eidx));
    end

    if ~any(A(:))
      continue
    end

    % Transform matrix into molecular frame representation
    ang = Sys.AFrame(nSp,eidx);
    if any(ang)
      R_M2A = erot(ang);  % mol frame -> A frame
      R_A2M = R_M2A.';    % A frame -> mol frame
      A = R_A2M*A*R_A2M.';
    else
      R_A2M=eye(3);
    end

    % preparing the derivatives (specific for each electron)
    dHhfx = sparse(nStates,nStates);
    dHhfy = sparse(nStates,nStates);
    dHhfz = sparse(nStates,nStates);

    % preparing the derivatives coefficient
    dAxM = R_A2M(:,1)*R_A2M(:,1).';  % rotate derivative wrt Ax to molecular frame
    dAyM = R_A2M(:,2)*R_A2M(:,2).';  % rotate derivative wrt Ay to molecular frame
    dAzM = R_A2M(:,3)*R_A2M(:,3).';  % rotate derivative wrt Az to molecular frame

    % Construct hyperfine Hamiltonian
    for c1 = 1:3
      for c2 = 1:3
        tempProduct=sop(SpinVec,[eSp c1; nElectrons+nSp c2],'sparse');
        Hhf = Hhf + A(c1,c2)*tempProduct;
        dHhfx = dHhfx + dAxM(c1,c2)*tempProduct;
        dHhfy = dHhfy + dAyM(c1,c2)*tempProduct;
        dHhfz = dHhfz + dAzM(c1,c2)*tempProduct;
      end
    end

  end % for all specified nuclei
  dHhf{eSp,nSp}={dHhfx,dHhfy,dHhfz}; % --> the indices in the derivative may not be needed?
end % for all specified electrons

Hhf = (Hhf+Hhf')/2; % hermitianise, e.g. guards against small imaginary remainders on the diagonal
if ~useSparseMatrices
  Hhf = full(Hhf); % convert sparse to full
end

end
