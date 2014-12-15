% hfine  Hyperfine interaction Hamiltonian 
%
%   F = hfine(System)
%   F = hfine(System,Spins)
%   F = hfine(System,Spins,'sparse')
%
%   Returns the hyperfine interaction Hamiltonian
%   [MHz] of the spins 'Spins' of the system
%   'SpinSystem'. Spins = 1 is the first electron.
%
%   If 'sparse' is given, the matrix is returned in sparse format.

function F = hfine(System,Spins,opt)

if (nargin==0), help(mfilename); return; end

if (nargin<2), Spins = []; end
if (nargin<3), opt = ''; end
if ~ischar(opt)
  error('Third input must be a string, ''sparse''.');
end
sparseResult = strcmp(opt,'sparse');

% Validate spin system.
[Sys,err] = validatespinsys(System);
error(err);

% Get spin vector and space dimension
SpinVec = Sys.Spins;
nStates = Sys.nStates;
nElectrons = Sys.nElectrons;
nNuclei = Sys.nNuclei;

F = sparse(nStates,nStates); % sparse zero matrix

if (nNuclei==0)
  if ~sparseResult
    F = full(F); % convert sparse to full
  end
  return
end

% Get spin list
if isempty(Spins) % no 'Spins' specified -> all nuclei and all electrons
  NucSpins = nElectrons+1:numel(SpinVec);
  ElSpins = 1:nElectrons;
else
  if any(Spins<1) || any(Spins>numel(SpinVec))
    error('Spins (2nd argument) contains out-of-range values!');
  end
  NucSpins = Spins;
  ElSpins = Spins;
  NucSpins(NucSpins<=nElectrons) = []; % remove electron spins
  ElSpins(ElSpins>nElectrons) = []; % remove nuclear spins
end

if (numel(NucSpins)<1) || (numel(ElSpins)<1)
  error('At least one electron and one nuclear spin must be specified!');
end


FullAMatrix = size(System.A,1)>nNuclei;

% Generate Hamiltonian for hyperfine interaction.
for edx = 1:length(ElSpins)
  eSp = ElSpins(edx);
  idx = (eSp-1)*3+(1:3);
  for ndx = 1:length(NucSpins)
    nSp = NucSpins(ndx);
    if Sys.I(nSp-nElectrons)<=0, continue; end
    % Construct full hyperfine matrix.
    if FullAMatrix
      A = Sys.A((nSp-nElectrons-1)*3+(1:3),idx);
    else
      A = diag(Sys.A(nSp-nElectrons,idx));
    end
    R_M2A = erot(Sys.AFrame(nSp-nElectrons,idx)); % mol frame -> A frame
    R_A2M = R_A2M.'; % A frame -> mol frame
    A = R_A2M*A*R_A2M.';
    A = A * Sys.Ascale(ndx);
    % Construct hyperfine Hamiltonian.
    for k = 1:3
      for q = 1:3
        F = F + A(k,q)*sop(SpinVec,[eSp,nSp],[k,q],'sparse');
      end
    end
  end % for all specified nuclei
end % for all specified electrons

F = (F+F')/2; % hermitianise, e.g. guards against small imaginary remainders on the diagonal
if ~sparseResult
  F = full(F); % convert sparse to full
end

return
