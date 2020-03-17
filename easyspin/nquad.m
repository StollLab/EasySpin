% nquad  Nuclear quadrupole interaction Hamiltonian
%
%   Hnq = nquad(SpinSystem)
%   Hnq = nquad(SpinSystem,Nuclei)
%   Hnq = nquad(SpinSystem,Nuclei,'sparse')
%
%   Returns the nuclear quadrupole interactions (NQI)
%   Hamiltonian, in MHz, of the spin system given in
%   SpinSystem. If the vector Nuclei is given, only
%   the NQI of the specified nuclei are computed. 1 is
%   the first nucleus, 2 the second, etc.
%
%   If 'sparse' is given, the matrix is returned in sparse format.

function Hnq = nquad(SpinSystem,Nuclei,Options)

if (nargin==0), help(mfilename); return; end

if (nargin<2), Nuclei = []; end
if (nargin<3), Options = ''; end
if ~ischar(Options)
  error('Third input must be a string, ''sparse''.');
end

sparseOutput = strcmp(Options,'sparse');

[Sys,err] = validatespinsys(SpinSystem);
error(err);

nNuclei = Sys.nNuclei;

if isempty(Nuclei)
  Nuclei = 1:nNuclei;
end

if any(Nuclei>nNuclei) || any(Nuclei<1)
  error('Nuclear spin index/indices (2nd argument) out of range!');
end

Hnq = sparse(Sys.nStates,Sys.nStates);
spvc = Sys.Spins;
nElectrons = Sys.nElectrons;

% Nuclear quadrupole interaction
%---------------------------------------------------------
for iNuc = Nuclei
  
  % Construct Q matrix.
  if Sys.fullQ
    Q = Sys.Q(3*(iNuc-1)+(1:3),:);
  else
    Q = diag(Sys.Q(iNuc,:));
  end
  
  if ~any(Q(:))
    continue
  end
  
  % Apply tensor frame transformation if given.
  if any(Sys.QFrame(iNuc,:))
    R_M2Q = erot(Sys.QFrame(iNuc,:)); % mol frame -> Q frame
    R_Q2M = R_M2Q.'; % Q frame -> mol frame
    Q = R_Q2M*Q*R_Q2M.';
  end
  
  % Construct NQI term of spin hamiltonian.
  for k = 1:3
    I{k} = sop(spvc,[nElectrons+iNuc,k],'sparse');
  end
  for k = 1:3
    for q = 1:3
      Hnq = Hnq + I{k}*Q(k,q)*I{q};
    end
  end
  
end % for all nuclear spins specified

Hnq = (Hnq+Hnq')/2; % Hermitianize

if ~sparseOutput
  Hnq = full(Hnq);
end

return
