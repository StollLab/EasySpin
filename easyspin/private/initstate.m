% initstate  Density matrix of the non-equilibrium state in Sys.initState
%
%   [rho,basis] = initstate(Sys)
%   [rho,basis] = initstate(Sys,H)
%
%   Input:
%     Sys    spin system processed by validatespinsys, with non-empty
%            Sys.initState = {state,basis}
%     H      Hamiltonian (MHz) at the magnetic field of interest; needed
%            to convert a state given in the 'eigen' basis
%
%   Output:
%     rho    density matrix, nStates x nStates
%     basis  'uncoupled', or 'eigen' if the state is given in the 'eigen'
%            basis and H is not provided
%
%   A state given for the electron spins only (dimension prod(2*S+1)) is
%   expanded to the full state space, with equal populations of all nuclear
%   sublevels. The trace is preserved.
%
%   Eigenstates ('eigen', 'zerofield') are ordered by increasing energy. If
%   degenerate eigenstates are assigned different populations or coherences
%   between them, the assignment is ambiguous and an error is issued.

function [rho,basis] = initstate(Sys,H)

state = Sys.initState{1};
basis = Sys.initState{2};

nStates = Sys.nStates;
nElStates = prod(2*Sys.S+1);
nRest = nStates/nElStates;

% Density matrix or population vector, expanded to the full state space
[sz1,sz2] = size(state);
if sz1==sz2
  if sz1~=nElStates && sz1~=nStates
    error('The density matrix in Sys.initState must have dimensions of nxn with n = %d or %d.',nElStates,nStates);
  end
  rho = state;
  if sz1~=nStates
    rho = kron(rho,eye(nRest)/nRest);
  end
elseif isvector(state)
  if numel(state)~=nElStates && numel(state)~=nStates
    error('The population vector in Sys.initState must have %d or %d elements.',nElStates,nStates);
  end
  pop = state(:);
  if numel(pop)~=nStates
    pop = kron(pop,ones(nRest,1)/nRest);
  end
  rho = diag(pop);
else
  error('Sys.initState must contain either a population vector or a density matrix.');
end

% Conversion to the uncoupled basis
switch basis
  case 'uncoupled'
    % nothing to do
  case 'zerofield'
    H0 = ham(Sys,zeros(1,3));
    V = sortedeigenstates(H0,rho,'at zero field');
    rho = V*rho*V';
    basis = 'uncoupled';
  case 'eigen'
    if nargin>1 && ~isempty(H)
      if ~isequal(size(H),[nStates nStates])
        error('The Hamiltonian must have dimensions of %dx%d.',nStates,nStates);
      end
      V = sortedeigenstates(H,rho,'at the given field');
      rho = V*rho*V';
      basis = 'uncoupled';
    end
  otherwise
    error('Unsupported basis ''%s'' in Sys.initState.',basis);
end

end


% Eigenstates of H, sorted by energy. Issues an error if rho (given in the
% eigenbasis) is not a multiple of the identity within each set of degenerate
% eigenstates, since then the result depends on the arbitrary choice of
% eigenvectors within that set.
function V = sortedeigenstates(H,rho,where)

H = full(H);
[V,E] = eig((H+H')/2);
[E,idx] = sort(real(diag(E)));
V = V(:,idx);

tolE = 1e-9*max(abs(E));
tolRho = 1e-9*max(abs(rho(:)));
groupStart = [true; diff(E)>tolE];
groupIdx = cumsum(groupStart);
for g = 1:groupIdx(end)
  idx = find(groupIdx==g);
  if numel(idx)<2, continue; end
  block = rho(idx,idx);
  if any(any(abs(block-block(1,1)*eye(numel(idx)))>tolRho))
    error(['Degenerate energy levels %s prevent an unambiguous assignment of the ' ...
      'state in Sys.initState to the eigenstates. Provide the state as a density ' ...
      'matrix in the uncoupled basis instead.'],where);
  end
end

end
