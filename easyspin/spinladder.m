% spinladder    Manifold Hamiltonian of a spin ladder
%
%    CSys = spinladder(Sys)
%    [CSys,En] = spinladder(Sys)
%    ... = spinladder(Sys,Temp)
%    spinladder(...)
%
%    Given a exchange-coupled two-electron-spin system in Sys,
%    this function computes the spin Hamiltonians for the
%    various spin manifolds in the coupled representation,
%    incl. g, A and D tensors, assuming the strong-exchange
%    limit.
%
%    CSys is a cell array that contains the coupled-spin
%         systems sorted by energy.
%    En   contains the center-of-gravity energies as
%         determined by the exchange coupling.
%
%    If no output is requested, spinladder prints some
%    information about the coupled manifolds.
%
%    If a temperature (in kelvin) is given, populations for
%    each manifold are computed and returned in the fields
%    CSys{}.weight.

% TODO:
% - expand to three and more spins
% - validate strong-exchange limit

% Expression taken from Bencini/Gatteschi book.
% Other relevant reference:
% - Scaringe et al, Molec. Phys. 35(3) 701-713 (1978)

function varargout = spinladder(Sys,Temperature)

switch (nargin)
  case 0, help(mfilename); return;
  case 1, Temperature = inf;
  case 2 % everything ok
  otherwise
    error('One or two input arguments expected.');
end

if ~isstruct(Sys)
  error('First input must be a spin system structure.');
end

if numel(Sys.S)~=2
  error('Only systems with two electron spins are supported.');
end

[Sys_,err] = validatespinsys(Sys);
nNuclei = Sys_.nNuclei;
error(err);

if isfield(Sys_,'ee2') && (Sys_.ee2~=0)
  error('Cannot compute spin ladder with Sys.ee2 present.');
end

% Coupling arithmetic for two spins
%-------------------------------------------------------------
Sa = Sys.S(1);
Sb = Sys.S(2);
S = abs(Sa-Sb):(Sa+Sb);

Sa2 = Sa*(Sa+1);
Sb2 = Sb*(Sb+1);
S2 = S.*(S+1);

% Bencini/Gatteschi p.54 Table 3.2
c = (Sa2-Sb2)./S2;
cp = (3*(Sa2-Sb2)^2+S2.*(3*S2-3-2*Sa2-2*Sb2))./((2*S+3).*(2*S-1).*S2);
cm = (4*S2*(Sa2-Sb2)-3*(Sa2-Sb2))./((2*S+3).*(2*S-1).*S2);

c(S==0) = 0;
cp(S==0) = 0;
cm(S==0) = 0;

% Bencini/Gatteschi p.53 Eq. (3.23)
c1 = (1+c)/2;
c2 = (1-c)/2;
d1 = (cp+cm)/2;
d2 = (cp-cm)/2;
d12 = (1-cp)/2;

% Hamiltonian parameters etc.
%-------------------------------------------------------------

% Bencini/Gatteschi p.50 Eq.(3.5)
full_ee = false;
if isfield(Sys,'ee')
  switch numel(Sys.ee)
    case 1
      J = Sys.ee;
    case 3
      J = mean(Sys.ee);
    case 9
      full_ee = true;
      J = mean(diag(Sys.ee));
    otherwise
      error('Sys.ee must have 1 or 3 elements.');
  end
elseif isfield(Sys,'J')
  J = Sys.J;
end
Energies = J/2*(S2-Sa2-Sb2);

% get site tensors/matrices
if isfield(Sys,'g')
  if (numel(Sys.g)==2)
    g1 = Sys.g(1);
    g2 = Sys.g(2);
  else
    g1 = Sys.g(1,:);
    g2 = Sys.g(2,:);
  end
else
  g1 = 2;
  g2 = 2;
end

if nNuclei>0
  if (numel(Sys.A)==4)
    A1 = Sys.A(:,1);
    A2 = Sys.A(:,2);
  else
    A1 = Sys.A(:,1:3);
    A2 = Sys.A(:,4:6);
  end
end

% Zero-field splitting tensors
%--------------------------------------------------
if ~isfield(Sys,'D')
  Sys.D = zeros(2,3);
end

if numel(Sys.D)==2
  % only D given
  D_(1,:) = Sys.D(1)*[-1,-1,+2]/3;
  D_(2,:) = Sys.D(2)*[-1,-1,+2]/3;
  Sys.D = D_;
elseif numel(Sys.D)==4
  % D and E given
  D_(1,:) = [-1/3 -1/3 +2/3]*Sys.D(1,1) + [+1 -1 0]*Sys.D(1,2);
  D_(2,:) = [-1/3 -1/3 +2/3]*Sys.D(2,1) + [+1 -1 0]*Sys.D(2,2);
  Sys.D = D_;
elseif numel(Sys.D)==6
  % D principal values given
end
D1 = Sys.D(1,:);
D2 = Sys.D(2,:);

if full_ee
  D12 = Sys.ee - J*eye(3);
  D1 = diag(D1);
  D2 = diag(D2);
else
  if isfield(Sys,'ee')
    D12 = Sys.ee - J;
  else
    D12 = 0;
  end
end

% Construct spin systems for manifolds
% (based on Bencini/Gatteschi p.53 Eqs. (3.20)-(3.22))
for iS = numel(S):-1:1
  Sys_ = Sys;
  Sys_.S = S(iS);
    
  % g tensors
  Sys_.g = c1(iS)*g1 + c2(iS)*g2;
  
  % Hyperfine couplings
  if nNuclei>0
    Sys_.A = c1(iS)*A1 + c2(iS)*A2;
  end
  
  % Zero-field splittings
  if Sys_.S>1/2
    Sys_.D = d1(iS)*D1 + d2(iS)*D2 + d12(iS)*D12;
  else
    if isfield(Sys_,'D'), Sys_ = rmfield(Sys_,'D'); end
  end
  
  % remove fields containing coupling constants
  if isfield(Sys_,'ee'), Sys_ = rmfield(Sys_,'ee'); end
  if isfield(Sys_,'J'), Sys_ = rmfield(Sys_,'J'); end
  
  CoupledSystems{iS} = Sys_;
end

% Compute fractional populations if temperature is given
if ~isinf(Temperature)
  Populations = exp(-(Energies-min(Energies))*1e6*planck/boltzm/Temperature);
  Populations = Populations.*(2*S+1);
  Populations = Populations/sum(Populations);
  for iS = numel(S):-1:1
    CoupledSystems{iS}.weight = Populations(iS);
  end
end

% Sort by energy
[Energies,idx] = sort(Energies);
CoupledSystems = CoupledSystems(idx);
if ~isinf(Temperature)
  Populations = Populations(idx);
end

switch nargout
  case 0
    fprintf('S1 = %g, S2 = %g (total %d electronic states)\n',Sa,Sb,sum(2*Sys.S+1));
    fprintf('%d manifolds (highest to lowest energy):\n',numel(S));
    for iS = numel(Energies):-1:1
      if Temperature~=inf
        popStr = sprintf('\n     population %0.3e/state, %0.3e total',...
          Populations(iS)/(2*CoupledSystems{iS}.S+1),Populations(iS));
      else
        popStr = '';
      end
      fprintf('  S = %g (%d electronic states)\n     energy %+g GHz%s\n',...
        CoupledSystems{iS}.S,2*CoupledSystems{iS}.S+1,Energies(iS)/1e3,popStr);
    end
  case 1
    varargout = {CoupledSystems};
  case 2
    varargout = {CoupledSystems,Energies};
end

return
