% nquad  Nuclear quadrupole interaction Hamiltonian 
%
%   F = nquad(SpinSystem)
%   F = nquad(SpinSystem,Nuclei)
%
%   Returns the nuclear quadrupole interactions
%   (NQI) Hamiltonian in MHz of the system
%   SpinSystem. If the vector Nuclei is given,
%   only the NQI of the specified nuclei are
%   computed. 1 is the first nucleus, etc.

function F = nquad(SpinSystem,Nuclei)

if (nargin==0), help(mfilename); return; end

[Sys,err] = validatespinsys(SpinSystem);
error(err);

nNuclei = Sys.nNuclei;

if (nargin==1)
  Nuclei = 1:nNuclei;
end

if any(Nuclei>nNuclei) | any(Nuclei<1)
  error('Nuclear spin index/indices (2nd argument) out of range!');
end

F = sparse(Sys.nStates,Sys.nStates);
nucspvc = Sys.I;
spvc = Sys.Spins;
nElectrons = Sys.nElectrons;

for iNuc = 1:length(Nuclei)
  
  idx = Nuclei(iNuc);

  % I < 1 -> no internal interactions possible -> go to next spin
  if nucspvc(idx)<1, continue; end
  
  % Nuclear quadrupole interaction
  %---------------------------------------------------------
  if any(Sys.Q(idx,:))
    if ~Sys.fullQ
      % Construct Q matrix.
      Rp = erot(Sys.Qpa(idx,:));
      Q = Rp*diag(Sys.Q(idx,:))*Rp.';
    else
      Q = Sys.Q(3*(idx-1)+(1:3),:);
    end
    Q = Q*Sys.Qscale(idx);
    % Construct NQI term.
    for k = 1:3
      so = sop(spvc,idx+nElectrons,k,'sparse');
      for q = 1:3
        F = F + so*Q(k,q)*sop(spvc,idx+nElectrons,q,'sparse');
      end
    end
  end
    
end % for all spins specified

F = full(F);
F = (F+F')/2; % Hermitianise

return
