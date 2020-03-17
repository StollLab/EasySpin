% crystalfield  Crystal-field Hamiltonian for the orbital angular momentu,
%
%   F = zfield(SpinSystem)
%   F = zfield(SpinSystem,OAMs)
%   F = zfield(SpinSystem,OAMs,'sparse')
%
%   Returns the crystal-field
%   Hamiltonian [MHz] of the system SpinSystem.
%
%   If the vector OAMs is given, the Crystal-field of only the
%   specified orbital angular momentums is returned (1 is the first, 2 the
%   second, etc). Otherwise, all orbital angular momentums are included.
%
%   If 'sparse' is given, the matrix is returned in sparse format.

function H = crystalfield(SpinSystem,idxOAM,opt)

if (nargin==0), help(mfilename); return; end

[Sys,err] = validatespinsys(SpinSystem);
error(err);

if (nargin<2), idxOAM = []; end
if (nargin<3), opt = ''; end
if ~ischar(opt)
  error('Third input must be a string, ''sparse''.');
end
sparseResult = strcmp(opt,'sparse');

if isempty(idxOAM)
  idxOAM = 1:Sys.nElectrons;
end

if any(idxOAM>Sys.nElectrons) || any(idxOAM<1),
  error('OAM index/indices (2nd argument) out of range!');
end

H = sparse(Sys.nStates,Sys.nStates);
% index of OAM starts with 1, but in Sys.Spins and in H, they appear after 
% electron spin and nuclei
offset = Sys.nElectrons+Sys.nNuclei;
% if no orbital angular momentum is defined, H remains all zero  
if numel(Sys.Spins)==offset
  if ~sparseResult
    H = full(H);
  end
  return;
end
OAM = Sys.Spins(1+offset:end);

for iOAM = idxOAM
  
  % L < 1: no zero-field interaction possible
  if OAM(iOAM)<1, continue; end

  % High-order terms in extended Stevens operator format
  %---------------------------------------------------------
  % CF1, CF2, CF3, ... (k = 1...12)
  
  % Run over all ranks k (Ck = C1, C2, C3, C4, ...)
  for k = 1:12
    fieldname = sprintf('CF%d',k);
    
    if ~isfield(Sys,fieldname), continue; end
    
    CFk = Sys.(fieldname);
    if isempty(CFk), continue; end
    
    if numel(CFk)~=Sys.nElectrons && numel(CFk)~=(2*k+1)*Sys.nElectrons
      error('Sys.%s has wrong size. It should contain 1 or %d elements',fieldname,2*k+1);
    end
    if all(CFk==0), continue; end
    
    q = k:-1:-k;
    for iq = find(CFk(iOAM,:)~=0)
      H = H + CFk(iOAM,iq)*stev(Sys.Spins,[k,q(iq),iOAM+offset]);
    end
    
  end % for all tensor ranks

end % for all electron spins specified

H = (H+H')/2; % Hermitianise
if ~sparseResult
  H = full(H);
end

return
