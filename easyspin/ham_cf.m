% ham_cf  Crystal-field Hamiltonian for the orbital angular momentum
%
%   H = ham_cf(Sys)
%   H = ham_cf(Sys,OAMs)
%   H = ham_cf(Sys,OAMs,'sparse')
%
%   Returns the crystal-field Hamiltonian (in MHz) of the spin system
%   Sys, utilizing the fields Sys.L, Sys.CF1, Sys.CF2, etc.
%
%   If the vector OAMs is given, the crystal field of only the
%   specified orbital angular momenta is returned (1 is the first in Sys.L,
%   2 the second, etc). Otherwise, all orbital angular momenta are included.
%
%   If 'sparse' is given, the matrix is returned in sparse format.

function H = ham_cf(SpinSystem,idxL,opt)

if nargin==0, help(mfilename); return; end

if nargin<2, idxL = []; end
if nargin<3, opt = ''; end

if ~ischar(opt)
  error('Third input must be a string, ''sparse''.');
end
useSparseMatrices = strcmp(opt,'sparse');

[Sys,err] = validatespinsys(SpinSystem);
error(err);

if isempty(idxL)
  idxL = 1:Sys.nL;
end
if any(idxL>Sys.nL) || any(idxL<1)
  error('OAM index/indices (2nd argument) out of range!');
end

H = sparse(Sys.nStates,Sys.nStates);

% if no orbital angular momentum is defined, H remains all zero  
if Sys.nL==0
  if ~useSparseMatrices
    H = full(H);
  end
  return
end

offset = Sys.nElectrons+Sys.nNuclei;  % skip over electrons and nuclei in stev()
L = Sys.L;
for iL = idxL
  
  % L < 1: no crystal-field interaction possible
  if L(iL)<1, continue; end

  % High-order terms in extended Stevens operator format
  %---------------------------------------------------------
  % CF1, CF2, CF3, ... (k = 1...12)
  
  % Run over all ranks k (Ck = C1, C2, C3, C4, ...)
  for k = 1:12
    fieldname = sprintf('CF%d',k);
    
    if ~isfield(Sys,fieldname), continue; end
    
    CFk = Sys.(fieldname);
    if isempty(CFk), continue; end
    
    if numel(CFk)~=(2*k+1)*Sys.nL
      error('Sys.%s has wrong size. It should contain %d elements for each L.',fieldname,2*k+1);
    end
    if all(CFk==0), continue; end
    
    q = k:-1:-k;
    for iq = find(CFk(iL,:)~=0)
      H = H + CFk(iL,iq)*stev(Sys.Spins,[k,q(iq),iL+offset]);
    end
    
  end % for all tensor ranks

end % for all electron spins specified

H = (H+H')/2; % Hermitianise

if ~useSparseMatrices
  H = full(H);
end

end
