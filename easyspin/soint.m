% soint  Spin-orbit interaction Hamiltonian 
%
%   H = soint(SpinSystem)
%   H = soint(SpinSystem,Spins)
%   H = soint(SpinSystem,Spins,'sparse')
%
%   Returns the spin orbit interaction (SOI)
%   Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. SOI
%       parameters are the spin-orbit coupling in soc 
%       and the orbital reduction facotr in orf.
%   - Spins: If given, specifies spins
%       for which the SOI should be computed. If
%       absent, all electrons are included.
%   - 'sparse': If given, the matrix is returned in sparse format.

%   Output:
%   - H: Hamiltonian matrix containing the SOI for
%        spins specified in eSpins.

function H = soint(System,Spins,opt)

if (nargin==0), help(mfilename); return; end

if (nargin<1) || (nargin>3), error('Wrong number of input arguments!'); end
if (nargout>1), error('Too many output arguments.'); end

if nargin<3, opt = ''; end
if nargin<2, Spins = []; end
if ~ischar(opt)
  error('Third input must be a string, ''sparse''.');
end
sparseResult = strcmp(opt,'sparse');

[System,err] = validatespinsys(System);
error(err);

shift = System.nElectrons+System.nNuclei;

% Special cases: no OAM, soc not given or all zero
H = sparse(System.nStates,System.nStates);
if (numel(System.Spins)== shift) || ~any(System.soc(:))
  if ~sparseResult
    H = full(H); % sparse -> full
  end
  return;
end


if isempty(Spins), Spins = 1:System.nElectrons; end

% Some error checking on the second input argument
if numel(Spins)<1
  error('Spins (2nd argument) must contain at least 1 values!'); 
end
if any(Spins<1) || any(Spins>System.nElectrons)
  error('Spins (2nd argument) must contain values between 1 and %d!',System.nElectrons);
end
if numel(unique(Spins))~=numel(Spins)
  error('Spins (2nd argument) contains double entries!');
end

% Compile list of wanted interactions
nPairs = length(Spins);

% Compile list of all spin pairs

% Bilinear coupling term S1*ee*S2
%----------------------------------------------------------------
for k = 1:nPairs
  if any(System.soc(Spins(k)))
    % Sum up Hamiltonian terms
    F = sparse(System.nStates,System.nStates);
    for c1 = 1:3
      %spin operator
      so = sop(System,[Spins(k),Spins(k)+shift],c1,'sparse');
      F = F + so;
    end
    for n=1:length(System.soc(Spins(k),:))
      H = H + System.soc(Spins(k),n)* (System.orf(Spins(k))*F)^n;
    end;
  end
  
end


H = (H+H')/2; % Hermitianise
if ~sparseResult
  H = full(H); % sparse -> full
end
