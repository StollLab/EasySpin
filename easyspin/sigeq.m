% sigeq  Thermal equilibrium spin density 
%
%   sigma = sigeq(Ham, Temp)
%   sigma = sigeq(Ham, Temp, 'pol')
%
%   Calculates the equilibrium density matrix
%   at temperature Temp of the system described by
%   the Hamiltonian Ham.
%
%   Input:
%   - Ham: Hamiltonian of the spin system (in MHz)
%   - Temp: temperature of the system (in K)
%   - 'pol': if present, only the polarization part
%       will be returned.
%
%   Output:
%   - sigma: density matrix at thermal equilibrium
%       in the same basis as Ham

function sigma = sigeq(Ham,Temp,varargin)

if (nargin==0), help(mfilename); return; end

if (nargin<2) || (nargin>3), error('Wrong number of input arguments!'); end

if ~isreal(Temp) | (numel(Temp)~=1) | (Temp<=0)
  error('Temperature (in K) must be a positive scalar!');
end

% Calculate thermal equilibrium density matrix.
sigma = expm(-1e6*planck/boltzm/Temp*Ham);
sigma = sigma/trace(sigma);
sigma = (sigma+sigma')/2;% Hermitianize

% If demanded, return polarization part only.
if (nargin==3)
  if strcmpi(varargin{1},'pol')
    N = length(sigma);
    sigma = sigma - eye(N)/N;
  else
    error('Invalid third argument!')
  end
end

return
