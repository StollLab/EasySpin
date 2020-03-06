% diptensor  Calculate dipolar coupling tensor between two spins
%
%   T = diptensor(spin1,spin2,rvec)
%
% Calculates the dipolar coupling tensor T (3x3 matrix, in MHz) between spin1
% and spin2, using the inter-spin distance vector rvec (in nm).
%
% The tensor T is for the Hamiltonian H = J1*T*J2, where J1 is the spin vector
% operator of spin1, and J2 is the spin vector operator for spin2.
%
% Inputs:
%   spin1  type of first spin: 'e', '1H', '14N', etc.
%   spin2  type of second spin: 'e', '1H', '14N', etc.
%   rvec   3-element vector pointing from spin1 to spin2, in nm
%
% Outputs:
%   T      dipolar coupling tensor (3x3 matrix), in MHz
%
% For electrons, g = gfree = 2.002319... is used.

function T = diptensor(spin1,spin2,rvec)

% Input checks
%-------------------------------------------------------------------------------
if nargin<3
  error('Three inputs are required: spin1, spin2, and rvec.');
end
if nargin>3
  error('Only three inputs are allowed: spin1, spin2, and rvec.');
end

if ~ischar(spin1)
  error('First input (spin1) must be a character array (''e'', ''1H'', etc.).');
end
if ~ischar(spin2)
  error('Second input (spin2) must be a character array (''e'', ''1H'', etc.).');
end
if ~isnumeric(rvec) || numel(rvec)~=3
  error('Third input (rvec) must be a 3-element vector.');
end


% Get gyromagnetic ratios
%-------------------------------------------------------------------------------
spins = {spin1,spin2};

for iSpin = 1:2
  if strcmp(spins{iSpin},'e')
    mug(iSpin) = +bmagn*gfree;
  else
    try
      gn = nucgval(spins{iSpin});
    catch
      error('Could not determine gn value of isotope ''%s''.',spins{iSpin});
    end
    if numel(gn)~=1
      error('Provide a single isotope in spin%d.',iSpin);
    end
    mug(iSpin) = -nmagn*gn;
  end
end


% Calculate dipolar tensor
%-------------------------------------------------------------------------------
rvec = rvec(:)*1e-9; % nm -> m
r = norm(rvec);
n = rvec(:)/r;
T = -(mu0/4/pi)*r^-3*mug(1)*mug(2)*(3*(n*n')-eye(3)); % J
T = T/planck/1e6; % J -> MHz

return

