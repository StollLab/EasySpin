% diptensor  Calculate dipolar coupling tensor between two spins
%
%   T = diptensor(g1,g2,rvec)
%   T = diptensor(nuc1,nuc2,rvec)
%   T = diptensor(g1,nuc2,rvec)
%   T = diptensor(nuc1,g2,rvec)
%
% Calculates the dipolar coupling tensor T (3x3 matrix, in MHz) between two
% spins using the inter-spin distance vector rvec (in nm). For electrons,
% provide g tensor information in g1 and g2. For nuclei, provide the
% isotope in nuc1 and nuc2.
%
% The tensor T is for the Hamiltonian H = S1*T*S2, where S1 are the spin operators
% of the two spins, in units of hbar.
%
% The tensor T is defined in the same frame as the vector rvec and the g
% tensors g1 and g2, which is typically the molecular frame.
%
% Inputs:
%   g1     g tensor of electron spin 1: isotropic (1 number) or full tensor
%          (3x3 matrix) or principal values plus Euler angles {gpv,gFrame}
%          Euler angles are in radians and defined as in Sys.gFrame
%   g2     same for spin 2
%   nuc1   nuclear isotope: '1H', '14N', etc.
%   nuc2   same for spin 2
%   rvec   3-element vector pointing from spin1 to spin2, in nm
%
% Outputs:
%   T      dipolar coupling tensor (3x3 matrix), in MHz

function T = diptensor(spin1,spin2,rvec)

if nargin==0
  help(mfilename);
  return
end

if nargout>2
  error('Too many output arguments. A maximum of two is possible.');
end

% Input checks
%-------------------------------------------------------------------------------
if nargin<3
  error('Three inputs are required: spin1, spin2, and rvec.');
end
if nargin>3
  error('Only three inputs are allowed: spin1, spin2, and rvec.');
end

if ~isnumeric(rvec) || numel(rvec)~=3
  error('Third input (rvec) must be a 3-element vector.');
end


% Get gyromagnetic ratios
%-------------------------------------------------------------------------------
spins = {spin1,spin2};

for iSpin = 2:-1:1
  input = spins{iSpin};

  if isnumeric(input)
    
    if numel(input)==1
      g = input*eye(3);
    elseif numel(input)==3
      g = diag(input);
    elseif all(size(input)==[3 3])
      g = input;
    else
      error('Input %d must be a 3-element array or a 3x3 array.',iSpin);
    end
    mug{iSpin} = +bmagn*g;
  
  elseif iscell(input)
    
    if numel(input)==2
      gpv = input{1};
      gFrame = input{2};
      R_M2g = erot(gFrame);
      R_g2M = R_M2g.';
      g = R_g2M*diag(gpv)*R_g2M.';
      mug{iSpin} = +bmagn*g;
    end

  else
    
    try
      gn = nucgval(input);
    catch
      error('Could not determine gn value of nuclear isotope ''%s''.',input);
    end
    if numel(gn)~=1
      error('Provide a single nuclear isotope in input %d.',iSpin);
    end
    mug{iSpin} = -nmagn*gn;

  end

end

% Calculate dipolar tensor
%-------------------------------------------------------------------------------
rvec = rvec(:)*1e-9; % nm -> m
r = norm(rvec);
n = rvec(:)/r;
d = 3*(n*n.') - eye(3);
T = -mu0/(4*pi)*r^-3*mug{1}.'*d*mug{2};  % J
T = T/planck/1e6;  % J -> MHz

end
