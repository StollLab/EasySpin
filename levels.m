% levels  Energy levels of a spin system 
%
%   En = levels(SpinSystem,Ori,B)
%   En = levels(SpinSystem,phi,theta,B)
%   [En,Ve] = levels(...);
%
%   Calculates energy levels of a spin system.
%
%   Input:
%   - SpinSystem: spin system specification structure
%   - phi,theta: vectors of orientation angles of
%     magnetic field [radians]
%   - Ori: 2xn or nx2 array of orientations (phi,theta).
%       or 'x', 'y', 'z'.
%   - B: vector of magnetic field magnitudes [mT]
%
%   Output:
%   - En: array containing all energy eigenvalues [MHz], sorted,
%     for all possible (phi,theta,B) or (Ori,B) combinations. Depending
%     on the dimensions of Ori, phi,theta and B, out can be up
%     to 4-dimensional. The dimensions are in the order phi,
%     theta (or Ori), field, level number.

function [Energies,Vectors] = levels(varargin)

error(chkmlver); % Error if Matlab too old.

switch nargin
case 0
  help(mfilename); return;
case 2
  error('Third input argument (magnetic field range) is missing.');
case 3
  OriList = 1;
  [SpinSystem,Ori,MagnField] = deal(varargin{:});
  if ischar(Ori)
    if numel(Ori)~=1
      error('Only ''x'', ''y'', or ''z'' are allowed.')
    end
    switch Ori
      case 'x', Ori = [0 pi/2];
      case 'y', Ori = [pi/2 pi/2];
      case 'z', Ori = [0 0];
      otherwise
        error('Only ''x'', ''y'', or ''z'' are allowed.')
    end
  end
case 4
  OriList = 0;
  [SpinSystem,phi,theta,MagnField] = deal(varargin{:});
otherwise
  error('Wrong number of input arguments!')
end

N = 200;
switch numel(MagnField)
  case 1
    MagnField = [0 MagnField];
    MagnField = linspace(MagnField(1),MagnField(2),N);
  case 2
    MagnField = linspace(MagnField(1),MagnField(2),N);
  otherwise
end


switch nargout
case 0, ComputeVectors = 0;
case 1, ComputeVectors = 0;
case 2, ComputeVectors = 1;
otherwise, error('Wrong number of output arguments');
end

% Check spin system and pre-calculate spin Hamiltonian components.
[Sys,err] = validatespinsys(SpinSystem);
error(err);

[F,Gx,Gy,Gz] = sham(Sys);

if OriList

  % Examine Ori array.
  sOri = size(Ori);
  if sOri(1)==2
    if sOri(2)==2
      warning('Orientations array ambiguous, taking orientations along columns.');
    end
    nOri = sOri(2);
  elseif sOri(2)==2
    nOri = sOri(1);
    Ori = Ori.';
  else
    error('Orientation array must be 2xn or nx2.');
  end
  
  % Pre-allocate results array
  Energies = zeros(nOri,numel(MagnField),length(F));
  if ComputeVectors
    Vectors = zeros(nOri,numel(MagnField),length(F),length(F));
  end
  
  % Loop over all parameter combinations
  v = ang2vec(Ori(1,:),Ori(2,:));
  for iOri = 1:nOri
    G = v(1,iOri)*Gx + v(2,iOri)*Gy + v(3,iOri)*Gz;
    for iField = 1:length(MagnField)
      H = F + MagnField(iField)*G;
      if ComputeVectors
        [Vectors(iOri,iField,:,:),E] = eig(H);
        E = diag(E);
      else
        E = eig(H);
        E = sort(E);
      end
      Energies(iOri,iField,:) = E;
    end
  end

else

  % Pre-calculate trigonometric functions
  cosphi = cos(phi);
  sinphi = sin(phi);
  costheta = cos(theta);
  sintheta = sin(theta);
  
  % Pre-allocate results array
  Energies = zeros(numel(phi),numel(theta),numel(MagnField),length(F));
  if ComputeVectors
    Vectors = zeros(numel(phi),numel(theta),numel(MagnField),length(F),length(F));
  end
  
  % Loop over all parameter combinations
  for iphi = 1:length(phi)
    Gplane = cosphi(iphi)*Gx + sinphi(iphi)*Gy;
    for itheta = 1:length(theta)
      G = sintheta(itheta)*Gplane + costheta(itheta)*Gz;
      for iField = 1:length(MagnField)
        if ComputeVectors
          [Vectors(iphi,itheta,iField,:,:),E] = eig(F + MagnField(iField)*G);
          E = diag(E);
        else
          E = sort(eig(F + MagnField(iField)*G));
        end
        Energies(iphi,itheta,iField,:) = E;
      end
    end
  end
end

% Remove singleton dimensions.
Energies = squeeze(Energies);
if ComputeVectors, Vectors = squeeze(Vectors); end

return
