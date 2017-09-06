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
%     magnetic field (in radians)
%   - Ori: orientations of the field in the molecular frame
%       a) nx2 array of Euler angles (phi,theta), or
%       b) nx3 array of Euler angles (phi,theta,chi), or
%       c) 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz' for
%          special directions
%   - B: vector of magnetic field magnitudes (mT)
%
%   Output:
%   - En: array containing all energy eigenvalues (in MHz), sorted,
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
  OriList = true;
  [SpinSystem,Ori,MagnField] = deal(varargin{:});
  if ischar(Ori)
    switch Ori
      case 'x', Ori = [0, pi/2];
      case 'y', Ori = [pi/2, pi/2];
      case 'z', Ori = [0, 0];
      case 'xy', Ori = [pi/4 pi/2];
      case 'xz', Ori = [0 pi/4];
      case 'yz', Ori = [pi/2 pi/4];
      case 'xyz', Ori = [pi/4 acos(1/sqrt(3))];
      otherwise
        error('Unknown value ''%s'' for orientation (2nd input argument).',Ori);
    end
  end
case 4
  OriList = false;
  [SpinSystem,phi,theta,MagnField] = deal(varargin{:});
otherwise
  error('Wrong number of input arguments!')
end

nFieldPoints = 200;
switch numel(MagnField)
  case 1
    MagnField = [0 MagnField];
    MagnField = linspace(MagnField(1),MagnField(2),nFieldPoints);
  case 2
    MagnField = linspace(MagnField(1),MagnField(2),nFieldPoints);
  otherwise
    % vector of magnetic field values given
end


switch nargout
case 0, computeVectors = false;
case 1, computeVectors = false;
case 2, computeVectors = true;
otherwise, error('Wrong number of output arguments');
end

% Check spin system and pre-calculate spin Hamiltonian components.
[Sys,err] = validatespinsys(SpinSystem);
error(err);

[F,Gx,Gy,Gz] = sham(Sys);

if OriList

  % Examine Ori array.
  if size(Ori,2)==2
    Ori(:,3) = 0;
  elseif size(Ori,2)==3
    % ok
  else
    error('Ori must be a string (''x'',''y'',''z'',''xy'',''xz'',''yz'',''xyz''), a two-element array ([phi theta]) or a three-element array [phi theta chi].');
  end
  nOri = size(Ori,1);
  
  % Pre-allocate results array
  Energies = zeros(nOri,numel(MagnField),length(F));
  if computeVectors
    Vectors = zeros(nOri,numel(MagnField),length(F),length(F));
  end
  
  % Loop over all parameter combinations
  v = ang2vec(Ori(:,1),Ori(:,2));
  for iOri = 1:nOri
    G = v(1,iOri)*Gx + v(2,iOri)*Gy + v(3,iOri)*Gz;
    for iField = 1:length(MagnField)
      H = F + MagnField(iField)*G;
      if computeVectors
        [V_,E] = eig(H);
        E = diag(E);
        [E,idx] = sort(E);
        Vectors(iOri,iField,:,:) = V_(:,idx);
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
  if computeVectors
    Vectors = zeros(numel(phi),numel(theta),numel(MagnField),length(F),length(F));
  end
  
  % Loop over all parameter combinations
  for iphi = 1:length(phi)
    Gplane = cosphi(iphi)*Gx + sinphi(iphi)*Gy;
    for itheta = 1:length(theta)
      G = sintheta(itheta)*Gplane + costheta(itheta)*Gz;
      for iField = 1:length(MagnField)
        if computeVectors
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
if computeVectors
  Vectors = squeeze(Vectors);
end

return
