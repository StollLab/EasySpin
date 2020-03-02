% levels  Energy levels of a spin system 
%
%   En = levels(SpinSystem)
%   En = levels(SpinSystem,Ori,B)
%   En = levels(SpinSystem,phi,theta,B)
%   [En,Ve] = levels(...);
%
%   Calculates energy levels of a spin system.
%
%   Input:
%   - SpinSystem: spin system specification structure
%   - phi,theta: vectors of orientation angles of
%     magnetic field (in radians); assumed zero if not given
%   - Ori: orientations of the field in the molecular frame
%       a) nx2 array of Euler angles (phi,theta), or
%       b) nx3 array of Euler angles (phi,theta,chi), or
%       c) 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz' for special directions
%   - B: array of magnetic field magnitudes (mT), assumed 0 if not given
%
%   Output:
%   - En: array containing all energy eigenvalues (in MHz), sorted,
%     for all possible (phi,theta,B) or (Ori,B) combinations. Depending
%     on the dimensions of Ori, phi,theta and B, out can be up
%     to 4-dimensional. The dimensions are in the order phi,
%     theta (or Ori), field, level number.

function [Energies,Vectors] = levels(varargin)

error(chkmlver); % Error if MATLAB is too old.

switch nargin
  case 0
    help(mfilename);
    return;
  case 1
    SpinSystem = varargin{1};
    phi = 0;
    theta = 0;
    MagnField = 0;
    OriList = false;
  case 2
    error('Third input argument (magnetic field range) is missing.');
  case 3
    OriList = true;
    [SpinSystem,Ori,MagnField] = deal(varargin{:});
    if ischar(Ori)
      n = letter2vec(Ori);
      Ori = vec2ang(n).';
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
    if MagnField>0
      MagnField = [0 MagnField];
      MagnField = linspace(MagnField(1),MagnField(2),nFieldPoints);
    end
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

[F,GxM,GyM,GzM] = sham(Sys);

if OriList

  % Examine Ori array.
  if size(Ori,2)==2
    Ori(:,3) = 0;
  elseif size(Ori,2)==3
    % ok
  else
    error('Ori must be a string (''x'',''y'',''z'',''xy'',''xz'',''yz'',''xyz''), a Nx2 array ([phi theta]) or a Nx3 array [phi theta chi].');
  end
  nOri = size(Ori,1);
  
  % Pre-allocate results array
  Energies = zeros(nOri,numel(MagnField),length(F));
  if computeVectors
    Vectors = zeros(nOri,numel(MagnField),length(F),length(F));
  end
  
  % Loop over all parameter combinations
  zL = ang2vec(Ori(:,1),Ori(:,2)); % z direction in lab frame
  for iOri = 1:nOri
    G = zL(1,iOri)*GxM + zL(2,iOri)*GyM + zL(3,iOri)*GzM;
    for iField = 1:length(MagnField)
      H = F + MagnField(iField)*G;
      if computeVectors
        [V_,E] = eig(H);
        E = diag(E);
        [E,idx] = sort(E);
        V_ = V_(:,idx);
        Vectors(iOri,iField,:,:) = V_;
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
    GxyM = cosphi(iphi)*GxM + sinphi(iphi)*GyM;
    for itheta = 1:length(theta)
      G = sintheta(itheta)*GxyM + costheta(itheta)*GzM;
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
