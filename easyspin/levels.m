% levels  Energy levels of a spin system 
%
%   E = levels(SpinSystem,Ori,B)
%   E = levels(SpinSystem,phi,theta,B)
%   [E,V] = levels(...)
%
%   Calculates energy levels of a spin system.
%
%   Input:
%   - SpinSystem: spin system structure
%   - phi,theta: vectors of polar angles (in radians) specifying orientation
%                of magnetic field in molecular frame
%   - Ori: orientations of the field in the molecular frame
%       a) nx2 array of Euler angles (phi,theta), or
%       b) nx3 array of Euler angles (phi,theta,chi), or
%       c) 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz' for special directions
%       phi and theta are the polar angles; the field direction is
%       independent of chi
%   - B: array of magnetic field magnitudes (mT)
%
%   Output:
%   - E: array containing all energy eigenvalues (in MHz), sorted,
%        for all possible (phi,theta,B) or (Ori,B) combinations. Depending
%        on the dimensions of Ori, phi, theta and B, E can be up
%        to 4-dimensional. The dimensions are in the order phi,
%        theta (or Ori), field, level number.
%   - V: array of eigenvectors
%
%   Example:
%
%    Sys = struct('S',5/2,'D',1000);
%    B = linspace(0,600,601);  % mT
%    E = levels(Sys,'z',B);
%    plot(B,E);

function [Energies,Vectors] = levels(varargin)

warning(chkmlver);  % Error if MATLAB is too old.

switch nargin
  case 0
    help(mfilename);
    return
  case 1
    error('At least three input arguments (Sys, Ori, B0) expected.');
  case 2
    error('Third input argument (magnetic field range) is missing.');
  case 3
    OriList = true;
    [SpinSystem,Ori,MagnField] = deal(varargin{:});
    if ischar(Ori)
      zL = letter2vec(Ori);
      Ori = vec2ang(zL).';
    end
  case 4
    OriList = false;
    [SpinSystem,phi,theta,MagnField] = deal(varargin{:});
  otherwise
    error('Wrong number of input arguments!')
end

nFieldPoints = 101;
switch numel(MagnField)
  case 1
    % single field value given
  case 2
    MagnField = linspace(MagnField(1),MagnField(2),nFieldPoints);
  otherwise
    % vector of magnetic field values given
end

switch nargout
case 0, computeVectors = false;
case 1, computeVectors = false;
case 2, computeVectors = true;
otherwise, error('Wrong number of output arguments!');
end

% Check spin system
[Sys,err] = validatespinsys(SpinSystem);
error(err);

% Pre-calculate spin Hamiltonian components
[H0,muxM,muyM,muzM] = ham(Sys);

if OriList

  % Examine Ori array
  if size(Ori,2)==2
    Ori(:,3) = 0;
  elseif size(Ori,2)==3
    % ok
  else
    error('Ori must be a string (''x'',''y'',''z'',''xy'',''xz'',''yz'',''xyz''), a Nx2 array ([phi theta]) or a Nx3 array [phi theta chi].');
  end
  nOri = size(Ori,1);
  
  % Pre-allocate results array
  Energies = zeros(nOri,numel(MagnField),length(H0));
  if computeVectors
    Vectors = zeros(nOri,numel(MagnField),length(H0),length(H0));
  end
  
  % Loop over all parameter combinations
  zL_M = ang2vec(Ori(:,1),Ori(:,2));  % lab frame z direction in molecular frame
  for iOri = 1:nOri
    muzL = zL_M(1,iOri)*muxM + zL_M(2,iOri)*muyM + zL_M(3,iOri)*muzM;
    for iField = 1:length(MagnField)
      H = H0 - MagnField(iField)*muzL;
      if computeVectors
        [V_,E_] = eig(H);
        E_ = diag(E_);
        Vectors(iOri,iField,:,:) = V_;
      else
        E_ = eig(H);
        E_ = sort(E_);
      end
      Energies(iOri,iField,:) = E_;
    end
  end

else

  % Pre-calculate trigonometric functions
  cosphi = cos(phi);
  sinphi = sin(phi);
  costheta = cos(theta);
  sintheta = sin(theta);
  
  % Pre-allocate output arrays
  Energies = zeros(numel(phi),numel(theta),numel(MagnField),length(H0));
  if computeVectors
    Vectors = zeros(numel(phi),numel(theta),numel(MagnField),length(H0),length(H0));
  end
  
  % Loop over all parameter combinations
  for iphi = 1:length(phi)
    muxyM = cosphi(iphi)*muxM + sinphi(iphi)*muyM;
    for itheta = 1:length(theta)
      muzL = sintheta(itheta)*muxyM + costheta(itheta)*muzM;
      for iField = 1:length(MagnField)
        H = H0 - MagnField(iField)*muzL;
        if computeVectors
          [V_,E_] = eig(H);
          E_ = diag(E_);
          Vectors(iphi,itheta,iField,:,:) = V_;
        else
          E_ = eig(H);
          E_ = sort(E_);
        end
        Energies(iphi,itheta,iField,:) = E_;
      end
    end
  end
  
end

% Remove singleton dimensions
Energies = squeeze(Energies);
if computeVectors
  Vectors = squeeze(Vectors);
end

end
