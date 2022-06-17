% photoselect Calculate orientation-dependent photo-selection weight
%
%   weight = photoselect(tdm,labFrame,k,alpha)
%
% This function calculates a selection weight between 0 and 1 for photo-excited
% spin centers, 1 indicating full excitation and 0 indicating no
% excitation. It needs information about the electric transition dipole moment
% orientation (tdm), the spin center orientation (labFrame), the light
% excitation direction (k) and polarization (alpha).
%
% Inputs:
%   tdm         Orientation of the transition dipole moment vector in the
%               molecular frame. There are three ways to input this:
%               - a letter or letter combination, e.g. 'x', 'z', 'xz', '-y', etc.
%               - a three-element vector [mx my mz]
%               - two spherical angles [phim thetam] (in radians)
%   labFrame    [phi theta chi] angles of lab frame in molecular frame;
%               phi and theta determine the direction of the lab z axis (zL, aligned
%               with B0), chi additionally determines the direction of lab x and y
%               (xL and yL)
%   k           Orientation of the propagation direction of the light excitation beam
%               in the lab frame
%               There are three ways to input this:
%               - a letter or letter combination, e.g. 'x', 'z', 'xz', '-y', etc.
%               - a three-element vector [kx ky kz]
%               - two spherical angles [phik thetak] (in radians)
%               The most common situation is k = 'y' (perpendicular to
%               B0) or k = [pi/2 pi/2] or k = [0; 1; 0]
%   alpha       Polarization angle (in radians) that indicates the orientation
%               of the E-field vector in the plane perpendicular to the propagation
%               direction k. Use alpha=NaN to indicate depolarized light.
%               For k = 'y', alpha = 0 puts the E-field along the lab z
%               axis (zL), parallel to B0; alpha = pi/2 puts the E-field along
%               the lab x axis (perpendicular to zL and B0)
%
% Outputs:
%   weight      photoselection weight, between 0 and 1
%               0 means no excitation at all, 1 indicates complete excitation
%
% Example:
%
%   k = 'y';  % propagation direction along lab y axis
%   photoselect(tdm,labFrame,k,pi/2)  % E-field along lab x axis
%   photoselect(tdm,labFrame,k,0)     % E-field along lab z axis

function weight = photoselect(tdm,labFrame,k,alpha)

% Parse tdm orientation, calculate tdm unit vector in lab frame representation
if ischar(tdm)
  tdm_mol = letter2vec(tdm);
elseif numel(tdm)==3
  tdm_mol = tdm(:)/norm(tdm(:));
elseif numel(tdm)==2
  tdm_mol = ang2vec(tdm(1),tdm(2));
else
  error('`tdmOri (first input) must be a letter designating a direction, a 3-vector, or an array with two angles.');
end

% Parse k orientation
if ischar(k)
  kOri = vec2ang(letter2vec(k));
elseif numel(k)==3
  kOri = vec2ang(k);
elseif numel(k)==2
  kOri = k;
else
  error('k (third input) must be a letter designating a direction, a 3-vector, or an array with two angles.');
end

% Parse polarization angle
if ~isfloat(alpha) || ~isscalar(alpha)
  error('Polarization angle alpha (4th input) must be a scalar.')
end
unpolarized = isnan(alpha);
if unpolarized
  alpha = 0; % to prevent NaN-issues with erot
end

% Calculate unit vectors of light frame (p,q,k)
% k = propagation direction
% p = E-field direction
lightFrame = [kOri(1) kOri(2) alpha];
[p_lab,~,k_lab] = erot(lightFrame,'rows');

% Calculate tdm unit vector in lab frame representation
R_M2L = erot(labFrame);
tdm_lab = R_M2L*tdm_mol; % mol frame -> lab frame

% Calculate photoselection weight
if unpolarized
  weight = (1 - abs(tdm_lab'*k_lab)^2)/2;
  %weight = (abs(tdm_lab'*p_lab)^2 + abs(tdm_lab'*q_lab)^2)/2; % equivalent
else
  weight = abs(tdm_lab'*p_lab)^2;
end

end
