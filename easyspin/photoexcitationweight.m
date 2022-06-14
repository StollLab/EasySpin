% photoexcitationweight Calculate photoexcitation weight
%
% weight = photoexcitationweight(tdmAngles,poldir,labFrame,pol)
%
% Inputs:
%   tdmAngles       [phi theta] angles of transition dipole moment vector
%                     in molecular frame
%   poldir          laser E-field direction in lab frame, specified as 
%                     a letter ('x' for perpendicular, 'z' for parallel, etc)
%                     or two Euler angles [phi theta]
%                     or one Euler angle [theta] (which sets phi=0)
%   labFrame        [phi theta chi] angles of lab frame in molecular frame;
%                     phi and theta determine the direction of B0 (lab z)
%                     chi additionally determines the direction of lab x and y
%   pol             degree of polarization, between 0 and 1
%
% Outputs:
%   weight          photoexcitation weight, between 0 and 1
%                   0 means no, 1 indicates maximal photoexcitation

function weight = photoexcitationweight(tdmAngles,poldir,labFrame,pol)

if isempty(poldir)
  weight = 1;
  return
end

% Assume full polarization if pol is not given
if nargin<4
  pol = 1;
end

% Calculate tdm unit vector in lab frame representation
tdm_mol = ang2vec(tdmAngles(1),tdmAngles(2));
R_M2L = erot(labFrame);
tdm_lab = R_M2L*tdm_mol; % mol frame -> lab frame

% Determine angle between laser E-field and B0 field
% Calculate E-field direction in lab frame representation
if ischar(poldir)
  Edir_lab = letter2vec(poldir);
elseif numel(poldir)==1
  alpha = 0;
  beta = poldir;
  Edir_lab = ang2vec([alpha beta]);
else
  alpha = poldir(1);
  beta = poldir(2);
  Edir_lab = ang2vec([alpha beta]);
end

% Calculate weight for polarized excitation
weight_pol = abs(tdm_lab'*Edir_lab)^2;

% Calculate total weight, including unpolarized contribution
weight_unpol = 1; % weight for completely unpolarized light
weight = pol*weight_pol + (1-pol)*weight_unpol;

end
