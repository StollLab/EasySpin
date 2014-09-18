% Determine excitation mode
%-----------------------------------------------------------------
% Exp.Polarization:
%   'linear' (default)
%   'circular' ( = 'circular+')
%   'circular+'
%   'circular-'
%   'unpolarized'
% Exp.Mode:
%   'perpendicular'
%   'parallel'
%   [k_tilt]
%   [k_tilt alpha_pol]

% B0 defines z(Lab) axis, B0 and k define z(Lab)y(Lab) plane, y(Lab) is such
%  that k points towards positive y(Lab)
% k_tilt: angle between B0 and k (wave vector), between 0 and pi
% alpha_pol: polarization angle: angle of B1 off the x(Lab) axis,
%   in the plane perpendicular to k; for alpha_pol = 0, B1 is along
%   x(Lab) axis.

% Polarization mode
%--------------------------------------------------------------------
if ~isfield(Exp,'Polarization'), Exp.Polarization = ''; end

unpolarizedMode = false;
circpolarizedMode = false;
linearpolarizedMode = false;
circSense = NaN;
if ~isempty(Exp.Polarization)
  if ~ischar(Exp.Polarization)
    error('Exp.Polarization must be a string.');
  end
  switch Exp.Polarization
    case 'linear', linearpolarizedMode = true; logstr = 'linear';
    case 'circular', circpolarizedMode = true; circSense = +1; logstr = 'circular+';
    case 'circular+', circpolarizedMode = true; circSense = +1; logstr = 'circular+';
    case 'circular-', circpolarizedMode = true; circSense = -1; logstr = 'circular-';
    case 'unpolarized', unpolarizedMode = true; logstr = 'unpolarized';
    otherwise
      error('Unrecognized Exp.Polarization: ''%s''.',Exp.Polarization);
  end
else
  linearpolarizedMode = true;
  logstr = 'linear';
end
logmsg(1,'  polarization mode: %s',logstr);

% Excitation mode
%--------------------------------------------------------------------
if ~isfield(Exp,'Mode'), Exp.Mode = ''; end

if isempty(Exp.Mode)
  if linearpolarizedMode
    Exp.Mode = [pi/2 pi/2]; % perpendicular
  else
    Exp.Mode = [0 0];
  end
end

ParallelMode = false;
k_tilt = 0;
alpha_pol = 0;
if ischar(Exp.Mode)
  switch Exp.Mode
    case 'perpendicular'
      k_tilt = pi/2; alpha_pol = pi/2; logstr = 'perpendicular';
    case 'parallel'
      k_tilt = pi/2; alpha_pol = 0; logstr = 'parallel'; ParallelMode = true;
    otherwise
      error('Unrecognized Exp.Mode: ''%s''.',Exp.Mode);
  end
  logmsg(1,'  resonator mode: %s',logstr);
else
  switch numel(Exp.Mode)
    case 0
    case 1
      k_tilt = Exp.Mode(1);
    case 2
      k_tilt = Exp.Mode(1);
      alpha_pol = Exp.Mode(2);
    otherwise
      error('Exp.Mode must contain 1 or 2 angles (k tilt angle, polarization angle');
  end
end
logmsg(1,'  mode angles: k tilt %g deg, polarization angle %g deg',k_tilt*180/pi,alpha_pol*180/pi);

if (k_tilt<0) || (k_tilt>pi)
  error('  The k tilt angle in Exp.Mode must be between 0 and pi.');
end

% Set up vectors for transition rate calculations
%---------------------------------------------------------------------
% all vectors are in lab frame coordinates
nB0 = [0; 0; 1]; % unit vector along B0, in lab coordinates

% Transformation from (nx,ny,nB0) lab frame to (nB1,*,nk) frame
R = erot(-pi/2,-k_tilt,alpha_pol);
nB1 = R(1,:).'; % unit vector along B1, in lab coordinates
nk = R(3,:).'; % unit vector along wave vector k, in lab coordinates

xi1 = nB1.'*nB0;
xik = nk.'*nB0;
