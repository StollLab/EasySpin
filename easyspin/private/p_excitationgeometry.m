function [xi1,xik,nB1,nk,nB0,mode] = p_excitationgeometry(mwMode)

% Parse microwave excitation mode
%===============================================================================
% Exp.mwMode:
%   'perpendicular' (default)
%   'parallel'
%   {k_tilt alpha_pol}  % linear
%   {k_tilt 'circular+'}
%   {k_tilt 'circular-'}
%   {k_tilt 'unpolarized'}

% B0 defines z(Lab) axis, B0 and k define z(Lab)y(Lab) plane, y(Lab) is such
%  that k points into the half-plane with y(Lab)>0
% k_tilt: angle between B0 and k (wave vector), between 0 and pi
% alpha_pol: polarization angle: angle of B1 off the x(Lab) axis,
%   in the plane perpendicular to k; for alpha_pol = 0, B1 is along
%   x(Lab) axis.

% Microwave excitation mode
%-------------------------------------------------------------------------------
if isempty(mwMode)
  mwMode = 'perpendicular';
end

unpolarizedMode = false;
circpolarizedMode = false;
linearpolarizedMode = false;
circSense = NaN;
parallelMode = false;

if ischar(mwMode)
  linearpolarizedMode = true;
  switch mwMode
    case 'perpendicular'
      k_tilt = pi/2; alpha_pol = pi/2; logstr = 'perpendicular';
    case 'parallel'
      k_tilt = pi/2; alpha_pol = 0; logstr = 'parallel';
      parallelMode = true;
    otherwise
      error('Unrecognized Exp.mwMode: ''%s''.',mwMode);
  end
  logmsg(1,'  resonator mode: %s',logstr);
elseif iscell(mwMode) && numel(mwMode)==2
  k_tilt = mwMode{1};
  Mode2 = mwMode{2};
  if isnumeric(Mode2)
    linearpolarizedMode = true;
    alpha_pol = Mode2;
  else
    switch Mode2
      case 'unpolarized', unpolarizedMode = true;
      case 'circular+', circpolarizedMode = true; circSense = +1;
      case 'circular-', circpolarizedMode = true; circSense = -1;
      otherwise
        error('Unrecognized 2nd element in Exp.mwMode: ''%s''',Mode2);
    end
    alpha_pol = 0;
  end
else
  error('Exp.mwMode must be ''perpendicular'' or ''parallel'' or a 2-element cell array.');
end

logmsg(1,'  mode angles: k tilt %g deg, polarization angle %g deg',k_tilt*180/pi,alpha_pol*180/pi);

if k_tilt<0 || k_tilt>pi
  error('  The k tilt angle in Exp.mwMode must be between 0 and pi.');
end

% Set up vectors for transition rate calculations
%-------------------------------------------------------------------------------
% all vectors are in lab frame coordinates
nB0 = [0; 0; 1]; % unit vector along B0, in lab coordinates

% Transformation from (nx,ny,nB0) lab frame to (nB1,*,nk) frame
% (1) using erot()
%{
R = erot(-pi/2,-k_tilt,alpha_pol);
nk = R(3,:).'; % unit vector along wave vector k, in lab coordinates
nB1 = R(1,:).'; % unit vector along B1, in lab coordinates
%}
% (2) manual
sk = sin(k_tilt);
ck = cos(k_tilt);
nk = [0; sk; ck]; % unit vector along wave vector k, in lab coordinates
ca = cos(alpha_pol);
sa = sin(alpha_pol);
nB1 = [sa; -ca*ck; ca*sk]; % unit vector along mw field vector B1, in lab coordinates

% projections of k and B1 unit vectors onto B0 direction
xik = nk.'*nB0;
xi1 = nB1.'*nB0;

mode.unpolarizedMode = unpolarizedMode;
mode.circpolarizedMode = circpolarizedMode;
mode.linearpolarizedMode = linearpolarizedMode;
mode.circSense = circSense;
mode.parallelMode = parallelMode;
