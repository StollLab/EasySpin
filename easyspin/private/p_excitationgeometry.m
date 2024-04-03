function [xi1,xik,nB1,nk,nB0,mode] = p_excitationgeometry(mwMode)

% Parse microwave excitation mode
%===============================================================================
% Exp.mwMode:
%   'perpendicular' (default)
%   'parallel'
%   {k alpha}  % linear
%   {k 'circular+'}
%   {k 'circular-'}
%   {k 'unpolarized'}

% B0 defines z(Lab) axis, B0 and k define z(Lab)y(Lab) plane, y(Lab) is such
%  that k points into the half-plane with y(Lab)>0
% k: letter, polar angles, or vector specifying microwave propagation direction
% alpha: polarization angle: angle of B1 off the x(Lab) axis,
%   in the plane perpendicular to k; for alpha_pol = 0, B1 is along
%   x(Lab) axis.

% Parse mwMode
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
      k = 'y'; alpha = -pi/2;  % gives B1 along xL
      logstr = 'perpendicular';
    case 'parallel'
      k = 'y'; alpha = pi;  % gives B1 along zL
      logstr = 'parallel';
      parallelMode = true;
    otherwise
      error('Unrecognized Exp.mwMode: ''%s''.',mwMode);
  end
elseif iscell(mwMode) && numel(mwMode)==2
  k = mwMode{1};
  Mode2 = mwMode{2};
  if isnumeric(Mode2)
    linearpolarizedMode = true;
    alpha = Mode2;
  else
    switch Mode2
      case 'unpolarized', unpolarizedMode = true;
      case 'circular+', circpolarizedMode = true; circSense = +1;
      case 'circular-', circpolarizedMode = true; circSense = -1;
      otherwise
        error('Unrecognized 2nd element in Exp.mwMode: ''%s''',Mode2);
    end
    alpha = 0;
  end
  logstr = 'user-specified';
else
  error('Exp.mwMode must be ''perpendicular'' or ''parallel'' or a 2-element cell array.');
end

logmsg(1,'  microwave excitation mode: %s',logstr);

% Parse k
%-------------------------------------------------------------------------------
% Parse k input, calculate polar angles of k vector, kOri: [phi theta]
if ischar(k)
  kOri = vec2ang(letter2vec(k));
elseif numel(k)==3
  kOri = vec2ang(k);
elseif numel(k)==2
  kOri = k;
elseif numel(k)==1
  kOri = [0 k];
else
  error('k (first element in mwMode) must be a letter designating a direction, a 3-vector, or an array with two angles.');
end

% Set up vectors for transition rate calculations
%-------------------------------------------------------------------------------
% all vectors are in lab frame coordinates

% Transformation from (nx,ny,nB0) lab frame to (nB1,*,nk) microwave frame
% nB1 = microwave B-field direction
% nk = microwave propagation direction
mwFrame = [kOri(1) kOri(2) alpha];
[nB1,~,nk] = erot(mwFrame,'rows');

% projections of k and B1 unit vectors onto B0 direction
nB0 = [0; 0; 1]; % unit vector along B0, in lab coordinates
xik = nk.'*nB0;
xi1 = nB1.'*nB0;

mode.unpolarizedMode = unpolarizedMode;
mode.circpolarizedMode = circpolarizedMode;
mode.linearpolarizedMode = linearpolarizedMode;
mode.circSense = circSense;
mode.parallelMode = parallelMode;

