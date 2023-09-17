% rescaledata    Rescaling of one data vector so that it fits a second
%
%   ynew = rescaledata(y,mode)
%   ynew = rescaledata(y,yref,mode)
%   [ynew, scale] = rescaledata(...)
%
% Shifts and rescales the data vector y. If given, yref serves
% as the reference. mode determines how the scale factor is calculated.
% Positive scaling is enforced, i.e. the rescaled data is never inverted.
%
% Inputs:
%   y           data to be rescaled, 1D array
%   yref        reference data (used by 'maxabs' and 'lsq'), 1D array
%   mode:
%     'maxabs'  scales y such that maximum magnitude of ynew is 1 (no shift)
%     'lsq'     least-squares fit of a*y to yref; yref is needed
%     'none'    no scaling
%
% Outputs:
%   ynew      rescaled data vector
%   scale     scaling factor

function varargout = rescaledata(varargin)

switch nargin
  case 0
    help(mfilename);
    return
  case 1
    y = varargin{1};
    yref = [];
    Mode = 'lsq';
  case 2
    y = varargin{1};
    in2 = varargin{2};
    if ischar(in2)
      yref = [];
      Mode = in2;
    else
      yref = in2;
      Mode = 'lsq';
    end
  case 3
    y = varargin{1};
    yref = varargin{2};
    Mode = varargin{3};
  otherwise
    error('Wrong number of input parameters.');
end    

% Make sure y is a vector
if ~isvector(y)
  error('y must be a row or column vector. rescaledata() does not work on 2D arrays.');
end
isRowVector = size(y,1)==1;

% Make sure yref is a vector
if ~isempty(yref) && ~isvector(yref)
  if isnumeric(yref)
    error('yref must be a row or column vector.');
  else
    error('Second input of three must be the reference vector.');
  end
end

N = numel(y);

Modes = {'maxabs','lsq','none'};
refNeeded = [false true false];
equalLengthNeeded = [false true false];

ModeID = find(strcmp(Mode,Modes));
if isempty(ModeID)
  error('Unknown scaling mode ''%s''',Mode);
end
if refNeeded(ModeID)
  if isempty(yref)
    error('For ''%s'' mode, a reference is needed.',Mode);
  end
  if equalLengthNeeded(ModeID)
    equalLength = numel(y)==numel(yref);
    if equalLengthNeeded(ModeID) && ~equalLength
      error('For least-squares rescaling, vectors must have same number of elements.');
    end
  end
end

y = y(:);
yref = yref(:);
yref_notnan = yref(~isnan(yref));
y_notnan = y(~isnan(y));

% Rescale data
switch ModeID
  case 1  % maxabs
    if ~isempty(yref)
      scalefactor = max(abs(yref_notnan))/max(abs(y_notnan));
    else
      scalefactor = 1/max(abs(y_notnan));
    end
    yscaled = scalefactor*y;
  case 2  % lsq
    notnan_both = ~isnan(y) & ~isnan(yref);
    % Rescale reference instead of signal (otherwise rmsd(scale) is wrong
    scalefactor = yref(notnan_both)\y(notnan_both);
    scalefactor = 1/scalefactor;
    yscaled = scalefactor*y;
  case 3  % no scaling
    yscaled = y;
    scalefactor = 1;
end

% Make sure signal is not inverted
if real(scalefactor(1))<0
  scalefactor(1) = abs(scalefactor(1));
  yscaled = y*scalefactor;
end

% Preserve row layout
if isRowVector
  yscaled = yscaled.';
end

switch nargout
  case 0
    x = 1:N;
    subplot(2,1,1);
    plot(x,y);
    legend('original');
    axis tight
    subplot(2,1,2);
    if isempty(yref)
      plot(x,yscaled);
      legend(sprintf('scaled by %g',scalefactor));
    else
      subplot(2,1,2);
      plot(x,yref,x,yscaled);
      legend('reference',sprintf('scaled by %g',scalefactor));
    end
    axis tight
    varargout = {};
  case 1
    varargout = {yscaled};
  case 2
    varargout = {yscaled, scalefactor};
  otherwise
    error('Wrong number of outputs.');
end

end
