% rescale     Rescaling of one data vector so that it fits a second
%
%   ynew = rescale(y,mode1)
%   ynew = rescale(y,yref,mode2)
%   [ynew, scalefactors] = rescale(...)
%
%   Shifts and rescales the data vector y. If given, ynew serves
%   as the reference. The rescaled y is returned in ynew.
%
%   y and yref need to be 1D vectors.
%
% Inputs:
%   mode1:
%     'minmax'  shifts and scales to minimum 0 and maximum 1
%     'maxabs'  scales to maximum abs 1, no shift
%     'none'    no scaling
%
%   mode2:
%     'minmax'  shifts&scales so that minimum and maximum of y
%               fit the minimum and maximum of yref
%     'maxabs'  scales y so that maximum absolute values of y fits yref
%     'lsq'     least-squares fit of a*y to yref
%     'lsq0'    least-squares fit of a*y+b to yref
%     'lsq1'    least-squares fit of a*y+b+c*x to yref
%     'lsq2'    least-squares fit of a*y+b+c*x+d*x^2 to yref
%     'lsq3'    least-squares fit of a*y+b+c*x+d*x^2+e*x^3 to yref
%     'none'    no scaling
%
%   Positive scaling is enforced, i.e. no inverting 
%   of the rescaled data  
%
% Output:
%   ynew          the new rescaled y vector
%   scalefactors  the scaling factors and polynomial coefficients

function varargout = rescale(varargin)

switch (nargin)
  case 0
    help(mfilename);
    return
  case 1
    error('A second input argument (scaling mode) is needed!');
  case 2
    y = varargin{1};
    in2 = varargin{2};
    if ischar(in2)
      yref = [];
      Mode = in2;
    else
      error('A third input argument (scaling mode) is needed!');
    end
  case 3
    y = varargin{1};
    yref = varargin{2};
    Mode = varargin{3};
  otherwise
    error('Wrong number of input parameters.');
end    

% Make sure y is a vector
N = numel(y);
if (length(y)~=N)
  error('y must be a row or column vector. rescale() does not work on 2D arrays.');
end
isRowVector = size(y,1)==1;

% Make sure yref is a vector
if (length(yref)~=numel(yref))
  error('yref must be a row or column vector. rescale() does not work on 2D arrays.');
end

if isempty(yref)
  
  % Rescaling without reference
  %----------------------------------------------------

  ModeID = find(strcmp(Mode,{'minmax','maxabs','none'}));
  if isempty(ModeID)
    error('Unknown scaling mode ''%s''',Mode);
  end
  
  y = y(:);
  
  switch ModeID
    case 1 % minmax
      mi = 0;
      ma = 1;
      ynew = mi + (ma-mi)/(max(y)-min(y))*(y-min(y));
    case 2 % maxabs
      ynew = y/max(abs(y));
    case 3 % no scaling
      ynew = y;
  end
  
else
  
  % Rescaling with reference
  %----------------------------------------------------
  ModeID = find(strcmp(Mode,{'maxabs','minmax','shift','lsq','lsq0','lsq1','lsq2','lsq3','none'}));
  equalLengthNeeded = [0 0 0 1 1 1 1 1 0];
  if isempty(ModeID)
    error('Unknown scaling mode ''%s''',Mode);
  end
  equalLength = numel(y)==numel(yref);
  if equalLengthNeeded(ModeID) && ~equalLength
    error('For least-squares rescaling, vectors must have same number of elements.');
  end
  
  y = y(:);
  yref = yref(:);
  yref_notnan = yref(~isnan(yref));
  y_notnan = y(~isnan(y));
  if equalLength
    notnan_both = ~isnan(y) & ~isnan(yref);
  end

  switch ModeID
    case 1 % maxabs
      scalefactor = max(abs(yref_notnan ))/max(abs(y_notnan));
      ynew = scalefactor*y;
    case 2 % minmax
      scalefactor = (max(yref_notnan)-min(yref_notnan))/(max(y_notnan)-min(y_notnan));
      ynew = scalefactor*(y-min(y_notnan)) + min(yref_notnan);
    case 3 % shift
      shift = mean(y_notnan) - mean(yref_notnan);
      ynew = y - shift;
      scalefactor = [1 shift];
    case 4 % lsq
      scalefactor = y(notnan_both)\yref(notnan_both);
      ynew = scalefactor*y;
    case 5 % lsq0
      D = [y ones(N,1)];
      scalefactor = D(notnan_both,:)\yref(notnan_both);
      ynew = D*scalefactor;
    case 6 % lsq1
      x = (1:N).'/N;
      D = [y ones(N,1) x];
      scalefactor = D(notnan_both,:)\yref(notnan_both);
      ynew = D*scalefactor;
    case 7 % lsq2
      x = (1:N).'/N;
      D = [y ones(N,1) x x.^2];
      scalefactor = D(notnan_both,:)\yref(notnan_both);
      ynew = D*scalefactor;
    case 8 % lsq3
      x = (1:N).'/N;
      D = [y ones(N,1) x x.^2 x.^3];
      scalefactor = D(notnan_both,:)\yref(notnan_both);
      ynew = D*scalefactor;
    case 9 % no scaling
      ynew = y;
  end
  if real(scalefactor(1))<0
    scalefactor(1) = abs(scalefactor(1));
    ynew = D*scalefactor;
  end
end

% Preserve row layout
if isRowVector
  ynew = ynew.';
end

switch nargout
  case 0
    x = 1:N;
    if isempty(yref)
      subplot(2,1,1);
      plot(x,y,'r'); legend('original');
      subplot(2,1,2);
      plot(x,ynew,'g'); legend('result');
    else
      plot(x,y,'r',x,yref,'k',x,ynew,'g');
      legend('original','reference','result');
    end
    varargout = {};
  case 1
    varargout = {ynew};
  case 2
    varargout = {ynew scalefactor};
  otherwise
    error('Wrong number of outputs.');
end

return
