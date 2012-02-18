% rescale     Rescaling of one spectrum so that it fits a second
%
%   ynew = rescale(y,mode1)
%   ynew = rescale(y,yref,mode2)
%
%   Shifts and rescales the spectrum y. If given, ynew serves
%   as the reference. The rescaled y is returned in ynew.
%
%   mode1:
%     'minmax'  shifts&scales to minimum 0 and maximum 1
%     'maxabs'  scales to maximum abs 1, no shift
%
%   mode2:
%     'minmax'  shifts&scales so that minimum and maximum fit yref
%     'maxabs'  scales so that maximum abs fits yref
%     'lsq'     least-squares fit of a*y
%     'lsq0'    least-squares fit of a*y+b
%     'lsq1'    least-squares fit of a*y+b+c*x
%     'lsq2'    least-squares fit of a*y+b+c*x+d*x^2

function varargout = rescale(varargin)

switch (nargin)
  case 0
    help(mfilename);
    return
  case 1
    y = varargin{1};
    yref = [];
    Mode = 'maxabs';
  case 2
    y = varargin{1};
    in2 = varargin{2};
    if ischar(in2)
      yref = [];
      Mode = in2;
    else
      yref = in2;
      Mode = 'maxabs';
    end
  case 3
    y = varargin{1};
    yref = varargin{2};
    Mode = varargin{3};
  otherwise
    error('Wrong number of input parameters.');
end    

N = numel(y);
if (length(y)~=N)
  error('y must be a vector.');
end
isRowVector = size(y,1)==1;

if isempty(yref)
  
  % Rescaling without reference
  %----------------------------------------------------

  ModeID = find(strcmp(Mode,{'minmax','maxabs'}));
  if isempty(ModeID)
    error('Unknown scaling mode ''%s''',Mode);
  end
  
  y = y(:);
  
  switch ModeID
    case 1 % minmax
      mi = 0; ma = 1;
      ynew = mi + (ma-mi)/(max(y)-min(y))*(y-min(y));
    case 2 % maxabs
      ynew = y/max(abs(y));
  end
  
else
  
  % Rescaling with reference
  %----------------------------------------------------
  
  ModeID = find(strcmp(Mode,{'maxabs','minmax','shift','lsq','lsq0','lsq1','lsq2'}));
  if isempty(ModeID)
    error('Unknown scaling mode ''%s''',Mode);
  end
  if (ModeID>=4)
    if numel(y)~=numel(yref)
      error('For least-squares rescaling, vectors must have same number of elements.');
    end
  end
  
  y = y(:);
  yref = yref(:);

  switch ModeID
    case 1 % maxabs
      scalefactor = max(abs(yref))/max(abs(y));
      ynew = scalefactor*y;
    case 2 % minmax
      scalefactor = (max(yref)-min(yref))/(max(y)-min(y));
      ynew = scalefactor*(y-min(y)) + min(yref);
    case 3 % shift
      nan = isnan(y) | isnan(yref);
      if any(nan)
        error('Cannot use shift scaling with data containing NaN.');
      end
      shift = mean(y)-mean(yref);
      ynew = y - shift;
    case 4 % lsq
      nan = isnan(y) | isnan(yref);
      y_ = y;
      yref_ = yref;
      y_(nan) = [];
      yref_(nan) = [];
      scalefactor = y_\yref_;
      scalefactor = abs(scalefactor);
      ynew = scalefactor*y;
    case 5 % lsq0
      nan = isnan(y) | isnan(yref);
      if any(nan)
        error('Cannot use lsq0 scaling with data containing NaN.');
      end
      D = [y ones(N,1)];
      params = D\yref;
      ynew = D*params;
    case 6 % lsq1
      nan = isnan(y) | isnan(yref);
      if any(nan)
        error('Cannot use lsq1 scaling with data containing NaN.');
      end
      x = (1:N).'/N;
      D = [y ones(N,1) x];
      params = D\yref;
      ynew = D*params;
    case 7 % lsq2
      nan = isnan(y) | isnan(yref);
      if any(nan)
        error('Cannot use lsq2 scaling with data containing NaN.');
      end
      x = (1:N).'/N;
      D = [y ones(N,1) x x.^2];
      params = D\yref;
      ynew = D*params;
  end
  
end

% Preserve row layout
if isRowVector, ynew = ynew.'; end

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
  otherwise
    error('Wrong number of outputs.');
end

return
%===========================================================

%x = linspace(-1,1,1001)*5;
%yref = 3*gaussian(x,0,1);
%y = gaussian(x,0,0.95)/10 + (randn(size(x))-0.5)/150 + x*0.1;
%rescale(y,yref,'lsq1');
