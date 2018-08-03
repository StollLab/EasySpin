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

N = numel(y);
if (length(y)~=N)
  error('y must be a row or column vector. rescale() does not work on 2D arrays.');
end
isRowVector = size(y,1)==1;

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
      mi = 0; ma = 1;
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
  IdenticalLengthNeeded = [0 0 0 1 1 1 1 1 0];
  if isempty(ModeID)
    error('Unknown scaling mode ''%s''',Mode);
  end
  if IdenticalLengthNeeded(ModeID)
    if numel(y)~=numel(yref)
      error('For least-squares rescaling, vectors must have same number of elements.');
    end
  end
  
  y = y(:);
  yref = yref(:);
  notnan = ~isnan(yref) & ~isnan(y);

  switch ModeID
    case 1 % maxabs
      scalefactor = max(abs(yref(notnan)))/max(abs(y(notnan)));
      ynew = scalefactor*y;
    case 2 % minmax
      scalefactor = (max(yref(notnan))-min(yref(notnan)))/(max(y(notnan))-min(y(notnan)));
      ynew = scalefactor*(y-min(y(notnan))) + min(yref(notnan));
    case 3 % shift
      shift = mean(y(notnan)) - mean(yref(notnan));
      ynew = y - shift;
    case 4 % lsq
      D = y;
      scalefactor = (y(notnan)'*y(notnan))\(y(notnan)'*yref(notnan));
      ynew = D*scalefactor;
    case 5 % lsq0
      D = [y ones(N,1)];
      scalefactor = D(notnan,:)\yref(notnan);
      ynew = D*scalefactor;
    case 6 % lsq1
      x = (1:N).'/N;
      D = [y ones(N,1) x];
      scalefactor = D(notnan,:)\yref(notnan);
      ynew = D*scalefactor;
    case 7 % lsq2
      x = (1:N).'/N;
      D = [y ones(N,1) x x.^2];
      scalefactor = D(notnan,:)\yref(notnan);
      ynew = D*scalefactor;
    case 8 % lsq3
      x = (1:N).'/N;
      D = [y ones(N,1) x x.^2 x.^3];
      scalefactor = D(notnan,:)\yref(notnan);
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
  case 2
    varargout = {ynew scalefactor};
  otherwise
    error('Wrong number of outputs.');
end

return
