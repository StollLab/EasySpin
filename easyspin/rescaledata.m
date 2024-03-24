% rescaledata    Rescaling of data
%
%   yscaled = rescaledata(y,mode)
%   yscaled = rescaledata(y,yref,mode)
%   yscaled = rescaledata(...,region)
%   [yscaled, scale] = rescaledata(...)
%
% Rescales the data in array y. If given, yref serves
% as the reference. mode determines how the scale factor is calculated.
% Positive scaling is enforced, i.e. the rescaled data is never inverted.
%
% Inputs:
%   y           data to be rescaled, 1D array
%   yref        reference data (used by 'maxabs' and 'lsq'), 1D array
%   mode:
%     'maxabs'  scales y such that maximum magnitude of yscaled is 1 (if
%               yref is not given) or max(abs(yref)) (if yref is given)
%     'lsq'     least-squares fit of a*y to yref; yref is needed
%     'int'     normalize integral (sum of datapoints) to 1
%     'dint'    normalize double integral (sum of cumsum of datapoints) to 1
%     'none'    no scaling
%
% Outputs:
%   yscaled   rescaled data vector
%   scale     scaling factor

function varargout = rescaledata(varargin)

switch nargin
  case 0
    help(mfilename);
    return
  case 1
    y = varargin{1};
    yref = [];
    Mode = 'maxabs';
    region = [];
  case 2
    y = varargin{1};
    input2 = varargin{2};
    if ischar(input2)
      yref = [];
      Mode = input2;
    else
      yref = input2;
      Mode = 'lsq';
    end
    region = [];
  case 3
    y = varargin{1};
    input2 = varargin{2};
    if ischar(input2)
      yref = [];
      Mode = input2;
      region = varargin{3};
    else
      yref = input2;
      Mode = varargin{3};
      region = []; 
    end
  case 4
    y = varargin{1};
    yref = varargin{2};
    Mode = varargin{3};
    region = varargin{4};
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
    error('Second input must be the reference vector.');
  end
end

N = numel(y);

Modes = {'maxabs','lsq','int','dint','none'};
refNeeded = [false true false false false];
equalLengthNeeded = [false true false false false];

ModeID = find(strcmp(Mode,Modes));
if isempty(ModeID)
  Modes_v5 = {'lsq0','lsq1','lsq2','minmax'};
  ModeID = strcmp(Mode,Modes_v5);
  if any(ModeID)
    error('Mode ''%s'' is obsolete. Use basecorr() instead.',Mode);
  else
    error('Unknown scaling mode ''%s''.',Mode);
  end
end
if refNeeded(ModeID)
  if isempty(yref)
    error('For ''%s'' mode, a reference is needed.',Mode);
  end
  if equalLengthNeeded(ModeID)
    equalLength = numel(y)==numel(yref);
    if equalLengthNeeded(ModeID) && ~equalLength
      error('For least-squares rescaling, y and yref must have same number of elements.');
    end
  end
end

y = y(:);
if ~isempty(yref)
  yref = yref(:);
end

% Rescale data
switch ModeID
  case 1  % maxabs
    idx = ~isnan(y);
    if ~isempty(region)
      idx = idx & region(:);
    end
    scalefactor = 1/max(abs(y(idx)));
    if ~isempty(yref)
      idx = ~isnan(yref);
      if ~isempty(region)
        idx = idx & region(:);
      end
      scalefactor = max(abs(yref(idx)))*scalefactor;
    end
  case 2  % lsq
    idx = ~isnan(y) & ~isnan(yref);
    if ~isempty(region)
      idx = idx & region(:);
    end
    % Rescale reference instead of signal (otherwise rmsd(scale) is wrong
    scalefactor = yref(idx)\y(idx);
    scalefactor = 1/scalefactor;
  case 3  % scaling by integral
    idx = ~isnan(y);
    if ~isempty(region)
      idx = idx & region(:);
    end
    scalefactor = 1/sum(y(idx));
  case 4  % scaling by double integral
    idx = ~isnan(y);
    if ~isempty(region)
      idx = idx & region(:);
    end
    scalefactor = 1/sum(cumsum(y(idx)));
  case 5  % no scaling
    scalefactor = 1;
end

% Make sure signal is not inverted
if real(scalefactor)<0
  scalefactor = abs(scalefactor);
end
yscaled = y*scalefactor;

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
    if ~isempty(region)
      masked = ~region(:);
      idx = find([masked(1); diff(masked)]);
      idx = [idx; numel(region)+1];
      ylims = get(gca,'YLim');
      for k = 1:2:numel(idx)-1
        patch([idx(k) idx(k) idx(k+1)-1 idx(k+1)-1],ylims([1 2 2 1]),[0.4902 0.4902 0.4902],'EdgeColor','none','FaceAlpha',0.3)
      end
    end
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
