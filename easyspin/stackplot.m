% stackplot    Stacked plot of data
%
%    stackplot(x,y)
%    stackplot(x,y,scale)
%    stackplot(x,y,scale,step)
%    stackplot(x,y,scale,step,sliceLabels)
%    stackplot(x,y,scale,step,sliceLabels,colors)
%
%  Plots the columns or rows of y stacked on top of each other. The length
%  of x determines whether columns or rows are plotted. 
%  Slices are rescaled as specified in scale ('maxabs', 'int', 'dint', or 
%  'none', see rescaledata.m for details).
%
% Input: x     - vector of x-axis values
%        y     - matrix of data to plot
%        scale - string defining the type of scaling, options are 'maxabs',
%                'int' (normalized integral), 'dint' (normalized double integral), 
%                or 'none' (default = 'maxabs')
%                an additional scaling factor can be defined by providing a
%                cell as input {scale scalefact}, e.g. {'maxabs',5} (useful if 
%                the shifts in the stacked plot correspond to a specific axis)
%        step  - relative step size for the stacked plot (as a fraction of
%                the maximum amplitude of the plotted data) or list of positions
%                for the different slices on the y-axis
%                (default = 1.5)
%        sliceLabels - list of values or cell array for labeling of the
%                      slices on the y-axis
%        colors      - single color, list of colors or string with colormap name
%                      for plotting of the different slices
%

function varargout = stackplot(x,y,scale,step,sliceLabels,colors)

if nargin==0, help(mfilename); return; end

if nargin<3
  scale = 'maxabs'; 
  scalefact = 1; 
end
if iscell(scale)
  scalefact = scale{2};
  scale = scale{1};
end
if ~ischar(scale)
  error('Third parameter (scale) must be a string (''maxabs'',''int'',''dint'' or ''none'').');
end
if isempty(scale)
  scale = 'maxabs';
end
if ~exist('scalefact','var')
  scalefact = 1;
end

if nargin<4, step = []; end

if nargout>1
  error('stackplot can provide at most one output.');
end

if size(y,2)==numel(x), y = y.'; end

nSlices = size(y,2);
if nargin<5, sliceLabels = []; end
if ~isempty(sliceLabels)
  if numel(sliceLabels)~=nSlices
    error('Number of y tick labels must equal number of slices.');
  end
end

if nargin<6, colors = parula(nSlices+1); colors = colors(1:nSlices,:); end
if ischar(colors)
  colors = eval([colors,'(',num2str(nSlices),')']);
end
if size(colors,1)==1
  colors = repmat(colors,nSlices,1);
elseif size(colors,1)~=nSlices
  error('Number of colors must equal number of slices.');
end

% Rescale data
for k = 1:nSlices
  yplot(:,k) = rescaledata(y(:,k),scale)*scalefact; %#ok
  maxval(k) = max(abs(yplot(:,k))); %#ok
end

% Shift data for stacked plot
if isempty(step)
  step = 1.5;
end
if numel(step)==1
  shift = ((1:nSlices)-1)*step*max(maxval);
elseif numel(step)==nSlices
  shift = step;
else
  error(['Fourth parameter (step) must be a single number or a vector with' ...
         'a number of values equal to the number of slices.'])
end
if shift(end)<shift(1)
  shift = shift(end:-1:1);
end
yplot = yplot + shift;

if min(size(x))==1
  x = x(:);
  x = repmat(x,1,size(y,2));
end

% Plot data
h = plot(x,yplot);
set(gca,'YTick',shift);
set(gca,'YTickLabel',sliceLabels);
colororder(colors)
axis tight

if nargout==1
  varargout{1} = h;
end
