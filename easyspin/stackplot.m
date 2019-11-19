% stackplot    Stacked plot of data
%
%    stackplot(x,y)
%    stackplot(x,y,scale)
%    stackplot(x,y,scale,step)
%
%  Plots the columns or rows of y stacked on top of each other. The length
%  of x determines whether columns or rows are plotted. Slices are rescaled
%  to max-min=scale and plotted at distance of step.
%
%  If step is missing, step = 1 is assumed.
%  If scale is missing, scale = 1 is assumed.
%  If scale is set to 0, no rescaling is done.

function varargout = stackplot(x,y,scale,step,SliceLabels)

if nargin==0, help(mfilename); return; end

if nargin<3, scale = 1; end
if nargin<4, step = 1; end

if nargout>1
  error('stackplot can provide at most one output.');
end

if ischar(scale) || numel(scale)>1
  error('Third parameter (scale) must be a number.');
end

if size(y,2)==numel(x), y = y.'; end

nSlices = size(y,2);
if nargin<5, SliceLabels = 1:nSlices; end
if numel(SliceLabels)~=nSlices
  error('Number of y tick labels must equal number of slices.');
end
for k = 1:nSlices
  yy = y(:,k);
  shift(k) = (k-1)*step;
  if scale<0
    y(:,k) = yy/max(abs((yy)))*abs(scale) + shift(k);
  elseif scale>0
    y(:,k) = (yy-min(yy))/(max(yy)-min(yy))*scale + shift(k);
  else
    y(:,k) = yy + shift(k);
  end
end

if min(size(x))==1
  x = x(:);
  x = repmat(x,1,size(y,2));
end

if shift(end)<shift(1)
  shift = shift(end:-1:1);
end

h = plot(x,y);
set(gca,'YTick',shift);
set(gca,'YTickLabel',SliceLabels);
axis tight

if nargout==1
  varargout={h};
end
