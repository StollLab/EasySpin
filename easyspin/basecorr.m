% basecorr  Baseline correction
%
%   datacorr = basecorr(data,dim,n)
%   datacorr = basecorr(data,dim,n,region)
%   [datacorr,baseline] = basecorr(...)
%
%   Performs polynomial baseline correction.
%
%   Inputs:
%     data       data to baseline correct, 1D or 2D
%     dim        dimension along which to baseline correct; if [], then
%                  a 2D function is fitted.
%     n          polynomial order, scalar for 1D data or row/column-wise 2D data,
%                and 2-element array for 2D; between 0 and 6
%     region     logic array, providing a mask for the region to include in
%                the baseline fit; true for each point to be included, and
%                false for each point to be excluded
%
%   Output:
%     datacorr   baseline-corrected data
%     baseline   fitted baseline
%
%   Examples:
%      z = basecorr(z,1,2);         % 2nd-order along dimension 1
%
%      region = x<10 | x>50;        % include x regions <10 and >50
%      z = basecorr(z,1,2,region);
%
%      z = basecorr(z,3,[]);        % two-dimensional 3rd-order surface
%

function varargout = basecorr(data,dim,n,region)

defaults.maxOrder = 6;

switch nargin
  case 0
    help(mfilename);
    return
  case 1
    error('Provide dimension dim (second input) and polynomial order n (third input).');
  case 2
    error('Provide polynomial order n (third input).');
  case 3
    region = [];
  case 4
    % ok
  otherwise
    error('basecorr() needs at least 3 inputs: basecorr(data,dim,n)');
end

if nargout>2, error('Too many output arguments.'); end


% Input argument checks
%-------------------------------------------------------------------------------
% Check data (first input)
if ~isnumeric(data) || isempty(data)
  error('Data must be a non-empty numerical array.');
end
if ~ismatrix(data)
  error('Only 1D or 2D data can be fitted. You provided a %dD data array.',ndims(data));
end

% Check dimension (second input)
if ~isempty(dim) && numel(dim)~=1 && dim~=1 && dim~=2
  error('dim (2nd input) must either 1 (fit columns), 2 (fit rows), or [] (2D fit).');
end
twoDimFit = isempty(dim);

% Convert row to column vector for 1D fit
rowVector = false;
if ~twoDimFit && dim==1 && size(data,1)==1
  rowVector = true;
  data = data(:);
end

% Check polyomial order(s) (third input)
if ~isnumeric(n) || any(n<0) || any(mod(n,1)) || any(n>defaults.maxOrder)
  error('Order must contain integers between 0 and %d!',defaults.maxOrder);
end

if twoDimFit
  if numel(n)~=2
    error('For 2D fit, the polynomial order n (3rd input) must have 2 elements, one for each dimension.');
  end
  for iDim = 1:2
    if n(iDim)>=size(data,iDim)
      error('The provided polynomial order n=%d must be smaller than the number of points, %d.',n,size(data,iDim));
    end
  end
else
  if numel(n)~=1
    error('Polynomial order (3rd input) must be a single number.');
  end
  if n>=size(data,dim)
    error('The provided polynomial order n=%d must be smaller than the number of points, %d.',n,size(data,dim));
  end
end

% Check region (fourth input)
if ~isempty(region)
  if ~islogical(region)
    error('region (4th input) must be a logical array.');
  end
  if twoDimFit
    warning('basecorr: Baseline correction with mask is only available for 1D fits.')
  end
end


% Baseline fit
%-------------------------------------------------------------------------------
if twoDimFit
  % Two-dimensional fit

  % Build the design matrix
  [r,c] = size(data);
  x = reshape(repmat(linspace(-1,1,r),c,1).',[],1);  % centered and scaled
  y = reshape(repmat(linspace(-1,1,c),r,1),[],1);  % centered and scaled
  col = 0;
  D = zeros(r*c,prod(n+1));
  for j = 0:n(2)
    for i = 0:n(1)
      col = col + 1;
      D(:,col) = x.^i .* y.^j;  % each column of D is a x^i*y^j monomial vector
    end
  end

  % Fit polynomial and evaluate baseline
  p = D\data(:);
  baseline = D*p;
  baseline = reshape(baseline,r,c);

else
  % One-dimensional fits along columns (dim=1) or rows (dim=2)

  % Transpose if needed to move dimension of interest along columns
  if dim==2
    data = data.';
  end

  % Build the design matrix
  x = linspace(-1,1,size(data,1));  % centered and scaled
  D = x(:).^(0:n);

  % Fit polynomials and evaluate baseline
  if isempty(region)
    p = D\data;
  else
    p = D(region(:),:)\data(region(:),:);
  end
  baseline = D*p;

  % Undo transposition
  if dim==2
    baseline = baseline.';
    data = data.';
  end

end


% Output, plotting
%-------------------------------------------------------------------------------
datacorr = data - baseline;
if rowVector
  datacorr = datacorr.';
  baseline = baseline.';
end
switch nargout
  case 0
    if isvector(data)
      x = 1:length(data);
      subplot(2,1,1)
      plot(x,data,x,baseline);
      if ~isempty(region)
        masked = ~region(:);
        idx = find([masked(1); diff(masked)]);
        idx = [idx; numel(region)+1];
        ylims = get(gca,'YLim');
        for k = 1:2:numel(idx)-1
          patch([idx(k) idx(k) idx(k+1)-1 idx(k+1)-1],ylims([1 2 2 1]),[0.4902 0.4902 0.4902],'EdgeColor','none','FaceAlpha',0.3)
        end
      end
      legend('data','fitted baseline');
      axis tight
      xlabel('index')
      subplot(2,1,2)
      plot(x,datacorr);
      yline(0)
      legend('baseline-corrected data');
      axis tight
      xlabel('index')
    else

    end
  case 1
    varargout = {datacorr};
  case 2
    varargout = {datacorr,baseline};
end

end
