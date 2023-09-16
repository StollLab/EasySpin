% basecorr  Baseline correction
%
%   corrdata = basecorr(data,dim,n)
%   corrdata = basecorr(data,dim,n,region)
%   [corrdata,baseline] = basecorr(...)
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
%     corrdata   baseline-corrected data
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
    error('Provide dimension (second input) and polynomial order (third input).');
  case 2
    error('Provide dimension (second input) and polynomial order (third input).');
  case 3
    region = [];
  case 4
    % ok
  otherwise
    error('basecorr() needs 3 inputs: basecorr(data,dimension,order)');
end

if nargout>2, error('Too many output arguments.'); end

if ~isnumeric(data) || isempty(data)
  error('Data must be a non-empty numerical array!');
end
if ~ismatrix(data)
  error('Only 1D or 2D data can be fitted. You rovided %dD data.',ndims(data));
end

if any(n<0) || any(mod(n,1)) || any(n>defaults.maxOrder)
  error('Order must contain integers between 0 and %d!',defaults.maxOrder);
end

twoDimFit = isempty(dim);
if twoDimFit
  if numel(data)==length(data)
    error('Multidimensional fitting not possible for 1D data!');
  end
  if ndims(data)>2  %#ok
    error('Multidimensional fitting for %dD arrays not implemented!',ndims(data));
  end
  if numel(n)~=ndims(data)
    error('n (3rd input) must have %d elements!',ndims(data));
  end
  if any(n>=size(data))
    error('A value in n is too large for the given data!');
  end
else
  if numel(n)~=numel(dim)
    error('dim (2nd input) and n (3rd input) must have the same number of elements.');
  end

  if any(dim>ndims(data)) || any(dim<1)
    error('Dimension out of range!');
  end
end

% Baseline fit
%-------------------------------------------------------------------------------
if twoDimFit
  % Two-dimensional baseline

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
  % 1D baseline corrections along rows or columns

  nPoints = size(data,dim);
  if n>=nPoints
    error('Dimension %d has %d points. Polynomial order %d not possible. Adjust order or dimension.',...
      dim,nPoints,n);
  end

  % Transpose if needed to move dimension of interest along columns
  if dim==2
    data = data.';
  end

  % Build the design matrix
  x = linspace(-1,1,nPoints);  % centered and scaled
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
corrdata = data - baseline;
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
        for k = 1:2:numel(idx)-1
          xregion(idx(k),idx(k+1)-1);
        end
      end
      legend('data','fitted baseline');
      xlabel('index')
      subplot(2,1,2)
      plot(x,corrdata);
      legend('baseline-corrected data');
      xlabel('index')
    else

    end
  case 1
    varargout = {corrdata};
  case 2
    varargout = {corrdata,baseline};
end

end
