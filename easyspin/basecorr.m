% basecorr  Baseline correction
%
%   corrdata = basecorr(data,order)
%   corrdata = basecorr(data,order,dimension)
%   corrdata = basecorr(data,order,dimension,region)
%   [corrdata,baseline] = basecorr(...)
%
%   Makes a polynomial baseline correction.
%
%   Inputs:
%     data       data to baseline correct, 1D or 2D
%     order      polynomial orders, scalar for 1D data, and 2-element array for 2D
%                all orders must be between 0 (offset) and 6
%     dimension  dimension(s) along which to baseline correct; if [], then
%                  a 2D function is fitted.
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
%      z = basecorr(z,3,[]);        % two-dimensional 3rd order surface
%      z = basecorr(z,[3 1],[1 2]); % 3rd-order along dim 1 and 1st
%                                   % order along dim 2
%
%      region = x<10 | x>50;        % include x regions <10 and >50
%      z = basecorr(z,1,2,region);

function varargout = basecorr(data,order,dimension,region)

defaults.maxOrder = 6;

switch nargin
  case 0
    help(mfilename);
    return
  case 1
    error('Provide the polynomial order as second input.');
  case 2
    if iscolumn(data)
      dimension = 1;
    elseif isrow(data)
      dimension = 2;
    else
      error('For 2D data, provide dimension(s) as third input.')
    end
    region = [];
  case 3
    region = [];
  case 4
  otherwise
    error('basecorr() needs 3 inputs: basecorr(Spectrum,Dimension,Order)');
end

if nargout>2, error('Too many output arguments.'); end

if ~isnumeric(data) || isempty(data)
  error('Spectrum must be a non-empty numerical array!');
end

if any(order<0) || any(mod(order,1)) || any(order>defaults.maxOrder)
  error('Order must contain integers between 0 and %d!',defaults.maxOrder);
end

twoDimFit = isempty(dimension);
if twoDimFit
  if numel(data)==length(data)
    error('Multidimensional fitting not possible for 1D data!');
  end
  if ndims(data)>2  %#ok
    error('Multidimensional fitting for %dD arrays not implemented!',ndims(data));
  end
  if numel(order)~=ndims(data)
    error('Order must have %d elements!',ndims(data));
  end
  if any(order>=size(data))
    error('A value in Order is too large for the given Spectrum!');
  end
else
  if numel(order)~=numel(dimension)
    error('order (2nd input) and dimension (3rd input) must have the same number of elements.');
  end

  if any(dimension>ndims(data)) || any(dimension<1)
    error('Dimension out of range!');
  end
end

% Baseline fit
%-------------------------------------------------------------------------------
if twoDimFit
  % Two-dimensional baseline

  % Build the design matrix
  [m,n] = size(data);
  x = reshape(repmat(1:m,n,1).',[],1);
  y = reshape(repmat(1:n,m,1),[],1);
  col = 0;
  D = zeros(m*n,prod(order+1));
  for j = 0:order(2) 
    for i = 0:order(1)
      col = col + 1;
      D(:,col) = x.^i .* y.^j;  % each column of D is a x^i*y^j monomial vector
    end
  end

  % Solve least-squares problem
  C = D\data(:);

  % Compute baseline
  baseline = reshape(D*C,m,n);

  % Correct data
  corrdata = data - baseline;

else
  % Multiple 1D baseline corrections


  baseline = zeros(size(data));
  corrdata = data;
  for col = 1:numel(dimension)
    d = dimension(col);
    if order(col)>=size(data,d)
      error('Dimension %d has %d points. Polynomial order %d not possible. Adjust order or dimension.',...
        d,size(data,d),order(col));
    end

    % Permute to move dimension of interest to first dimension
    spectrum_ = permute(corrdata,[d 1:d-1 d+1:ndims(data)]);
    if order(col)==0
      if ~isempty(region)
        p = mean(spectrum_(region));
      else
        p = mean(spectrum_);
      end
      baseline_ = zeros(size(spectrum_));
      for c = 1:size(spectrum_,2)
        baseline_(:,c) = p(c);
      end
    else
      x = (1:size(spectrum_,1)).';
      if ~isempty(region)
        [p,~,mu] = polyfitvec(x(region),spectrum_(region),order(col));
      else
        [p,~,mu] = polyfitvec(x,spectrum_,order(col));
      end
      baseline_ = zeros(size(spectrum_));
      for c = size(spectrum_,2):-1:1
        baseline_(:,c) = baseline_(:,c) + polyval(p(c,:),(x-mu(1))/mu(2));
      end
    end
    % Undo opermutation
    baseline_ = permute(baseline_,[2:d 1 d+1:ndims(data)]);

    % Add baselines, correct data
    baseline = baseline + baseline_;
    corrdata = corrdata - baseline_;
  end

end

% Output, plotting
%-------------------------------------------------------------------------------
switch nargout
  case 0
    if isvector(data)
      x = 1:length(data);
      subplot(2,1,1)
      plot(x,data,x,baseline);
      if ~isempty(region)
        idx = find([1 diff(region)]);
        idx = [idx numel(region)+1];
        for k = 1:2:numel(idx)-1
          xregion(idx(k),idx(k+1)-1,'FaceColor',[0.7 1 0.8]);
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
