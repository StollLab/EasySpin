% basecorr  Baseline correction 
%
%   CorrSpectrum = basecorr(Spectrum,Dimension,Order)
%   [CorrSpectrum,BaseLine] = basecorr(...)
%
%   Makes a polynomial baseline correction for the array Spectrum.
%   Dimension gives a list of dimensions to which successive 1D
%   corrections should be applied. If Dimensions=[], then a multi-
%   dimensional (hyper)surface fitting is performed along all
%   dimensions.
%   Order gives a list of the polynomial orders for all dimensions
%   specified. All orders must be between 0 and 4.
%
%   Examples:
%      z = basecorr(z,1,2); % 2nd order along dimension 1
%      z = basecorr(z,[],3); % twodimensional 3rd order surface
%      z = basecorr(z,[1 2],[3 1]); % 3rd order along dim 1 and 1st
%                                   % order along dim 2.


function varargout = basecorr(Spectrum,Dimension,Order)

Defaults.maxOrder = 6;

if (nargin==0), help(mfilename); return; end

switch nargin
  case 0, help(mfilename); return;
  case 1, Dimension = 1; if size(Spectrum,1)==1, Dimension = 2; end; Order = 1;
  case 2, Order = 1;
  case 3,
  otherwise
    error('basecorr() needs 3 inputs: basecorr(Spectrum,Dimension,Order)');
end

if (nargout<1), error('Not enough output arguments.'); end
if (nargout>2), error('Too many output arguments.'); end

if ~isnumeric(Spectrum) || isempty(Spectrum),
  error('Spectrum must be a non-empty numerical array!');
end

if any(Order<0) || any(mod(Order,1)) || any(Order>Defaults.maxOrder)
  error('Order must contain integers between 0 and %d!',Defaults.maxOrder);
end

if isempty(Dimension)
  
  % Multidimensional fitting
  %-----------------------------------------------------------------------
  
  if numel(Spectrum)==length(Spectrum)
    error('Multidimensional fitting not possible for 1D data!');
  end
  if ndims(Spectrum)>2
    error('Multidimensional fitting for %dD arrays not implemented!',ndims(Spectrum));
  end
  if numel(Order)~=ndims(Spectrum)
    error('Order must have %d elements!',ndims(Spectrum));
  end
  if any(Order>=size(Spectrum))
    error('A value in Order is too large for the given Spectrum!');
  end
  [m,n] = size(Spectrum);
  x = repmat(1:m,n,1).'; x = x(:);
  y = repmat(1:n,m,1); y = y(:);
  % Build the design matrix
  q = 0;
  for j = 0:Order(2),  % each column a x^i*y^j monomial vector
    for i = 0:Order(1),
      q = q+1;
      D(:,q) = x.^i .* y.^j;
    end
  end
  % Solve least squares problem
  C = D\Spectrum(:);
  % Compute baseline
  BaseLine = zeros(size(x));
  for q = 1:size(D,2)
    BaseLine = BaseLine + C(q)*D(:,q);
  end
  BaseLine = reshape(BaseLine,m,n);
  CorrSpectrum = Spectrum - BaseLine;
  
else
  
  % Multiple 1D baseline corrections
  %-------------------------------------------------------------------------------
  
  if numel(Order)~=numel(Dimension)
    error('Order and Dimension must have the same number of elements!');
  end
  
  if any(Dimension>ndims(Spectrum)) || any(Dimension<1)
    error('Dimension out of range!');
  end
  
  BaseLine = zeros(size(Spectrum));
  CorrSpectrum = Spectrum;
  for q = 1:numel(Dimension)
    d = Dimension(q);
    if Order(q)>=size(Spectrum,d),
      error('Dimension %d has %d points. Polynomial order %d not possible. Check Order and Dimension!',...
        d,size(Spectrum,d),Order(q));
    end
    thisSpectrum = permute(CorrSpectrum,[d 1:d-1 d+1:ndims(Spectrum)]);
    x = (1:size(thisSpectrum,1)).';
    thisBaseLine = zeros(size(thisSpectrum));
    if Order(q)==0
      p = mean(thisSpectrum);
      for c=size(thisSpectrum,2):-1:1
        thisBaseLine(:,c) = p(c);
      end
    else
      [p,dummy,mu] = polyfitvec(x,thisSpectrum,Order(q));
      for c = size(thisSpectrum,2):-1:1
        %[p,dummy,mu] = polyfit(x,thisSpectrum(:,c),Order(q));
        %thisBaseLine(:,c) = thisBaseLine(:,c) + polyval(p,(x-mu(1))/mu(2));
        thisBaseLine(:,c) = thisBaseLine(:,c) + polyval(p(c,:),(x-mu(1))/mu(2));
      end
    end
    thisBaseLine = permute(thisBaseLine,[2:d 1 d+1:ndims(Spectrum)]);
    BaseLine = BaseLine + thisBaseLine;
    CorrSpectrum = CorrSpectrum - thisBaseLine;
  end
  
end

% Prepare output.
varargout = {CorrSpectrum, BaseLine};
varargout = varargout(1:nargout);

return
